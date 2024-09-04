import logging
import typing

import hpotk

from phenopackets.schema.v2.phenopackets_pb2 import Phenopacket
import phenopackets.schema.v2.core.individual_pb2 as ppi
from phenopackets.schema.v2.core.disease_pb2 import Disease as PPDisease
from phenopackets.schema.v2.core.interpretation_pb2 import GenomicInterpretation
from phenopackets.vrsatile.v1.vrsatile_pb2 import VcfRecord, VariationDescriptor
from phenopackets.vrs.v1.vrs_pb2 import Variation

from gpsea.model import SampleLabels, Patient, Sex, Disease
from gpsea.model import (
    VariantClass,
    VariantCoordinates,
    ImpreciseSvInfo,
    VariantInfo,
    Variant,
    Genotype,
    Genotypes,
)
from gpsea.model.genome import GenomeBuild, GenomicRegion, Strand
from ._api import (
    VariantCoordinateFinder,
    FunctionalAnnotator,
    ImpreciseSvFunctionalAnnotator,
)
from ._audit import Notepad
from ._patient import PatientCreator
from ._phenotype import PhenotypeCreator


class PhenopacketGenotypeParser:
    """
    `PhenopacketGenotypeParser` tries to extract :class:`Genotype` from `GenomicInterpretation`.
    """

    def find_genotype(
        self,
        item: GenomicInterpretation,
    ) -> typing.Optional[Genotype]:
        if item.HasField("variant_interpretation"):
            variant_interpretation = item.variant_interpretation
            if variant_interpretation.HasField("variation_descriptor"):
                variation_descriptor = variant_interpretation.variation_descriptor
                if variation_descriptor.HasField("allelic_state"):
                    genotype = variation_descriptor.allelic_state.label
                    return self._map_geno_genotype_label(genotype)

        return None

    @staticmethod
    def _map_geno_genotype_label(genotype: str) -> Genotype:
        """
        Mapping from labels of the relevant GENO terms that is valid as of Oct 2nd, 2023.
        """
        if genotype in ("heterozygous", "compound heterozygous", "simple heterozygous"):
            return Genotype.HETEROZYGOUS
        elif genotype == "homozygous":
            return Genotype.HOMOZYGOUS_ALTERNATE
        elif genotype in ("hemizygous", "hemizygous X-linked", "hemizygous Y-linked"):
            return Genotype.HEMIZYGOUS
        else:
            raise ValueError(f"Unknown genotype {genotype}")


class PhenopacketVariantCoordinateFinder(
    VariantCoordinateFinder[GenomicInterpretation]
):
    """
    `PhenopacketVariantCoordinateFinder` figures out :class:`~gpsea.model.VariantCoordinates`
    and :class:`~gpsea.model.Genotype` from `GenomicInterpretation` element of Phenopacket Schema.

    :param build: genome build to use in `VariantCoordinates`
    :param hgvs_coordinate_finder: the coordinate finder to use for parsing HGVS expressions
    """

    def __init__(
        self,
        build: GenomeBuild,
        hgvs_coordinate_finder: VariantCoordinateFinder[str],
    ):
        self._logger = logging.getLogger(__name__)
        self._build = hpotk.util.validate_instance(build, GenomeBuild, "build")
        self._hgvs_finder = hpotk.util.validate_instance(
            hgvs_coordinate_finder,
            VariantCoordinateFinder,
            "hgvs_coordinate_finder",
        )

    def find_coordinates(
        self,
        item: GenomicInterpretation,
    ) -> typing.Optional[VariantCoordinates]:
        """
        Tries to extract the variant coordinates from the `GenomicInterpretation`.

        Args:
            item (GenomicInterpretation): a genomic interpretation element from Phenopacket Schema

        Returns:
            typing.Optional[VariantCoordinates]: variant coordinates
        """
        if not isinstance(item, GenomicInterpretation):
            raise ValueError(
                f"item must be a Phenopacket GenomicInterpretation but was type {type(item)}"
            )

        variation_descriptor = item.variant_interpretation.variation_descriptor

        if self._vcf_is_available(variation_descriptor.vcf_record):
            # We have a VCF record.
            if not self._check_assembly(
                variation_descriptor.vcf_record.genome_assembly
            ):
                raise ValueError(
                    f"Variant id {variation_descriptor.id} for patient {item.subject_or_biosample_id} "
                    "has a different Genome Assembly than what was given. "
                    f"{variation_descriptor.vcf_record.genome_assembly} is not {self._build.identifier}."
                )
            contig = self._build.contig_by_name(variation_descriptor.vcf_record.chrom)
            start = int(variation_descriptor.vcf_record.pos) - 1
            ref = variation_descriptor.vcf_record.ref
            alt = variation_descriptor.vcf_record.alt
            end = start + len(ref)
            change_length = len(alt) - len(ref)

            region = GenomicRegion(contig, start, end, Strand.POSITIVE)
            return VariantCoordinates(region, ref, alt, change_length)
        elif self._cnv_is_available(variation_descriptor.variation):
            # We have a CNV.
            variation = variation_descriptor.variation
            seq_location = variation.copy_number.allele.sequence_location
            refseq_contig_name = seq_location.sequence_id.split(":")[1]
            contig = self._build.contig_by_name(refseq_contig_name)

            # Assuming SV coordinates are 1-based (VCF style),
            # so we subtract 1 to transform to 0-based coordinate system
            start = int(seq_location.sequence_interval.start_number.value) - 1
            end = int(seq_location.sequence_interval.end_number.value)
            ref = "N"
            number = int(variation.copy_number.number.value)
            if number == 1:
                alt = "<DEL>"
                change_length = -(end - start)
            elif number == 3:
                alt = "<DUP>"
                change_length = end - start
            else:
                raise ValueError(
                    f"The copy number of {number} is not supported. Supported values: {{1, 3}}"
                )

            region = GenomicRegion(contig, start, end, Strand.POSITIVE)
            return VariantCoordinates(region, ref, alt, change_length)
        elif len(variation_descriptor.expressions) > 0:
            # We have some expressions. Let's try to find the 1st expression with `hgvs.c` syntax.
            for expression in variation_descriptor.expressions:
                if expression.syntax == "hgvs.c":
                    return self._hgvs_finder.find_coordinates(expression.value)
        elif self._looks_like_large_sv(variation_descriptor):
            # We cannot extract exact variant coordinates from a variation descriptor in this format.
            return None
        else:
            raise ValueError("Unable to find variant coordinates.")

    def _check_assembly(self, genome_assembly: str) -> bool:
        if "38" in genome_assembly and self._build.identifier == "GRCh38.p13":
            return True
        elif (
            "37" in genome_assembly or "19" in genome_assembly
        ) and self._build.identifier == "GRCh37.p13":
            return True
        else:
            return False

    @staticmethod
    def _vcf_is_available(vcf_record: VcfRecord) -> bool:
        """
        Check if we can parse data out of VCF record.
        """
        return (
            vcf_record.genome_assembly != ""
            and vcf_record.chrom != ""
            and vcf_record.pos >= 0
            and vcf_record.ref != ""
            and vcf_record.alt != ""
        )

    @staticmethod
    def _cnv_is_available(variation: Variation):
        seq_location = variation.copy_number.allele.sequence_location
        return (
            seq_location.sequence_id != ""
            and seq_location.sequence_interval.start_number.value >= 0
            and seq_location.sequence_interval.end_number.value >= 0
            and variation.copy_number.number.value != ""
        )

    @staticmethod
    def _looks_like_large_sv(
        variation_descriptor: VariationDescriptor,
    ) -> bool:
        structural_type = (
            variation_descriptor.structural_type
            if variation_descriptor.HasField("structural_type")
            else None
        )
        gene_context = (
            variation_descriptor.gene_context
            if variation_descriptor.HasField("gene_context")
            else None
        )

        # If we have these fields, we seem to have all information
        # to parse the variation descriptor elsewhere.
        return structural_type is not None and gene_context is not None


class PhenopacketPatientCreator(PatientCreator[Phenopacket]):
    """
    `PhenopacketPatientCreator` transforms `Phenopacket` into :class:`~gpsea.model.Patient`.
    """

    def __init__(
        self,
        build: GenomeBuild,
        phenotype_creator: PhenotypeCreator,
        functional_annotator: FunctionalAnnotator,
        imprecise_sv_functional_annotator: ImpreciseSvFunctionalAnnotator,
        hgvs_coordinate_finder: VariantCoordinateFinder[str],
        assume_karyotypic_sex: bool = True,
    ):
        """
        Create an instance.

        :param build: the genome build to use for variants.
        :param phenotype_creator: a phenotype creator for creating phenotypes.
        :param functional_annotator: for computing functional annotations.
        :param imprecise_sv_functional_annotator: for getting info about imprecise variants.
        :param hgvs_coordinate_finder: for finding chromosomal coordinates for HGVS variant descriptions.
        :param assume_karyotypic_sex: `True` if it is OK to assume that `FEMALE` has the `XX` karyotype
            and `MALE` has `XY`.
        """
        self._logger = logging.getLogger(__name__)
        # Violates DI, but it is specific to this class, so I'll leave it "as is".
        self._coord_finder = PhenopacketVariantCoordinateFinder(
            build, hgvs_coordinate_finder
        )
        self._gt_parser = PhenopacketGenotypeParser()
        self._phenotype_creator = hpotk.util.validate_instance(
            phenotype_creator, PhenotypeCreator, "phenotype_creator"
        )
        self._functional_annotator = hpotk.util.validate_instance(
            functional_annotator, FunctionalAnnotator, "functional_annotator"
        )
        self._imprecise_sv_functional_annotator = hpotk.util.validate_instance(
            imprecise_sv_functional_annotator,
            ImpreciseSvFunctionalAnnotator,
            "imprecise_sv_functional_annotator",
        )

        # Set of sequence ontology IDs that we will treat as a deletion (`DEL`)
        # for the purpose of assigning imprecise SV info with a variant class.
        self._so_deletions = {
            "1000029",  # chromosomal deletion: An incomplete chromosome.
            # transcript ablation: A feature ablation whereby the deleted region includes a transcript feature.
            "0001893",
            # feature_ablation: A sequence variant, caused by an alteration of the genomic sequence,
            # where the deletion, is greater than the extent of the underlying genomic features.
            "0001879",
        }
        self._so_duplications = {
            "1000037",  # chromosomal_duplication
        }
        self._assume_karyotypic_sex = assume_karyotypic_sex

    def process(self, pp: Phenopacket, notepad: Notepad) -> Patient:
        """Creates a Patient from the data in a given Phenopacket

        Args:
            pp (Phenopacket): A Phenopacket object
            notepad (Notepad): notepad to write down the issues
        Returns:
            Patient: A Patient object
        """
        sample_id = SampleLabels(
            label=pp.subject.id,
            meta_label=pp.id if len(pp.id) > 0 else None,
        )

        # Extract karyotypic sex
        indi = notepad.add_subsection("individual")
        sex = self._extract_sex(pp.subject, indi)

        # Check phenotypes
        pfs = notepad.add_subsection("phenotype-features")
        phenotypes = self._phenotype_creator.process(
            ((pf.type.id, not pf.excluded) for pf in pp.phenotypic_features), pfs
        )

        # Check diseases
        dip = notepad.add_subsection("diseases")
        diseases = self._add_diseases(pp.diseases, dip)

        # Check variants
        vs = notepad.add_subsection("variants")
        variants = self._add_variants(sample_id, pp, vs)

        return Patient.from_raw_parts(
            sample_id,
            sex=sex,
            phenotypes=phenotypes,
            variants=variants,
            diseases=diseases,
        )

    def _add_diseases(
        self, diseases: typing.Sequence[PPDisease], notepad: Notepad
    ) -> typing.Sequence[Disease]:
        """Creates a list of Disease objects from the data in a given Phenopacket

        Args:
            diseases (Sequence[Disease]): A sequence of Phenopacket Disease objects
            notepad (Notepad): notepad to write down the issues
        Returns:
            Sequence[Dis]: A list of Disease objects
        """
        if len(diseases) == 0:
            notepad.add_warning("No diseases found")
            return ()

        final_diseases = []
        for i, dis in enumerate(diseases):
            if not dis.HasField("term"):
                notepad.add_error(f"#{i} has no `term`")
                continue
            else:
                term_id = hpotk.TermId.from_curie(dis.term.id)

            # Do not include excluded diseases if we decide to assume excluded if not included
            final_diseases.append(Disease(term_id, dis.term.label, not dis.excluded))

        return final_diseases

    def _extract_sex(
        self,
        individual: ppi.Individual,
        notepad: Notepad,
    ) -> typing.Optional[Sex]:
        # Let's use the phenotypic sex as fallback
        sex = individual.sex
        if sex == ppi.FEMALE:
            return Sex.FEMALE
        elif sex == ppi.MALE:
            return Sex.MALE
        elif sex == ppi.OTHER_SEX or sex == ppi.UNKNOWN_SEX:
            return Sex.UNKNOWN_SEX
        else:
            notepad.add_warning(f'Unknown sex type: {sex}')
            return Sex.UNKNOWN_SEX

    def _add_variants(
        self,
        sample_id: SampleLabels,
        pp: Phenopacket,
        notepad: Notepad,
    ) -> typing.Sequence[Variant]:
        """Creates a list of Variant objects from the data in a given Phenopacket

        Args:
            pp (Phenopacket): A Phenopacket object
            notepad (Notepad): notepad to write down the issues
        Returns:
            Sequence[Variant]: A list of Variant objects
        """
        variants = []

        for i, interpretation in enumerate(pp.interpretations):
            sub_note = notepad.add_subsection(f"#{i}")
            if interpretation.HasField("diagnosis"):
                for (
                    genomic_interpretation
                ) in interpretation.diagnosis.genomic_interpretations:
                    gt = self._gt_parser.find_genotype(genomic_interpretation)
                    if gt is None:
                        sub_note.add_warning(
                            "Could not extract genotype from genomic interpretation",
                            "Remove variant from testing",
                        )
                        continue

                    variant_info = self._extract_variant_info(
                        sample_id, genomic_interpretation, sub_note
                    )
                    if variant_info is None:
                        # We already complained in the extract function
                        continue

                    if variant_info.has_variant_coordinates():
                        try:
                            tx_annotations = self._functional_annotator.annotate(
                                variant_coordinates=variant_info.variant_coordinates,
                            )
                        except ValueError as error:
                            sub_note.add_warning(
                                f"Patient {pp.id} has an error with variant {variant_info.variant_key}",
                                f"Try again or remove variant form testing... {error}",
                            )
                            continue
                    elif variant_info.has_sv_info():
                        try:
                            tx_annotations = (
                                self._imprecise_sv_functional_annotator.annotate(
                                    item=variant_info.sv_info,
                                )
                            )
                        except ValueError as error:
                            sub_note.add_warning(
                                f"Patient {pp.id} has an error with variant {variant_info.variant_key}",
                                f"Try again or remove variant form testing... {error}",
                            )
                            continue
                    else:
                        raise ValueError(
                            "VariantInfo should have either the coordinates or the SV info, but had neither!"
                        )

                    if len(tx_annotations) == 0:
                        sub_note.add_warning(
                            f"Patient {pp.id} has an error with variant {variant_info.variant_key}",
                            "Remove variant from testing... tx_anno == 0",
                        )
                        continue

                    genotype = Genotypes.single(sample_id, gt)
                    variants.append(
                        Variant(
                            variant_info=variant_info,
                            tx_annotations=tx_annotations,
                            genotypes=genotype,
                        )
                    )

        if len(variants) == 0:
            notepad.add_warning(f"Patient {pp.id} has no variants to work with")

        return variants

    def _extract_variant_info(
        self,
        sample_id: SampleLabels,
        genomic_interpretation: GenomicInterpretation,
        notepad: Notepad,
    ) -> typing.Optional[VariantInfo]:
        variant_coordinates = None
        sv_info = None

        try:
            variant_coordinates = self._coord_finder.find_coordinates(
                genomic_interpretation
            )
        except ValueError:
            notepad.add_warning(
                "Expected a VCF record, a VRS CNV, or an expression with `hgvs.c`"
                f"but had an error retrieving any from patient {sample_id}",
                "Remove variant from testing",
            )
            return None

        if variant_coordinates is None:
            sv_info = self._map_to_imprecise_sv(genomic_interpretation)
            if sv_info is None:
                notepad.add_warning(
                    "Could not extract the information for large SV annotation",
                    "Remove variant from testing",
                )
                return None

        return VariantInfo(
            variant_coordinates=variant_coordinates,
            sv_info=sv_info,
        )

    def _map_to_imprecise_sv(
        self,
        genomic_interpretation: GenomicInterpretation,
    ) -> typing.Optional[ImpreciseSvInfo]:
        if genomic_interpretation.HasField("variant_interpretation"):
            variant_interpretation = genomic_interpretation.variant_interpretation
            if variant_interpretation.HasField("variation_descriptor"):
                variation_descriptor = variant_interpretation.variation_descriptor

                structural_type = (
                    variation_descriptor.structural_type
                    if variation_descriptor.HasField("structural_type")
                    else None
                )
                gene_context = (
                    variation_descriptor.gene_context
                    if variation_descriptor.HasField("gene_context")
                    else None
                )

                if structural_type is not None and gene_context is not None:
                    st = hpotk.TermId.from_curie(curie=structural_type.id)
                    variant_class = self._map_structural_type_to_variant_class(st)
                    return ImpreciseSvInfo(
                        structural_type=st,
                        variant_class=variant_class,
                        gene_id=gene_context.value_id,
                        gene_symbol=gene_context.symbol,
                    )

        return None

    def _map_structural_type_to_variant_class(
        self,
        structural_type: hpotk.TermId,
    ) -> VariantClass:
        # This method is most likely incomplete.
        # Please open a ticket if you receive a `ValueError`
        # for a structural type, that is not mapped at the moment,
        # to help us enhance the mapping.
        if structural_type.prefix == "SO":
            if structural_type.id in self._so_deletions:
                return VariantClass.DEL
            elif structural_type.id in self._so_duplications:
                return VariantClass.DUP
            else:
                raise ValueError(f"Unknown structural type {structural_type}")
        else:
            raise ValueError(f"Unknown structural type {structural_type}")
