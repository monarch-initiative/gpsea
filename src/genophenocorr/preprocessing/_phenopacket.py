import logging
import os
import typing

import hpotk

# pyright: reportGeneralTypeIssues=false
from google.protobuf.json_format import Parse
from phenopackets import GenomicInterpretation, Phenopacket
from tqdm import tqdm

from genophenocorr.model import Patient, Cohort, SampleLabels
from genophenocorr.model import Phenotype, ProteinMetadata, VariantCoordinates, Variant, Genotype, Genotypes
from genophenocorr.model.genome import GenomeBuild, GenomicRegion, Strand
from ._api import VariantCoordinateFinder, FunctionalAnnotator
from ._audit import AuditReport
from ._patient import PatientCreator
from ._phenotype import PhenotypeCreator


class PhenopacketVariantCoordinateFinder(VariantCoordinateFinder[GenomicInterpretation]):
    """
    `PhenopacketVariantCoordinateFinder` figures out :class:`genophenocorr.model.VariantCoordinates`
    and :class:`genophenocorr.model.Genotype` from `GenomicInterpretation` element of Phenopacket Schema.

    :param build: genome build to use in `VariantCoordinates
    :param hgvs_coordinate_finder: the coordinate finder to use for parsing HGVS expressions
    """

    def __init__(self, build: GenomeBuild,
                 hgvs_coordinate_finder: VariantCoordinateFinder[str]):
        self._logger = logging.getLogger(__name__)
        self._build = hpotk.util.validate_instance(build, GenomeBuild, 'build')
        self._hgvs_finder = hpotk.util.validate_instance(hgvs_coordinate_finder, VariantCoordinateFinder,
                                                         'hgvs_coordinate_finder')

    def find_coordinates(self, item: GenomicInterpretation) -> typing.Tuple[VariantCoordinates, Genotype]:
        """Creates a VariantCoordinates object from the data in a given Phenopacket

        Args:
            item (GenomicInterpretation): a genomic interpretation element from Phenopacket Schema
        Returns:
            tuple[VariantCall, Genotype]: A tuple of :class:`VariantCoordinates` and :class:`Genotype`
        """
        if not isinstance(item, GenomicInterpretation):
            raise ValueError(f"item must be a Phenopacket GenomicInterpretation but was type {type(item)}")

        variant_descriptor = item.variant_interpretation.variation_descriptor

        vc = None
        if self._vcf_is_available(variant_descriptor.vcf_record):
            # We have a VCF record.
            if not self._check_assembly(variant_descriptor.vcf_record.genome_assembly):
                raise ValueError(f"Variant id {variant_descriptor.id} for patient {item.subject_or_biosample_id} has a different Genome Assembly than what was given. " \
                                 + f"{variant_descriptor.vcf_record.genome_assembly} is not {self._build.identifier}.")
            contig = self._build.contig_by_name(variant_descriptor.vcf_record.chrom)
            start = int(variant_descriptor.vcf_record.pos) - 1
            ref = variant_descriptor.vcf_record.ref
            alt = variant_descriptor.vcf_record.alt
            end = start + len(ref)
            change_length = end - start

            region = GenomicRegion(contig, start, end, Strand.POSITIVE)
            vc = VariantCoordinates(region, ref, alt, change_length)
        elif self._cnv_is_available(variant_descriptor.variation):
            # We have a CNV.
            variation = variant_descriptor.variation
            seq_location = variation.copy_number.allele.sequence_location
            refseq_contig_name = seq_location.sequence_id.split(':')[1]
            contig = self._build.contig_by_name(refseq_contig_name)

            # Assuming SV coordinates are 1-based (VCF style),
            # so we subtract 1 to transform to 0-based coordinate system
            start = int(seq_location.sequence_interval.start_number.value) - 1
            end = int(seq_location.sequence_interval.end_number.value)
            ref = 'N'
            number = int(variation.copy_number.number.value)
            if number > 2:
                alt = '<DUP>'
            else:
                alt = '<DEL>'
            change_length = end - start

            region = GenomicRegion(contig, start, end, Strand.POSITIVE)
            vc = VariantCoordinates(region, ref, alt, change_length)
        elif len(variant_descriptor.expressions) > 0:
            # We have some expressions. Let's try to find the 1st expression with `hgvs.c` syntax.
            for expression in variant_descriptor.expressions:
                if expression.syntax == 'hgvs.c':
                    vc = self._hgvs_finder.find_coordinates(expression.value)
                    break

        # Last, parse genotype.
        genotype = variant_descriptor.allelic_state.label
        gt = self._map_geno_genotype_label(genotype)

        return vc, gt

    def _check_assembly(self, genome_assembly:str) -> bool:
        if '38' in genome_assembly and self._build.identifier == 'GRCh38.p13':
            return True
        elif ('37' in genome_assembly or '19' in genome_assembly) and self._build.identifier == 'GRCh37.p13':
            return True
        else:
            return False


    @staticmethod
    def _vcf_is_available(vcf_record) -> bool:
        """
        Check if we can parse data out of VCF record.
        """
        return (vcf_record.genome_assembly != ''
                and vcf_record.chrom != ''
                and vcf_record.pos >= 0
                and vcf_record.ref != ''
                and vcf_record.alt != '')

    @staticmethod
    def _cnv_is_available(variation):
        seq_location = variation.copy_number.allele.sequence_location
        return (seq_location.sequence_id != ''
                and seq_location.sequence_interval.start_number.value >= 0
                and seq_location.sequence_interval.end_number.value >= 0
                and variation.copy_number.number.value != '')

    @staticmethod
    def _map_geno_genotype_label(genotype: str) -> Genotype:
        """
        Mapping from labels of the relevant GENO terms that is valid as of Oct 2nd, 2023.
        """
        if genotype in ('heterozygous', 'compound heterozygous', 'simple heterozygous'):
            return Genotype.HETEROZYGOUS
        elif genotype == 'homozygous':
            return Genotype.HOMOZYGOUS_ALTERNATE
        elif genotype in ('hemizygous', 'hemizygous X-linked', 'hemizygous Y-linked'):
            return Genotype.HEMIZYGOUS
        else:
            raise ValueError(f'Unknown genotype {genotype}')


class PhenopacketPatientCreator(PatientCreator[Phenopacket]):
    """
    `PhenopacketPatientCreator` transforms `Phenopacket` into :class:`genophenocorr.model.Patient`.
    """

    def __init__(self, build: GenomeBuild,
                 phenotype_creator: PhenotypeCreator,
                 var_func_ann: FunctionalAnnotator,
                 hgvs_coordinate_finder: VariantCoordinateFinder[str]):
        self._logger = logging.getLogger(__name__)
        # Violates DI, but it is specific to this class, so I'll leave it "as is".
        self._coord_finder = PhenopacketVariantCoordinateFinder(build, hgvs_coordinate_finder)
        self._phenotype_creator = hpotk.util.validate_instance(phenotype_creator, PhenotypeCreator, 'phenotype_creator')
        self._func_ann = hpotk.util.validate_instance(var_func_ann, FunctionalAnnotator, 'var_func_ann')

    def process(self, inputs: Phenopacket) -> AuditReport[Patient]:
        """Creates a Patient from the data in a given Phenopacket

        Args:
            inputs (Phenopacket): A Phenopacket object
        Returns:
            Patient: A Patient object
        """
        sample_id = SampleLabels(label=inputs.subject.id, meta_label=inputs.id if len(inputs.id) > 0 else None)

        # Validation - relay the errors found by `PhenotypeCreator`
        # we report issues in any case
        phenotypes = self._add_phenotypes(inputs)

        # Validate
        #  - there is >=1 variant in the subject
        #    - mitigate: skip, warning
        variants = self._add_variants(sample_id, inputs)

        # TODO: clean up when the protein metadata are removed from the patient.
        issues = []

        patient = Patient(sample_id, phenotypes=phenotypes.outcome, variants=variants.outcome, proteins=())
        return AuditReport(patient, issues)

    def _add_variants(self, sample_id: SampleLabels, pp: Phenopacket) -> AuditReport[typing.Sequence[Variant]]:
        """Creates a list of Variant objects from the data in a given Phenopacket

        Args:
            pp (Phenopacket): A Phenopacket object
        Returns:
            Sequence[Variant]: A list of Variant objects
        """
        # Variant validation:
        #
        #  - for each variant:
        #    - genomic interpretation must be a VCF record or a VRS CNV or an HGVS expression
        #      - mitigate: skip, warning
        #    - Lauren check the `vc.alt == 'N'` if it ever happens. I think we should not care if we do not have this
        #      in our phenopackets. Maybe a regex check.
        #    - there is at least 1 tx annotation. This can happen, since the `FunctionalAnnotator` is allowed to return
        #      an empty sequence. This has pros and cons, but I judge that it is better to give the annotator some wiggle
        #      room than to be really really really strict.
        #      - mitigate: ADD the variant but emit a warning. We can still use the variant in the `VariantPredicate`
        variants = []
        issues = []
        for i, interp in enumerate(pp.interpretations):
            for genomic_interp in interp.diagnosis.genomic_interpretations:
                vc, gt = self._coord_finder.find_coordinates(genomic_interp)
                if vc is None:
                    self._logger.warning('Expected a VCF record, a VRS CNV, or an expression with `hgvs.c` '
                                         'but did not find one in patient %s', sample_id)
                    continue

                if "N" in vc.alt:
                    self._logger.warning(
                        'Patient %s has unknown alternative variant %s, this variant will not be included.', pp.id,
                        vc.variant_key)
                    continue

                tx_annotations = self._func_ann.annotate(vc)
                if len(tx_annotations) == 0:
                    self._logger.warning("Patient %s has an error with variant %s, this variant will not be included.",
                                         pp.id, vc.variant_key)
                    continue

                genotype = Genotypes.single(sample_id, gt)
                variants.append(Variant(vc, tx_annotations, genotype))

        return AuditReport(variants, issues)

    def _add_phenotypes(self, pp: Phenopacket) -> AuditReport[typing.Sequence[Phenotype]]:
        """Creates a list of Phenotype objects from the data in a given Phenopacket

        Args:
            pp (Phenopacket): A Phenopacket object
        Returns:
            Sequence[Phenotype]: A list of Phenotype objects
        """
        hpo_id_list = []
        for hpo_id in pp.phenotypic_features:
            hpo_id_list.append((hpo_id.type.id, not hpo_id.excluded))

        return self._phenotype_creator.process(hpo_id_list)

    def _add_protein_data(self, variants: typing.Sequence[Variant]) -> typing.Collection[ProteinMetadata]:
        """Creates a list of ProteinMetadata objects from a given list of Variant objects

        Args:
            variants (Sequence[Variant]): A list of Variant objects
        Returns:
            Collection[ProteinMetadata]: A list of ProteinMetadata objects
        """
        final_prots = set()
        for var in variants:
            for trans in var.tx_annotations:
                for prot in trans.protein_affected:
                    final_prots.add(prot)
        return final_prots


def load_phenopacket_folder(pp_directory: str,
                            patient_creator: PhenopacketPatientCreator,
                            # policy:str = "strict/lenient/etc",
                            include_patients_with_no_HPO: bool = False) -> Cohort:
    """
    Creates a Patient object for each phenopacket formatted JSON file in the given directory `pp_directory`.

    :param pp_directory: path to a folder with phenopacket JSON files. An error is raised if the path does not point to
      a directory with at least one phenopacket.
    :param patient_creator: patient creator for turning a phenopacket into a :class:`genophenocorr.Patient`
    :return: a cohort made of the phenopackets
    """
    ## NOTE: Strict would stop code at error vs lenient would just print error but continue.
    if not os.path.isdir(pp_directory):
        raise ValueError("Could not find directory of Phenopackets.")
    hpotk.util.validate_instance(patient_creator, PhenopacketPatientCreator, 'patient_creator')

    # load Phenopackets
    pps = _load_phenopacket_dir(pp_directory)
    if len(pps) == 0:
        raise ValueError(f"No JSON Phenopackets were found in {pp_directory}")

    # turn phenopackets into patients using patient creator
    patients = []
    for pp in tqdm(pps, desc='Patients Created'):
        patient = patient_creator.create_patient(pp)
        if len(patient.phenotypes) == 0:
            patient_creator._logger.warning('Patient %s has no phenotypes listed and will not be included in this analysis.', patient.patient_id)
        if len(patient.variants) == 0:
            patient_creator._logger.warning('Patient %s has no variants listed and will not be included in this analysis.', patient.patient_id)
        if patient is not None:
            patients.append(patient)

    # create cohort from patients
    return Cohort.from_patients(patients, include_patients_with_no_HPO)


def _load_phenopacket_dir(pp_dir: str) -> typing.Sequence[Phenopacket]:
    patients = []
    for patient_file in os.listdir(pp_dir):
        if patient_file.endswith('.json'):
            phenopacket_path = os.path.join(pp_dir, patient_file)
            pp = load_phenopacket(phenopacket_path)
            patients.append(pp)
    return patients


def load_phenopacket(phenopacket_path: str) -> Phenopacket:
    """
    Load phenopacket JSON file.

    :param phenopacket_path: a `str` pointing to phenopacket JSON file.
    """
    with open(phenopacket_path) as f:
        return Parse(f.read(), Phenopacket())

