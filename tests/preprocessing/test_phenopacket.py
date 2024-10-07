import os

import hpotk
import pytest

from google.protobuf.json_format import Parse
from phenopackets.schema.v2.core.interpretation_pb2 import GenomicInterpretation
from phenopackets.schema.v2.phenopackets_pb2 import Phenopacket

from gpsea.model.genome import GenomeBuild, Strand

from gpsea.preprocessing import VVHgvsVariantCoordinateFinder
from gpsea.preprocessing import (
    VariantCoordinateFinder,
    PhenopacketVariantCoordinateFinder,
)
from gpsea.preprocessing import (
    FunctionalAnnotator,
    VarCachingFunctionalAnnotator,
    VepFunctionalAnnotator,
    ImpreciseSvFunctionalAnnotator,
    DefaultImpreciseSvFunctionalAnnotator,
)
from gpsea.preprocessing import PhenopacketPatientCreator, PhenotypeCreator
from gpsea.preprocessing import VVMultiCoordinateService


class TestPhenopacketVariantCoordinateFinder:

    @pytest.fixture(scope="class")
    def fpath_test_genomic_interpretations(
        self,
        fpath_preprocessing_data_dir: str,
    ) -> str:
        return os.path.join(fpath_preprocessing_data_dir, "pp_genomic_interpretations")

    @pytest.fixture(scope="class")
    def hgvs_vc_finder(
        self,
        genome_build: GenomeBuild,
    ) -> VariantCoordinateFinder:
        return VVHgvsVariantCoordinateFinder(genome_build)

    @pytest.fixture(scope="class")
    def pp_vc_finder(
        self,
        genome_build: GenomeBuild,
        hgvs_vc_finder: VariantCoordinateFinder,
    ) -> PhenopacketVariantCoordinateFinder:
        return PhenopacketVariantCoordinateFinder(genome_build, hgvs_vc_finder)

    @pytest.mark.online
    @pytest.mark.parametrize(
        "pp_name, contig, start, end, ref, alt, change_length",
        [
            (
                "deletion_test.json",
                "16", 89284128, 89284134, "CTTTTT", "C", -5,
            ),
            (
                "insertion_test.json",
                "16", 89280828, 89280829, "C", "CA", 1,
            ),
            (
                "missense_test.json",
                "16", 89279134, 89279135, "G", "C", 0,
            ),
            (
                "missense_hgvs_test.json",
                "16", 89279134, 89279135, "G", "C", 0,
            ),
            (
                "duplication_test.json",
                "16", 89279849, 89279850, "G", "GC", 1,
            ),
            (
                "delinsert_test.json",
                "16", 89284600, 89284602, "GG", "A", -1,
            ),
            (
                "CVDup_test.json",
                "16", 89_284_522, 89_373_231, "N", "<DUP>", 88_709,
            ),
            (
                "CVDel_test.json",
                "16", 89_217_280, 89_506_042, "N", "<DEL>", -288_762,
            ),
        ],
    )
    def test_find_coordinates(
        self,
        pp_name: str,
        contig: str,
        start: int,
        end: int,
        ref: str,
        alt: str,
        change_length: int,
        fpath_test_genomic_interpretations: str,
        pp_vc_finder: PhenopacketVariantCoordinateFinder,
    ):
        fpath_pp = os.path.join(fpath_test_genomic_interpretations, pp_name)
        gi = read_genomic_interpretation_json(fpath_pp)

        vc = pp_vc_finder.find_coordinates(gi)

        assert vc is not None

        assert vc.chrom == contig
        assert vc.start == start
        assert vc.end == end
        assert vc.region.strand == Strand.POSITIVE
        assert vc.ref == ref
        assert vc.alt == alt
        assert vc.change_length == change_length

    def test_find_large_structural(
        self,
        fpath_test_genomic_interpretations: str,
        pp_vc_finder: PhenopacketVariantCoordinateFinder,
    ):
        fpath_pp = os.path.join(
            fpath_test_genomic_interpretations, "chromosomal_deletion.ANKRD11.json"
        )
        gi = read_genomic_interpretation_json(fpath_pp)

        vc = pp_vc_finder.find_coordinates(gi)
        assert vc is None


def read_genomic_interpretation_json(fpath: str) -> GenomicInterpretation:
    with open(fpath) as fh:
        return Parse(fh.read(), GenomicInterpretation())


class TestPhenopacketPatientCreator:

    @pytest.fixture
    def phenotype_creator(
        self,
        hpo: hpotk.MinimalOntology,
        validation_runner: hpotk.validate.ValidationRunner,
    ) -> PhenotypeCreator:
        return PhenotypeCreator(
            hpo=hpo,
            validator=validation_runner,
        )

    @pytest.fixture
    def functional_annotator(
        self,
        fpath_project_dir: str,
    ) -> FunctionalAnnotator:
        fpath_cache_dir = os.path.join(fpath_project_dir, ".gpsea_cache")
        fpath_variant_cache_dir = os.path.join(fpath_cache_dir, "variant_cache")
        os.makedirs(fpath_variant_cache_dir, exist_ok=True)

        return VarCachingFunctionalAnnotator.with_cache_folder(
            fpath_cache_dir=fpath_variant_cache_dir,
            fallback=VepFunctionalAnnotator(
                timeout=20,
            ),
        )

    @pytest.fixture
    def imprecise_sv_functional_annotator(
        self,
        genome_build: GenomeBuild,
    ) -> ImpreciseSvFunctionalAnnotator:
        return DefaultImpreciseSvFunctionalAnnotator(
            gene_coordinate_service=VVMultiCoordinateService(
                genome_build=genome_build,
            ),
        )

    @pytest.fixture
    def variant_coordinate_finder(
        self,
        genome_build: GenomeBuild,
    ) -> VariantCoordinateFinder:
        return VVHgvsVariantCoordinateFinder(
            genome_build=genome_build,
        )

    @pytest.fixture
    def patient_creator(
        self,
        genome_build: GenomeBuild,
        phenotype_creator: PhenotypeCreator,
        functional_annotator: FunctionalAnnotator,
        imprecise_sv_functional_annotator: ImpreciseSvFunctionalAnnotator,
        variant_coordinate_finder: VariantCoordinateFinder,
    ) -> PhenopacketPatientCreator:
        return PhenopacketPatientCreator(
            build=genome_build,
            phenotype_creator=phenotype_creator,
            functional_annotator=functional_annotator,
            imprecise_sv_functional_annotator=imprecise_sv_functional_annotator,
            hgvs_coordinate_finder=variant_coordinate_finder,
        )

    @pytest.fixture
    def phenopacket(
        self,
        fpath_phenopacket_dir: str,
    ) -> Phenopacket:
        fpath_pp = os.path.join(
            fpath_phenopacket_dir, "PMID_30968594_individual_1.json"
        )

        with open(fpath_pp) as fh:
            return Parse(fh.read(), Phenopacket())

    def test_phenopacket_patient_creator(
        self,
        phenopacket: Phenopacket,
        patient_creator: PhenopacketPatientCreator,
    ):
        notepad = patient_creator.prepare_notepad("A phenopacket")
        patient = patient_creator.process(phenopacket, notepad)

        # No issues
        assert not notepad.has_errors_or_warnings(include_subsections=True)

        # Individual credentials are OK
        assert (
            patient.labels.label_summary() == "individual 1[PMID_30968594_individual_1]"
        )
        assert patient.sex.is_male()

        # 6 Phenotype features
        assert len(patient.phenotypes) == 6
        assert tuple(p.identifier.value for p in patient.phenotypes) == (
            "HP:0000953",
            "HP:0030088",
            "HP:0003154",
            "HP:0008163",
            "HP:0000870",
            "HP:0025133",
        )
        assert tuple(p.is_present for p in patient.phenotypes) == (
            True,
            True,
            True,
            True,
            True,
            False,
        )

        # 6 measurements
        assert len(patient.measurements) == 6
        assert tuple(m.identifier.value for m in patient.measurements) == (
            "LOINC:1668-3",
            "LOINC:2986-8",
            "LOINC:2141-0",
            "LOINC:2143-6",
            "LOINC:2842-3",
            "LOINC:2243-4",
        )

        units = tuple(m.unit.value for m in patient.measurements)
        assert units == (
            "UCUM:ng/dL",
            "UCUM:ng/dL",
            "UCUM:ng/dL",
            "UCUM:ng/dL",
            "UCUM:ng/dL",
            "UCUM:ng/dL",
        )

        values = tuple(m.test_result for m in patient.measurements)
        assert values == (800.0, 127.0, 180.2, 116.6, 52.93, 23.71)

    # TODO: test disease and variants/interpretations
