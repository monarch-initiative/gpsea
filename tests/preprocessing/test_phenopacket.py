import os

import pytest

from google.protobuf.json_format import Parse
from phenopackets.schema.v2.core.interpretation_pb2 import GenomicInterpretation

from genophenocorr.model.genome import GenomeBuild

from genophenocorr.preprocessing import VVHgvsVariantCoordinateFinder
from genophenocorr.preprocessing import (
    VariantCoordinateFinder,
    PhenopacketVariantCoordinateFinder,
)


class TestPhenopacketVariantCoordinateFinder:


    @pytest.fixture(scope='class')
    def fpath_test_genomic_interpretations(
        self,
        fpath_preprocessing_data_dir: str,
    ) -> str:
        return os.path.join(fpath_preprocessing_data_dir, 'pp_genomic_interpretations')

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
        "pp_name, expected",
        [
            ("deletion_test.json", "16_89284129_89284134_CTTTTT_C"),
            ("insertion_test.json", "16_89280829_89280829_C_CA"),
            ("missense_test.json", "16_89279135_89279135_G_C"),
            ("missense_hgvs_test.json", "16_89279135_89279135_G_C"),
            ("duplication_test.json", "16_89279850_89279850_G_GC"),
            ("delinsert_test.json", "16_89284601_89284602_GG_A"),
            ("CVDup_test.json", "16_89284523_89373231_DUP"),
            ("CVDel_test.json", "16_89217281_89506042_DEL"),
        ],
    )
    def test_find_coordinates(
        self, 
        pp_name : str,
        fpath_test_genomic_interpretations: str,
        expected: str, 
        pp_vc_finder: PhenopacketVariantCoordinateFinder,
    ):
        fpath_pp = os.path.join(fpath_test_genomic_interpretations, pp_name)
        gi = read_genomic_interpretation_json(fpath_pp)

        vc, _ = pp_vc_finder.find_coordinates(gi)

        assert vc.variant_key == expected


def read_genomic_interpretation_json(fpath: str) -> GenomicInterpretation:
    with open(fpath) as fh:
        return Parse(fh.read(), GenomicInterpretation())
