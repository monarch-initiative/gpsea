import os

import pytest

from google.protobuf.json_format import Parse
from phenopackets.schema.v2.core.interpretation_pb2 import GenomicInterpretation

from gpsea.model.genome import GenomeBuild, Strand

from gpsea.preprocessing import VVHgvsVariantCoordinateFinder
from gpsea.preprocessing import (
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
        fpath_pp = os.path.join(fpath_test_genomic_interpretations, 'chromosomal_deletion.ANKRD11.json')
        gi = read_genomic_interpretation_json(fpath_pp)

        vc = pp_vc_finder.find_coordinates(gi)
        assert vc is None


def read_genomic_interpretation_json(fpath: str) -> GenomicInterpretation:
    with open(fpath) as fh:
        return Parse(fh.read(), GenomicInterpretation())
