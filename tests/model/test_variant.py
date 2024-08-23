import pytest

import typing

from gpsea.model import Variant, Cohort, VariantCoordinates
from gpsea.model.genome import GenomeBuild, GenomicRegion, Strand


class TestVariant:

    @pytest.fixture
    def some_variant(
        self,
        suox_cohort: Cohort,
    ) -> Variant:
        return suox_cohort.get_variant_by_key('12_56004525_56004525_A_G')

    @pytest.mark.parametrize(
            "tx_id, expected", 
            [
                ("NM_001032386.2", 'NM_001032386.2:c.1136A>G'),
                ("Whatever", None),
            ])
    def test_get_hgvs_cdna_by_tx(
        self,
        some_variant: Variant,
        tx_id: str,
        expected: typing.Optional[str],
    ):
        hgvs = some_variant.get_hgvs_cdna_by_tx_id(transcript_id=tx_id)

        assert hgvs == expected


class TestVariantCoordinates:

    @pytest.mark.parametrize(
        "contig_name, start, end, ref, alt, change_length, expected",
        [
            ("chrX", 100, 101, "C", "T", 0, "X_101_101_C_T"),
            ("chrY", 150, 152, "G", "GG", 1, "Y_151_152_G_GG"),
            ("chr16", 1000, 1040, "A", "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG", 39, "16_1001_1040_A_--40bp--"),
            ("chr2", 200, 301, "N", "<DEL>", 100, "2_201_301_DEL")
        ]
    )
    def test_variant_key(
        self,
        genome_build: GenomeBuild,
        contig_name: str,
        start: int, end: int,
        ref: str, alt: str,
        change_length: int,
        expected: str,
    ):
        contig = genome_build.contig_by_name(contig_name)
        assert contig is not None

        vc = VariantCoordinates(
            region=GenomicRegion(
                contig=contig,
                start=start, end=end, strand=Strand.POSITIVE,
            ),
            ref=ref, alt=alt, change_length=change_length,
        )

        assert vc.variant_key == expected
