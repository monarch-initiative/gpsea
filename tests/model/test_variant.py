import pytest

import typing

from genophenocorr.model import Variant, Cohort
from genophenocorr.view import VariantFormatter


class TestVariant:
    
    @pytest.fixture
    def test_cohort(
        self, 
        suox_cohort: Cohort
    ) -> Cohort:
        return suox_cohort

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
        
    @pytest.mark.parametrize(
        "variant, expected",
        [
            ('12_56004525_56004525_A_G', "NM_001032386.2:c.1136A>G"),
        ]
    )
    def test_variant_formatter(
        self, 
        variant: str,
        expected: str,
        test_cohort: Cohort
    ):
        var = test_cohort.get_variant_by_key(variant)
        formatter = VariantFormatter()
        assert formatter.format_as_string(var) == expected
