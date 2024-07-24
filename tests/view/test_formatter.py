import pytest

import typing

from genophenocorr.model import Cohort
from genophenocorr.view import VariantFormatter

class TestFormatter:
    
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
        suox_cohort: Cohort
    ):
        var = suox_cohort.get_variant_by_key(variant)
        formatter = VariantFormatter()
        assert formatter.format_as_string(var) == expected