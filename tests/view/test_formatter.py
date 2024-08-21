import pytest

from genophenocorr.model import Cohort
from genophenocorr.view import VariantFormatter


class TestVariantFormatter:
    
    @pytest.fixture(scope='class')
    def formatter(self) -> VariantFormatter:
        return VariantFormatter()

    @pytest.mark.parametrize(
        "variant, expected",
        [
            ('12_56004525_56004525_A_G', "NM_001032386.2:c.1136A>G"),
            # ('16_78386858_78425054_--38197bp--_A', 'NM_016373.4:c.517_791del')
        ]
    )
    def test_variant_formatter(
        self,
        variant: str,
        expected: str,
        formatter: VariantFormatter,
        suox_cohort: Cohort
    ):
        var = suox_cohort.get_variant_by_key(variant)
        
        assert formatter.format_as_string(var) == expected
