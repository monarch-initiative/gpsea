import pytest

import typing

from genophenocorr.model import Variant, Cohort

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
        
    @pytest.mark.parametrize(
        "search, expected",
        [
            (56002673, "12_56002674_56002674_T_C"),
            (78386857, "16_78386858_78425054_--38197bp--_A"),
            ('SO:1000029', 'SO:1000037_HGNC:21316_ANKRD11')
        ] # Couldn't find an example that produces a variant key similar to this: 22_10001_20000_INV
          # That is an example in the documentation though, so is it possible? (specifically the "INV"/"DEL"/etc.)
    )
    def test_variant_key(self, suox_cohort: Cohort, search, expected: typing.Optional[str]):
        
            for var in suox_cohort.all_variants():
                if var.variant_info.variant_coordinates is not None:
                    if var.variant_info.variant_coordinates.start == search:
                        assert var.variant_info.variant_key == expected
                if var.variant_info.variant_coordinates is None:
                    if var.variant_info.sv_info.structural_type == search:
                        assert var.variant_info.variant_key == expected
                
        