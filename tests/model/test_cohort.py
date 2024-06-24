import pytest

from genophenocorr.model import Cohort, VariantEffect


class TestCohort:

    def test_all_transcript_ids(
        self,
        suox_cohort: Cohort,
    ):
        # The SUOX cohort includes variants that affect the following transcripts:
        assert suox_cohort.all_transcript_ids == {
            "NM_001032386.2",  # MANE
            "NM_001351091.2",            
            "NM_001032387.2",
            "NM_001351089.2",
            "NM_000456.3",
        }

    def test_variant_effect_count_by_tx(
        self,
        suox_cohort: Cohort,
    ):
        counts = suox_cohort.variant_effect_count_by_tx()

        # We should get variant effect counts wrt. all transcripts
        assert set(counts.keys()) == suox_cohort.all_transcript_ids

        # Let's check counts on the SUOX MANE transcript
        suox_mane = 'NM_001032386.2'
        assert suox_mane in counts
        
        mane_counts = counts[suox_mane]

        # print(mane_counts)
        assert mane_counts == {'MISSENSE_VARIANT': 29, 'STOP_GAINED': 10, 'FRAMESHIFT_VARIANT': 9}

    def test_variant_effect_count_by_tx__singular(
            self,
            suox_cohort: Cohort,
    ):
        suox_mane = 'NM_001032386.2'
        counts = suox_cohort.variant_effect_count_by_tx(tx_id=suox_mane)

        assert suox_mane in counts
        assert len(counts) == 1, 'The counts should only have one item'

        assert counts[suox_mane] == {'MISSENSE_VARIANT': 29, 'STOP_GAINED': 10, 'FRAMESHIFT_VARIANT': 9}
