import pytest

from genophenocorr.analysis import CohortAnalysisConfiguration, MTC_Strategy


class TestCohortAnalysisConfiguration:

    def test_default_values(self):
        config = CohortAnalysisConfiguration()

        assert config.missing_implies_excluded == False
        assert config.pval_correction == 'bonferroni'
        assert config.min_perc_patients_w_hpo == pytest.approx(.2)
        assert config.include_sv == False
        assert config.mtc_alpha == pytest.approx(.05)
        assert config.mtc_strategy == MTC_Strategy.ALL_HPO_TERMS
        assert config.terms_to_test is None
    

    @pytest.mark.parametrize(
        'value',
        [
            -.0,
            1.01,
        ]
    )
    def test_cannot_set_invalid_threshold_in_heuristic_strategy(
        self,
        value: float,
    ):
        config = CohortAnalysisConfiguration()
        
        with pytest.raises(ValueError) as e:
            config.heuristic_strategy(value)

        assert e.value.args[0] == f'`threshold_HPO_observed_frequency` must be in range (0, 1] but was {value}'
