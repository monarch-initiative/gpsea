import pytest

from genophenocorr.analysis import CohortAnalysisConfiguration, MtcStrategy


class TestCohortAnalysisConfiguration:

    def test_default_values(self):
        config = CohortAnalysisConfiguration()

        assert config.missing_implies_excluded is False
        assert config.pval_correction == 'bonferroni'
        assert config.min_patients_w_hpo is None
        assert config.include_sv is False
        assert config.mtc_alpha == pytest.approx(.05)
        assert config.mtc_strategy == MtcStrategy.ALL_PHENOTYPE_TERMS
        assert config.terms_to_test is None
    
    def test_set_all_terms_strategy(self):
        config = CohortAnalysisConfiguration()
        assert config.mtc_strategy == MtcStrategy.ALL_PHENOTYPE_TERMS

        config.specify_terms_strategy(('HP:0001250', 'HP:0001166'))

        assert config.mtc_strategy == MtcStrategy.SPECIFY_TERMS
        assert config.terms_to_test == ('HP:0001250', 'HP:0001166')

    def test_set_hpo_mtc_strategy(self):
        config = CohortAnalysisConfiguration()
        assert config.mtc_strategy == MtcStrategy.ALL_PHENOTYPE_TERMS

        config.hpo_mtc_strategy()

        assert config.mtc_strategy == MtcStrategy.HPO_MTC
        assert config.min_patients_w_hpo == pytest.approx(0.2)

    @pytest.mark.parametrize(
        'value',
        [
            -.0,
            1.01,
        ]
    )
    def test_cannot_set_invalid_threshold_in_hpo_mtc_strategy(
        self,
        value: float,
    ):
        config = CohortAnalysisConfiguration()
        
        with pytest.raises(ValueError) as e:
            config.hpo_mtc_strategy(value)

        assert e.value.args[0] == f'`min_patients_w_hpo` must be in range (0, 1] but was {value}'
