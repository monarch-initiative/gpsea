import json
import os

import hpotk
import pytest

from gpsea.model import Cohort
from gpsea.io import GpseaJSONDecoder
from gpsea.analysis.predicate.genotype import VariantPredicates, monoallelic_predicate
from gpsea.analysis.tempo import SurvivalAnalysis
from gpsea.analysis.tempo.endpoint import hpo_onset
from gpsea.analysis.tempo.stats import LogRankTest


@pytest.fixture
def umod_cohort(
    fpath_cohort_data_dir: str,
) -> Cohort:
    fpath_cohort = os.path.join(fpath_cohort_data_dir, "UMOD.0.1.20.json")
    with open(fpath_cohort) as fh:
        return json.load(fh, cls=GpseaJSONDecoder)


class TestSurvivalAnalysis:

    @pytest.fixture(scope="class")
    def survival_analysis(self) -> SurvivalAnalysis:
        return SurvivalAnalysis(statistic=LogRankTest())

    def test_compare_genotype_vs_survival(
        self,
        survival_analysis: SurvivalAnalysis,
        hpo: hpotk.MinimalOntology,
        umod_cohort: Cohort,
    ):
        in_exon_3 = VariantPredicates.exon(3, tx_id="NM_003361.4")
        gt_predicate = monoallelic_predicate(
            a_predicate=in_exon_3,
            b_predicate=~in_exon_3,
            names=("Exon 3", "Other exon")
        )
        endpoint = hpo_onset(
            hpo=hpo,
            term_id="HP:0003774",  # Stage 5 chronic kidney disease
        )
        
        result = survival_analysis.compare_genotype_vs_survival(
            cohort=umod_cohort,
            gt_predicate=gt_predicate,
            endpoint=endpoint,
        )

        print(result)
