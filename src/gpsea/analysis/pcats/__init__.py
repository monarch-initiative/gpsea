"""
The `gpsea.analysis.pcats` tests the association between genotype and phenotype classes,
if the classess can be defined in terms of discrete, unique, and non-overlapping categories.

Each individual is assigned into a genotype and phenotype class
using :class:`~gpsea.analysis.clf.GenotypeClassifier`
and :class:`~gpsea.analysis.clf.PhenotypeClassifier` respectively.

A contingency matrix with group counts is prepared
and the counts are tested for association using :class:`~gpsea.analysis.pcats.stats.CountStatistic`,
such as :class:`~gpsea.analysis.pcats.stats.FisherExactTest`.

It is typical to test several phenotype groups at the same time.
Therefore, we must correct for multiple testing to prevent false positive findings.
See :ref:`MTC section <mtc>` for more info.

The results are provided as :class:`~gpsea.analysis.MultiPhenotypeAnalysisResult`
(or more specific :class:`~gpsea.analysis.pcats.HpoTermAnalysisResult`
for :class:`~gpsea.analysis.pcats.HpoTermAnalysis`).

Use :func:`~gpsea.analysis.pcats.configure_hpo_term_analysis` to configure the HPO term analysis with the default parameters.
"""

from ._impl import MultiPhenotypeAnalysis, MultiPhenotypeAnalysisResult
from ._impl import DiseaseAnalysis
from ._impl import HpoTermAnalysis, HpoTermAnalysisResult
from ._impl import apply_classifiers_on_individuals
from ._config import configure_hpo_term_analysis

__all__ = [
    "MultiPhenotypeAnalysis",
    "MultiPhenotypeAnalysisResult",
    "DiseaseAnalysis",
    "HpoTermAnalysis",
    "HpoTermAnalysisResult",
    "apply_classifiers_on_individuals",
    "configure_hpo_term_analysis",
]
