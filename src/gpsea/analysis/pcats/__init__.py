"""
The `gpsea.analysis.pcats` tests the association between genotype and phenotype groups,
if the groups can be defined in terms of discrete, unique, and non-overlapping categories.

Each individual is assigned into a genotype and phenotype group
using :class:`~gpsea.analysis.predicate.genotype.GenotypePolyPredicate`
and :class:`~gpsea.analysis.predicate.phenotype.PhenotypePolyPredicate` respectively.

A contingency matrix with group counts is prepared
and the counts are tested for association using :class:`~gpsea.analysis.pcats.stats.CountStatistic`,
such as :class:`~gpsea.analysis.pcats.stats.ScipyFisherExact`.

It is typical to test several phenotype groups at the same time.
Therefore, we must correct for multiple testing to prevent false positive findings.
See :ref:`MTC section <mtc>` for more info.

The results are provided as :class:`MultiPhenotypeAnalysisResult`
(or more specific :class:`HpoTermAnalysisResult` for :class:`HpoTermAnalysis`).
"""

from ._impl import MultiPhenotypeAnalysis, MultiPhenotypeAnalysisResult
from ._impl import DiseaseAnalysis
from ._impl import HpoTermAnalysis, HpoTermAnalysisResult
from ._impl import apply_predicates_on_patients

__all__ = [
    'MultiPhenotypeAnalysis', 'MultiPhenotypeAnalysisResult',
    'DiseaseAnalysis',
    'HpoTermAnalysis', 'HpoTermAnalysisResult',
    'apply_predicates_on_patients',
]
