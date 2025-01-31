import typing

import hpotk

from ..mtc_filter import IfHpoFilter
from ._impl import HpoTermAnalysis
from .stats import CountStatistic, FisherExactTest


def configure_hpo_term_analysis(
    hpo: hpotk.MinimalOntology,
    count_statistic: CountStatistic = FisherExactTest(),
    mtc_correction: typing.Optional[str] = "fdr_bh",
    mtc_alpha: float = 0.05,
) -> HpoTermAnalysis:
    """
    Configure HPO term analysis with default parameters.

    The default analysis will pre-filter HPO terms with :class:`~gpsea.analysis.mtc_filter.IfHpoFilter`,
    then compute nominal p values using `count_statistic` (default Fisher exact test),
    and apply multiple testing correction (default Benjamini/Hochberg (`fdr_bh`))
    with target `mtc_alpha` (default `0.05`).
    """
    return HpoTermAnalysis(
        mtc_filter=IfHpoFilter.default_filter(hpo),
        count_statistic=count_statistic,
        mtc_correction=mtc_correction,
        mtc_alpha=mtc_alpha,
    )
