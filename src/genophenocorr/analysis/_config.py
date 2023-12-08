import logging
import typing

import hpotk

from genophenocorr.model import Cohort

from ._api import CohortAnalysis
from ._commie import CommunistCohortAnalysis


P_VAL_OPTIONS = ['bonferroni', 'b', 
                 'sidak', 's', 
                 'holm-sidak', 'hs', 
                 'holm', 'h', 
                 'simes-hochberg', 'sh', 
                 'hommel', 'ho', 
                 'fdr_bh', 
                 'fdr_by', 
                 'fdr_tsbh', 
                 'fdr_tsbky',
                 'fdr_gbs',
                 None]

class CohortAnalysisConfiguration:
    """
    `CohortAnalysisConfiguration` is a value class for storing :class:`genophenocorr.analysis.CohortAnalysis`
    configuration options.

    Default values
    ^^^^^^^^^^^^^^

    ==============================  ===========  ====================
        Option                        Type           Default value
    ==============================  ===========  ====================
     ``missing_implies_excluded``    `bool`      `False`
     ``pval_correction``             `str`       `bonferroni`
     ``min_perc_patients_w_hpo``     `float`     `0.1`
     ``include_sv``                  `bool`      `False`
    ==============================  ===========  ====================

    """

    def __init__(self, missing_implies_excluded: bool,
                 pval_correction: typing.Optional[str],
                 min_perc_patients_w_hpo: float,
                 include_sv: bool):
        self._missing_implies_excluded = missing_implies_excluded
        self._pval_correction = pval_correction
        self._min_perc_patients_w_hpo = min_perc_patients_w_hpo
        self._include_sv = include_sv

    @staticmethod
    def builder():
        """
        Create a new :class:`CohortAnalysisConfigurationBuilder` initialized with default configuration options.
        """
        return CohortAnalysisConfigurationBuilder()

    @property
    def missing_implies_excluded(self) -> bool:
        """
        `True` if we assume that a patient without a specific phenotype listed *does not* have the phenotype.
        Otherwise, the only excluded phenotypes are those that are excluded explicitly.
        """
        return self._missing_implies_excluded

    @property
    def pval_correction(self) -> typing.Optional[str]:
        """
        Get method for multiple testing p value correction. Default: `bonferroni`.
        """
        return self._pval_correction

    @property
    def min_perc_patients_w_hpo(self) -> float:
        """
        A threshold for removing rare HPO terms, only the terms that are observed in at least this fraction of patients
        will be retained for the analysis.
        """
        return self._min_perc_patients_w_hpo

    @property
    def include_sv(self) -> bool:
        """
        `True` if we want to include structural variants in the analysis
        (i.e. the variants that use symbolic VCF notation).
        """
        return self._include_sv


class CohortAnalysisConfigurationBuilder:
    """
    `CohortAnalysisConfigurationBuilder` builds :class:`CohortAnalysisConfiguration`. The configuration options
    can be customized or set with default values.

    If an invalid option is passed, a warning is logged and the previous value is retained.
    Consequently, it is *impossible* to mis-configure the analysis.
    """

    def __init__(self):
        self._logger = logging.getLogger(__name__)
        self._missing_implies_excluded = False
        self._pval_correction = 'bonferroni'
        self._min_perc_patients_w_hpo = .1
        self._include_sv = False

    def missing_implies_excluded(self, missing_implies_excluded: bool):
        """
        Set `missing_implies_excluded` option.
        """
        if isinstance(missing_implies_excluded, bool):
            self._missing_implies_excluded = missing_implies_excluded
        else:
            self._logger.warning('Ignoring invalid `missing_implies_excluded` value %s',
                                 missing_implies_excluded)

        return self

    def pval_correction(self, pval_correction: typing.Optional[str]):
        """
        Set `pval_correction` option.
        """
        if pval_correction in P_VAL_OPTIONS:
            self._pval_correction = pval_correction
        else:
            self._logger.warning('Ignoring invalid `pval_correction` value %s. Using default "bonferroni" correction.', pval_correction)
        return self

    def min_perc_patients_w_hpo(self, min_perc_patients_w_hpo: float):
        """
        Set `min_perc_patients_w_hpo` option.
        """
        if not isinstance(min_perc_patients_w_hpo, float):
            try:
                min_perc_patients_w_hpo = float(min_perc_patients_w_hpo)
            except ValueError:
                self._logger.warning("min_perc_patients_w_hpo must be a number. Using default of 0.1")
                min_perc_patients_w_hpo = 0.1
        if min_perc_patients_w_hpo > 1 or min_perc_patients_w_hpo <= 0:
            self._logger.warning("min_perc_patients_w_hpo must be greater than 0 and at most 1. Using default of 0.1")
        else:
            self._min_perc_patients_w_hpo = min_perc_patients_w_hpo
        return self

    def include_sv(self, include_sv: bool):
        """
        Set `include_sv` option.
        """
        if isinstance(include_sv, bool):
            self._include_sv = include_sv
        else:
            self._logger.warning('Ignoring invalid `include_sv` value %s. Defaulting to not include large structural variants.', include_sv)
        return self

    def build(self) -> CohortAnalysisConfiguration:
        """
        Build the configuration.
        """
        return CohortAnalysisConfiguration(self._missing_implies_excluded,
                                           self._pval_correction,
                                           self._min_perc_patients_w_hpo,
                                           self._include_sv)


def configure_cohort_analysis(cohort: Cohort,
                              hpo: hpotk.MinimalOntology,
                              config: typing.Optional[CohortAnalysisConfiguration] = None) -> CohortAnalysis:
    """
    Configure :class:`genophenocorr.analysis.CohortAnalysis` for given `cohort`.

    :param cohort: a :class:`genophenocorr.model.Cohort` to analyze
    :param hpo: a :class:`hpotk.MinimalOntology` with HPO to use in the analysis
    :param config: an optional :class:`CohortAnalysisConfiguration` to parameterize the analysis.
     The default parameters will be used if `None`.
    """
    if config is None:
        config = CohortAnalysisConfiguration.builder().build()

    return CommunistCohortAnalysis(cohort, hpo,
                                   missing_implies_excluded=config.missing_implies_excluded,
                                   include_sv=config.include_sv,
                                   p_val_correction=config.pval_correction,
                                   min_perc_patients_w_hpo=config.min_perc_patients_w_hpo)
