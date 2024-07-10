import logging
import os
import typing
import enum
import hpotk

from genophenocorr.model import Cohort
from genophenocorr.preprocessing import ProteinMetadataService, UniprotProteinMetadataService, ProteinAnnotationCache, \
    ProtCachingMetadataService
from ._api import CohortAnalysis, HpoMtcFilter
from ._filter import SimplePhenotypeFilter
from ._gp_analysis import FisherExactAnalyzer, IdentityTermMtcFilter, HeuristicMtcFilter, SpecifiedTermsMtcFilter
from ._gp_impl import GpCohortAnalysis

P_VAL_OPTIONS = (
    'bonferroni', 'b',
    'sidak', 's',
    'holm-sidak', 'hs',
    'holm', 'h',
    'simes-hochberg', 'sh',
    'hommel', 'ho',
    'fdr_bh', 'fdr_by',
    'fdr_tsbh', 'fdr_tsbky',
    'fdr_gbs',
    None,
)

class MTC_Strategy(enum.Enum):
    ALL_HPO_TERMS = 0
    HEURISTIC_SAMPLER = 1
    SPECIFY_TERMS = 2




class CohortAnalysisConfiguration:
    """
    `CohortAnalysisConfiguration` is a value class for storing :class:`genophenocorr.analysis.CohortAnalysis`
    configuration options.

    The class contains the default values upon creation and the configuration option values can be set as properties.

    If an invalid value option is passed to the property setter, a warning is logged and the previous value is retained.
    Therefore, it is *impossible* to mis-configure the analysis.

    Default values
    ^^^^^^^^^^^^^^

    ==============================  =======================  =====================================
        Option                        Type                    Default value
    ==============================  =======================  =====================================
     ``missing_implies_excluded``    `bool`                   `False`
     ``pval_correction``             `str`                    `bonferroni`
     ``min_perc_patients_w_hpo``     `float`                  `0.1`
     ``mtc_alpha``                   `float`                  `0.05`
     ``include_sv``                  `bool`                   `False`
     ``mtc_strategy``                :class:`MTC_Strategy`    :class:`MTC_Strategy.ALL_HPO_TERMS`
    ==============================  =======================  =====================================

    """

    def __init__(self):
        self._logger = logging.getLogger(__name__)
        self._missing_implies_excluded = False
        self._pval_correction = 'bonferroni'
        self._mtc_alpha = .05
        self._min_perc_patients_w_hpo = .2
        self._include_sv = False
        self._mtc_strategy = MTC_Strategy.ALL_HPO_TERMS
        self._terms_to_test = set() ## only relevant for SPECIFIED_TERMS strategy

    @property
    def missing_implies_excluded(self) -> bool:
        """
        `True` if we assume that a patient without a specific phenotype listed *does not* have the phenotype.
        Otherwise, the only excluded phenotypes are those that are excluded explicitly.
        """
        return self._missing_implies_excluded

    @missing_implies_excluded.setter
    def missing_implies_excluded(self, missing_implies_excluded: bool):
        """
        :param missing_implies_excluded: a `bool` to dictate is the missing phenotypic features should be considered
          as excluded.
        """
        if isinstance(missing_implies_excluded, bool):
            self._missing_implies_excluded = missing_implies_excluded
        else:
            self._logger.warning('Ignoring invalid `missing_implies_excluded` value %s. Using %s.',
                                 missing_implies_excluded, self._missing_implies_excluded)

    @property
    def pval_correction(self) -> typing.Optional[str]:
        """
        Get method for multiple testing p value correction. Default: `bonferroni`.
        """
        return self._pval_correction

    @pval_correction.setter
    def pval_correction(self, pval_correction: typing.Optional[str]):
        """
        :param pval_correction: a `str` with the name of the desired multiple testing correction method or `None`
          if no MTC should be applied.
        """
        if pval_correction in P_VAL_OPTIONS:
            self._pval_correction = pval_correction
        else:
            self._logger.warning('Ignoring invalid `pval_correction` value %s. Using %s correction.', pval_correction,
                                 self._pval_correction)

    @property
    def min_perc_patients_w_hpo(self) -> float:
        """
        A threshold for removing rare HPO terms, only the terms that are observed in at least this fraction of patients
        will be retained for the analysis (default 0.2).
        """
        return self._min_perc_patients_w_hpo

    @min_perc_patients_w_hpo.setter
    def min_perc_patients_w_hpo(self, value: float):
        """
        Set new multiple testing correction alpha value.

        :param mtc_alpha: a `float` in range :math:`(0,1]`.
        """
        if isinstance(value, float) and 0. < value <= 1.:
            self._min_perc_patients_w_hpo = value
        else:
            self._logger.warning(
                '`value` should be a `float` in range `(0, 1]` but was %s. Keeping the previous value %s',
                value, self._min_perc_patients_w_hpo
            )


    @property
    def mtc_alpha(self) -> float:
        """
        The alpha value for multiple testing correction.
        """
        return self._mtc_alpha

    @mtc_alpha.setter
    def mtc_alpha(self, mtc_alpha: float):
        """
        Set new multiple testing correction alpha value.

        :param mtc_alpha: a `float` in range :math:`(0,1]`.
        """
        if isinstance(mtc_alpha, float) and 0. < mtc_alpha <= 1.:
            self._mtc_alpha = mtc_alpha
        else:
            self._logger.warning('`mtc_alpha` should be a `float` in range `(0, 1]` but was %s. Keeping the previous value %s',
                                 mtc_alpha, self._mtc_alpha)

    @property
    def include_sv(self) -> bool:
        """
        `True` if we want to include structural variants in the analysis
        (i.e. the variants that use symbolic VCF notation).
        """
        return self._include_sv

    @include_sv.setter
    def include_sv(self, include_sv: bool):
        """
        Set `include_sv` option.

        :param include_sv: a `bool` with the value.
        """
        if isinstance(include_sv, bool):
            self._include_sv = include_sv
        else:
            self._logger.warning('Ignoring invalid `include_sv` value %s. Using %s', include_sv, self._include_sv)

    @property
    def mtc_strategy(self) -> MTC_Strategy:
        """
        TODO
        """
        return self._mtc_strategy

    @mtc_strategy.setter
    def mtc_strategy(self, strategy: MTC_Strategy):
        """
        TODO
        """
        # TODO: this needs to be removed to prevent inconsistent configuration.
        # If we use:
        # ```
        # config.mtc_strategy = MTC_Strategy.SPECIFY_TERMS
        # ````
        # then the user will get an errror because of missing terms to test.
        if isinstance(strategy, MTC_Strategy):
            self._mtc_strategy = strategy
        else:
            self._logger.warning('TODO DOCUMENT `include_sv` value %s. Using %s', strategy, self._mtc_strategy)

    def heuristic_strategy(self, threshold_HPO_observed_frequency: float=0.2):
        self.mtc_strategy = MTC_Strategy.HEURISTIC_SAMPLER
        if not isinstance(threshold_HPO_observed_frequency, float):
            try:
                threshold_HPO_observed_frequency = float(threshold_HPO_observed_frequency)
            except ValueError:
                self._logger.warning("min_perc_patients_w_hpo must be a number, but was %s. Using %f",
                                     threshold_HPO_observed_frequency, self._min_perc_patients_w_hpo)
        if threshold_HPO_observed_frequency > 1 or threshold_HPO_observed_frequency <= 0:
            errmsg = f"min_perc_patients_w_hpo must be greater than 0 and at most 1, but was {threshold_HPO_observed_frequency}. Fix this before continuingf",
            self._logger.warning(errmsg)
            raise ValueError(errmsg)
        else:
            self._min_perc_patients_w_hpo = threshold_HPO_observed_frequency
        self._min_perc_patients_w_hpo = threshold_HPO_observed_frequency

    def specify_terms_strategy(self,  specified_term_set: typing.Iterable[typing.Union[str, hpotk.TermId]]):
        self.mtc_strategy = MTC_Strategy.SPECIFY_TERMS
        self._terms_to_test = specified_term_set

    def get_terms_to_test(self):
        return self._terms_to_test


def configure_cohort_analysis(cohort: Cohort,
                              hpo: hpotk.MinimalOntology,
                              protein_source: str = 'UNIPROT',
                              cache_dir: typing.Optional[str] = None,
                              config: typing.Optional[CohortAnalysisConfiguration] = None) -> CohortAnalysis:
    """
    Configure :class:`genophenocorr.analysis.CohortAnalysis` for given `cohort`.

    :param cohort: a :class:`genophenocorr.model.Cohort` to analyze
    :param hpo: a :class:`hpotk.MinimalOntology` with HPO to use in the analysis
    :param protein_source: the resource to retrieve protein annotations from if we cannot find the annotations locally.
     Choose from ``{'UNIPROT'}`` (just one fallback implementation is available at the moment).
    :param config: an optional :class:`CohortAnalysisConfiguration` to parameterize the analysis.
     The default parameters will be used if `None`.
    """
    if config is None:
        config = CohortAnalysisConfiguration()  # Use the default config
    if cache_dir is None:
        cache_dir = os.path.join(os.getcwd(), '.genophenocorr_cache')
    protein_metadata_service = _configure_protein_service(protein_source, cache_dir)

    # Phenotype filter defines how we select the HPO terms of interest.
    # For instance, we may choose to only investigate the terms
    # which annotate at least 3 cohort members, etc.
    phenotype_filter = SimplePhenotypeFilter(
        hpo,
        config.min_perc_patients_w_hpo,
    )

    mtc_filter: HpoMtcFilter
    if config.mtc_strategy == MTC_Strategy.HEURISTIC_SAMPLER:
        mtc_filter = HeuristicMtcFilter(hpo=hpo, hpo_term_frequency_filter = config._min_perc_patients_w_hpo)
    elif config.mtc_strategy == MTC_Strategy.SPECIFY_TERMS:
        validated_terms_to_test = _validate_terms_to_test(hpo, config.get_terms_to_test())
        mtc_filter = SpecifiedTermsMtcFilter(hpo=hpo, terms_to_test=validated_terms_to_test)
    elif config.mtc_strategy == MTC_Strategy.ALL_HPO_TERMS:
        mtc_filter = IdentityTermMtcFilter()
    else:
        raise ValueError(f"Did not recognize MtcStrategy {config.mtc_strategy}")

    # Choosing a simple Fisher's exact test for now.
    gp_analyzer = FisherExactAnalyzer(
        hpo_mtc_filter=mtc_filter,
        p_val_correction=config.pval_correction,
        mtc_alpha=config.mtc_alpha,
    )

    return GpCohortAnalysis(
        cohort=cohort,
        hpo=hpo,
        protein_service=protein_metadata_service,
        phenotype_filter=phenotype_filter,
        gp_analyzer=gp_analyzer,
        missing_implies_excluded=config.missing_implies_excluded,
        include_sv=config.include_sv,
    )

def _validate_terms_to_test(
        hpo: hpotk.MinimalOntology,
        terms_to_test: typing.Iterable[typing.Union[hpotk.TermId, str]],
    ) -> typing.Iterable[hpotk.TermId]:
    """
    Check that:
     * all terms to test are valid TermIds/CURIES,
     * the term IDs are in used HPO, and
     * there is at least one term to test
    """
    validated_terms_to_test = set()

    for term in terms_to_test:
        if isinstance(term, hpotk.TermId):
            pass
        if isinstance(term, str):
            term = hpotk.TermId.from_curie(term)
        else:
            raise ValueError(f'{term} is neither a TermId nor a CURIE `str`!')

        if term not in hpo:
            raise ValueError(f"HPO ID {term} not in HPO ontology {hpo.version}")
        validated_terms_to_test.add(term)
    if len(validated_terms_to_test) == 0:
        raise ValueError('Cannot run use {MTC_Strategy.SPECIFY_TERMS} with no HPO terms!')

    return validated_terms_to_test


def _configure_protein_service(protein_fallback: str, cache_dir) -> ProteinMetadataService:
    # (1) ProteinMetadataService
    # Setup fallback
    protein_fallback = _configure_fallback_protein_service(protein_fallback)
    # Setup protein metadata cache
    prot_cache_dir = os.path.join(cache_dir, 'protein_cache')
    os.makedirs(prot_cache_dir, exist_ok=True)
    prot_cache = ProteinAnnotationCache(prot_cache_dir)
    # Assemble the final protein metadata service
    protein_metadata_service = ProtCachingMetadataService(prot_cache, protein_fallback)
    return protein_metadata_service


def _configure_fallback_protein_service(protein_fallback: str) -> ProteinMetadataService:
    if protein_fallback == 'UNIPROT':
        fallback1 = UniprotProteinMetadataService()
    else:
        raise ValueError(f'Unknown protein fallback annotator type {protein_fallback}')
    return fallback1
