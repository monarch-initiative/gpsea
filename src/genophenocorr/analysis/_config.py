import logging
import os
import typing
import enum
import hpotk

from genophenocorr.model import Cohort
from genophenocorr.preprocessing import ProteinMetadataService, UniprotProteinMetadataService, ProteinAnnotationCache, \
    ProtCachingMetadataService
from ._api import CohortAnalysis
from ._mtc_filter import PhenotypeMtcFilter, UseAllTermsMtcFilter, SpecifiedTermsMtcFilter, HpoMtcFilter
from ._gp_analysis import FisherExactAnalyzer
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


class MtcStrategy(enum.Enum):
    """
    A strategy for mitigating the multiple testing correction (MTC) burden.
    """
    
    ALL_PHENOTYPE_TERMS = 0
    """
    All phenotype terms (HPO or disease IDs) will be tested.
    """

    SPECIFY_TERMS = 1
    """
    Only the manually provided HPO terms will be tested.
    """

    HPO_MTC = 2
    """
    Only HPO terms present in at least a certain fraction of patients will be tested.
    """


class CohortAnalysisConfiguration:
    """
    `CohortAnalysisConfiguration` is a value class for storing :class:`~genophenocorr.analysis.CohortAnalysis`
    configuration options.

    The class contains the default values upon creation and the configuration option values can be set as properties.

    If an invalid value option is passed to the property setter, a warning is logged and the previous value is retained.
    Therefore, it is *impossible* to mis-configure the analysis.

    Default values
    ^^^^^^^^^^^^^^

    ==============================  =======================  =========================================
        Option                        Type                    Default value
    ==============================  =======================  =========================================
     ``missing_implies_excluded``    `bool`                   `False`
     ``pval_correction``             `str`                    `bonferroni`
     ``min_n_patients_with_term``    `int`                    `2`
     ``mtc_alpha``                   `float`                  `0.05`
     ``include_sv``                  `bool`                   `False`
     ``mtc_strategy``                :class:`MtcStrategy`    :class:`MtcStrategy.ALL_PHENOTYPE_TERMS`
    ==============================  =======================  =========================================

    """

    def __init__(self):
        self._logger = logging.getLogger(__name__)
        self._missing_implies_excluded = False
        self._pval_correction = 'bonferroni'
        self._mtc_alpha = .05
        self._min_n_patients_with_term = 2
        self._include_sv = False
        self._mtc_strategy = MtcStrategy.ALL_PHENOTYPE_TERMS
        self._terms_to_test = None  # # only relevant for SPECIFIED_TERMS strategy
        self._min_patients_w_hpo = None  # # only relevant for HPO_MTC strategy

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
        Set the method for p value correction. 
        See Statsmodels' `documentation <https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html>`_ 
        for the acceptable values.

        :param pval_correction: a `str` with the name of the desired multiple testing correction method or `None`
          if no MTC should be applied.
        """
        if pval_correction in P_VAL_OPTIONS:
            self._pval_correction = pval_correction
        else:
            self._logger.warning('Ignoring invalid `pval_correction` value %s. Using %s correction.', pval_correction,
                                 self._pval_correction)

    @property
    def min_n_patients_with_term(self) -> int:
        """
        Get the minimum number of patients that must be annotated with an HPO term
        for including the term in the analysis.
        """
        return self._min_n_patients_with_term

    @min_n_patients_with_term.setter
    def min_n_patients_with_term(self, value: int):
        if isinstance(value, int) and value >= 0:
            self._min_n_patients_with_term = value
        else:
            self._logger.warning(
                'Ignoring invalid `min_n_patients_with_term` value %s. Using %s', 
                value, self._min_n_patients_with_term,
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
    def mtc_strategy(self) -> MtcStrategy:
        """
        Get the MTC filtering strategy to be used.
        """
        return self._mtc_strategy

    def all_terms_strategy(self):
        """
        Test all phenotype terms.

        See :ref:`use-all-terms-strategy` for more info.
        """
        self._mtc_strategy = MtcStrategy.ALL_PHENOTYPE_TERMS
        self._min_patients_w_hpo = None
        self._terms_to_test = None

    def hpo_mtc_strategy(
        self,
        min_patients_w_hpo: float = 0.2,
    ):
        """
        Only test the HPO terms that pass all rules of the HPO filter strategy.

        See :ref:`hpo-mtc-filter-strategy` section for more info on the rules.

        :param threshold_HPO_observed_frequency: a float in range :math:`(0, 1]` to represent
          the minimum fraction of patients for an HPO term to be included.
        """
        if not isinstance(min_patients_w_hpo, float):
            raise ValueError(f'`min_patients_w_hpo` is not a `float`: {min_patients_w_hpo}')
        if not 0 < min_patients_w_hpo <= 1:
            raise ValueError(f'`min_patients_w_hpo` must be in range (0, 1] but was {min_patients_w_hpo}')
       
        self._mtc_strategy = MtcStrategy.HPO_MTC
        self._min_patients_w_hpo = min_patients_w_hpo
        self._terms_to_test = None

    def specify_terms_strategy(
        self,
        terms_to_test: typing.Iterable[typing.Union[str, hpotk.TermId]],
    ):
        """
        Mitigate the MTC burden by only testing the specified HPO terms.

        The HPO terms are validated before running the analysis,
        to point out invalid CURIE (e.g. `WHATEVER`) values,
        or HPO term IDs that are not in the currently used HPO.

        Calling this method will clear any previously specified terms.

        See :ref:`specify-terms-strategy` for more info.

        :param terms_to_test: an iterable with CURIEs (e.g. `HP:0001250`)
          or :class:`hpotk.TermId` instances representing the terms to test.
        """
        self._mtc_strategy = MtcStrategy.SPECIFY_TERMS
        self._min_patients_w_hpo = None
        self._terms_to_test = tuple(terms_to_test)

    @property
    def terms_to_test(self) -> typing.Optional[typing.Iterable[typing.Union[str, hpotk.TermId]]]:
        """
        Get the ids of the terms to be tested in `specify_terms_strategy`
        or `None` if :class:`MtcStrategy.SPECIFY_TERMS` will *not* be used.
        """
        return self._terms_to_test
    
    @property
    def min_patients_w_hpo(self) -> typing.Optional[float]:
        """
        Get the minimum fraction of patients needed to be annotated with an HPO term
        to justify being tested or `None` if :class:`MtcStrategy.HPO_MTC` will *not* be used.
        """
        return self._min_patients_w_hpo


def configure_cohort_analysis(
    cohort: Cohort,
    hpo: hpotk.MinimalOntology,
    protein_source: str = 'UNIPROT',
    cache_dir: typing.Optional[str] = None,
    config: typing.Optional[CohortAnalysisConfiguration] = None,
) -> CohortAnalysis:
    """
    Configure :class:`~genophenocorr.analysis.CohortAnalysis` for given `cohort`.

    :param cohort: a :class:`~genophenocorr.model.Cohort` to analyze
    :param hpo: a :class:`~hpotk.MinimalOntology` with HPO to use in the analysis
    :param protein_source: the resource to retrieve protein annotations from if we cannot find the annotations locally.
     Choose from ``{'UNIPROT'}`` (just one fallback implementation is available at the moment).
    :param config: an optional :class:`CohortAnalysisConfiguration` to parameterize the analysis.
     The default parameters will be used if `None`.
    """
    if config is None:
        config = CohortAnalysisConfiguration()  # Use the default config
    cache_dir = _configure_cache_dir(cache_dir)
    protein_metadata_service = configure_default_protein_metadata_service(protein_source, cache_dir)

    mtc_filter: PhenotypeMtcFilter
    if config.mtc_strategy == MtcStrategy.HPO_MTC:
        assert config.min_patients_w_hpo is not None, '`min_patients_w_hpo` must be set if using `HPO_MTC` strategy'
        mtc_filter = HpoMtcFilter.default_filter(
            hpo=hpo,
            term_frequency_threshold=config.min_patients_w_hpo,
        )
    elif config.mtc_strategy == MtcStrategy.SPECIFY_TERMS:
        assert config.terms_to_test is not None, '`terms_to_test` must be set if using `SPECIFY_TERMS` strategy'
        validated_terms_to_test = _validate_terms_to_test(hpo, config.terms_to_test)
        mtc_filter = SpecifiedTermsMtcFilter(hpo=hpo, terms_to_test=validated_terms_to_test)
    elif config.mtc_strategy == MtcStrategy.ALL_PHENOTYPE_TERMS:
        mtc_filter = UseAllTermsMtcFilter()
    else:
        raise ValueError(f"Did not recognize MtcStrategy {config.mtc_strategy}")

    # Choosing a simple Fisher's exact test for now.
    gp_analyzer = FisherExactAnalyzer(
        mtc_filter=mtc_filter,
        p_val_correction=config.pval_correction,
        mtc_alpha=config.mtc_alpha,
    )

    return GpCohortAnalysis(
        cohort=cohort,
        hpo=hpo,
        protein_service=protein_metadata_service,
        gp_analyzer=gp_analyzer,
        missing_implies_excluded=config.missing_implies_excluded,
        min_n_of_patients_with_term=config.min_n_patients_with_term,
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
        elif isinstance(term, str):
            term = hpotk.TermId.from_curie(term)
        else:
            raise ValueError(f'{term} is neither a TermId nor a CURIE `str`!')

        if term not in hpo:
            raise ValueError(f"HPO ID {term} not in HPO ontology {hpo.version}")
        validated_terms_to_test.add(term)
    if len(validated_terms_to_test) == 0:
        raise ValueError('Cannot run use {MTC_Strategy.SPECIFY_TERMS} with no HPO terms!')

    return validated_terms_to_test


def configure_default_protein_metadata_service(
    protein_source: str = 'UNIPROT',
    cache_dir: typing.Optional[str] = None,
) -> ProteinMetadataService:
    """
    Create default protein metadata service that will cache the protein metadata 
    in current working directory under `.genophenocorr_cache/protein_cache` 
    and reach out to UNIPROT REST API if a cache entry is missing.
    """
    cache_dir = _configure_cache_dir(cache_dir)
    return _configure_protein_service(protein_fallback=protein_source, cache_dir=cache_dir)


def _configure_protein_service(
        protein_fallback: str,
        cache_dir: str,
) -> ProteinMetadataService:
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


def _configure_cache_dir(cache_dir: typing.Optional[str]) -> str:
    if cache_dir is None:
        cache_dir = os.path.join(os.getcwd(), '.genophenocorr_cache')
    return cache_dir
