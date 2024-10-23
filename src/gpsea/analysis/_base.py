import abc
import typing


from .predicate.genotype import GenotypePolyPredicate


class Statistic(metaclass=abc.ABCMeta):
    """
    Mixin for classes that are used to compute a nominal p value for a genotype-phenotype association.
    """
    
    def __init__(
        self,
        name: str,
    ):
        self._name = name

    @property
    def name(self) -> str:
        """
        Get the name of the statistic (e.g. `Fisher Exact Test`, `Logrank test`).
        """
        return self._name
    
    def __eq__(self, value: object) -> bool:
        if isinstance(value, Statistic):
            return self._name == value._name
        return NotImplemented
    
    def __hash__(self) -> int:
        return hash((self._name,))


class AnalysisResult(metaclass=abc.ABCMeta):
    """
    `AnalysisResult` includes the common parts of results of all analyses.
    """

    def __init__(
        self,
        gt_predicate: GenotypePolyPredicate,
        statistic: Statistic,
        mtc_correction: typing.Optional[str]
    ):
        assert isinstance(gt_predicate, GenotypePolyPredicate)
        self._gt_predicate = gt_predicate

        assert isinstance(statistic, Statistic)
        self._statistic = statistic

        if mtc_correction is not None:
            assert isinstance(mtc_correction, str)
        self._mtc_correction = mtc_correction
    
    @property
    def gt_predicate(self) -> GenotypePolyPredicate:
        """
        Get the genotype predicate used in the survival analysis that produced this result.
        """
        return self._gt_predicate
    
    @property
    def statistic(self) -> Statistic:
        """
        Get the statistic which computed the (nominal) p values for this result.
        """
        return self._statistic
    
    @property
    def mtc_correction(self) -> typing.Optional[str]:
        """
        Get name/code of the used multiple testing correction
        (e.g. `fdr_bh` for Benjamini-Hochberg) or `None` if no correction was applied.
        """
        return self._mtc_correction

    # test, p values, corrected values (optional), correction

    def __eq__(self, value: object) -> bool:
        return isinstance(value, AnalysisResult) \
            and self._gt_predicate == value._gt_predicate \
            and self._statistic == value._statistic \
            and self._mtc_correction == value._mtc_correction
    
    def __hash__(self) -> int:
        return hash((
            self._gt_predicate,
            self._statistic,
            self._mtc_correction,
        ))
