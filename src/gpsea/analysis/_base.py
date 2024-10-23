import abc


from .predicate.genotype import GenotypePolyPredicate


class AnalysisResult(metaclass=abc.ABCMeta):
    """
    `AnalysisResult` includes the common parts of results of all analyses.
    """

    def __init__(
        self,
        gt_predicate: GenotypePolyPredicate,
    ):
        assert isinstance(gt_predicate, GenotypePolyPredicate)
        self._gt_predicate = gt_predicate
    
    @property
    def gt_predicate(self) -> GenotypePolyPredicate:
        """
        Get the genotype predicate used in the survival analysis that produced this result.
        """
        return self._gt_predicate

    def __eq__(self, value: object) -> bool:
        return isinstance(value, AnalysisResult) \
            and self._gt_predicate == value._gt_predicate
    
    def __hash__(self) -> int:
        return hash((self._gt_predicate,))
