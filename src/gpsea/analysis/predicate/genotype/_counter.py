from gpsea.model import Patient, Genotype
from ._api import VariantPredicate


class AlleleCounter:
    """
    `AlleleCounter` counts the number of alleles of all variants that pass the selection with a given `predicate`.

    :param predicate: a :class:`VariantPredicate` for selecting the target variants.
    """
    # TODO: this class should probably be an implementation detail, 
    #   and not a public member of the package.

    def __init__(
        self,
        predicate: VariantPredicate,
    ):
        self._predicate = predicate

    def get_question(self) -> str:
        """
        Get the question tested by the predicate.

        Returns:
            str: the question tested by the predicate
        """
        return self._predicate.get_question()

    def count(
        self,
        patient: Patient,
    ) -> int:
        """
        Count the number of alleles of all variants that pass the predicate.
        Args:
            patient: the patient to test

        Returns:
            int: the count of the passing alleles
        """
        count = 0

        for var in patient.variants:
            if self._predicate.test(var):
                genotype = var.genotypes.for_sample(patient.labels)
                if genotype == Genotype.HOMOZYGOUS_ALTERNATE:
                    count += 2
                elif (
                    genotype == Genotype.HETEROZYGOUS or genotype == Genotype.HEMIZYGOUS
                ):
                    count += 1
        
        return count

    def __eq__(self, value: object) -> bool:
        return isinstance(value, AlleleCounter) and self._predicate == value._predicate
    
    def __hash__(self) -> int:
        return hash((self._predicate,))
    
    def __str__(self) -> str:
        return f'AlleleCounter(predicate={self._predicate})'
    
    def __repr__(self) -> str:
        return str(self)
