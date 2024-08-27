import abc

from gpsea.model import Patient


class PhenotypeScorer(metaclass=abc.ABCMeta):
    """
    `PhenotypeScorer` assigns the patient with a phenotype score.
    """
    
    def score(self, patient: Patient) -> float:
        """
        Compute the score for the `patient`.
        """
        pass
