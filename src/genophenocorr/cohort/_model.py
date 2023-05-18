import abc
from ._cohort_data import Cohort


class CohortCreator(metaclass=abc.ABCMeta):
    @abc.abstractmethod
    def create_cohort(self, item) -> Cohort:
        pass
