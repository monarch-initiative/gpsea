import abc
import typing

from collections import namedtuple

#from genophenocorr.cohort import Cohort
from genophenocorr.patient import Patient


class SimplePredicate(metaclass=abc.ABCMeta):
    """
    Predicate for testing if a :class:`Patient` belongs to a group.
    """
    def __init__(self) -> None:
        pass
        

    @abc.abstractmethod
    def test(self, patient: Patient) -> bool:
        pass


PatientCategory = namedtuple('PatientCategory', field_names=['cat_id', 'name'])


class PolyPredicate(metaclass=abc.ABCMeta):
    """
    Predicate for categorizing a :class:`Patient` into one of discrete groups.
    """

    @property
    @abc.abstractmethod
    def categories(self) -> typing.Sequence[PatientCategory]:
        """
        Get a sequence of all categories which the `PolyPredicate` can produce.
        """
        pass

    @abc.abstractmethod
    def test(self, patient: Patient) -> typing.Optional[PatientCategory]:
        """
        Assign patient into a category.

        Return `None` if the patient cannot be assigned into any meaningful category.
        """
        pass
