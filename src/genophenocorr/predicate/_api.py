import abc
import typing

import hpotk.util

from genophenocorr.patient import Patient

T = typing.TypeVar('T')

class SimplePredicate(metaclass=abc.ABCMeta):
    """
    Predicate for testing if a :class:`Patient` belongs to a group.
    """
    def __init__(self) -> None:
        pass
        

    @abc.abstractmethod
    def test(self, patient: Patient) -> bool:
        pass


class PatientCategory:

    def __init__(self, cat_id: int,
                 name: str,
                 description: typing.Optional[str] = None):
        self._cat_id = hpotk.util.validate_instance(cat_id, int, 'cat_id')
        self._name = hpotk.util.validate_instance(name, str, 'name')
        self._description = hpotk.util.validate_optional_instance(description, str, 'description')

    @property
    def cat_id(self) -> int:
        return self._cat_id

    @property
    def name(self) -> str:
        return self._name

    @property
    def description(self) -> typing.Optional[str]:
        return self._description

    def __str__(self) -> str:
        return f"PatientCategory(cat_id={self.cat_id}, " \
               f"name={self.name}, " \
               f"description={self.description})"

    def __repr__(self) -> str:
        return str(self)

    def __eq__(self, other) -> bool:
        return isinstance(other, PatientCategory) \
            and self.cat_id == other.cat_id \
            and self.name == other.name \
            and self.description == other.description

    def __hash__(self) -> int:
        return hash((self.cat_id, self.name, self.description))


class PolyPredicate(typing.Generic[T], metaclass=abc.ABCMeta):
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
    def test(self, patient: Patient, query: T) -> typing.Optional[PatientCategory]:
        """
        Assign patient into a category.

        Return `None` if the patient cannot be assigned into any meaningful category.
        """
        pass
