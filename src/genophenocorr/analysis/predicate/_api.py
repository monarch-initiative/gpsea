import abc
import typing

import hpotk.util

from genophenocorr.model import Patient


class PatientCategory:
    """
    `PatientCategory` represents one of several exclusive discrete groups.

    Patient category has :attr:`cat_id`, a unique numeric identifier of the group,
    :attr:`name` with human-readable group name, and :attr:`description` with an optional verbose description.
    """

    def __init__(self, cat_id: int,
                 name: str,
                 description: typing.Optional[str] = None):
        self._cat_id = hpotk.util.validate_instance(cat_id, int, 'cat_id')
        self._name = hpotk.util.validate_instance(name, str, 'name')
        self._description = hpotk.util.validate_optional_instance(description, str, 'description')

    @property
    def cat_id(self) -> int:
        """
        Get an `int` with the unique numeric identifier of the group.
        """
        return self._cat_id

    @property
    def name(self) -> str:
        """
        Get a `str` with a human-readable name of the group.
        """
        return self._name

    @property
    def description(self) -> typing.Optional[str]:
        """
        Get a `str` with an optional detailed group description.
        """
        return self._description

    def __repr__(self) -> str:
        return f"PatientCategory(cat_id={self.cat_id}, " \
               f"name={self.name}, " \
               f"description={self.description})"

    def __str__(self) -> str:
        return self.name

    def __eq__(self, other) -> bool:
        return isinstance(other, PatientCategory) \
            and self.cat_id == other.cat_id \
            and self.name == other.name \
            and self.description == other.description

    def __hash__(self) -> int:
        return hash((self.cat_id, self.name, self.description))


class PolyPredicate(metaclass=abc.ABCMeta):
    """
    `PolyPredicate` bins a :class:`genophenocorr.model.Patient` into one of several discrete groups represented
    by :class:`PatientCategory`.

    The groups must be *exclusive* - the patient can be binned into one and only one group,
    and *exhaustive* - the groups must cover all possible scenarios.

    However, if the patient cannot be assigned into any meaningful category, `None` can be returned.

    Note, that `PolyPredicate` must be used on a happy path - testing a patient must inherently make sense.
    Predicate will *not* check if, for instance, the patient variants are compatible with a certain mode of inheritance.
    """

    @staticmethod
    @abc.abstractmethod
    def get_categories() -> typing.Sequence[PatientCategory]:
        """
        Get a sequence of all categories which the `PolyPredicate` can produce.
        """
        pass

    @abc.abstractmethod
    def get_question(self) -> str:
        """
        Prepare a `str` with the question the predicate can answer.
        """
        pass

    @abc.abstractmethod
    def test(self, patient: Patient) -> typing.Optional[PatientCategory]:
        """
        Assign a `patient` into a category.

        Return `None` if the patient cannot be assigned into any meaningful category.
        """
        pass

    @classmethod
    def n_categories(cls) -> int:
        """
        Get the number of categories the predicate can produce.
        """
        return len(cls.get_categories())

    @staticmethod
    def _check_patient(patient: Patient):
        """
        Check if the `patient` meets the predicate requirements.
        """
        if not isinstance(patient, Patient):
            raise ValueError(f"patient must be type Patient but was type {type(patient)}")


class BooleanPredicate(PolyPredicate, metaclass=abc.ABCMeta):
    """
    `BooleanPredicate` tests if a :class:`genophenocorr.model.Patient` belongs to a group and returns a boolean binning.
    """

    NO = PatientCategory(0, 'No', 'The patient does not belong to the group.')
    """
    Category for a patient who does *not* belong to the tested group.
    """

    YES = PatientCategory(1, 'Yes', 'The patient belongs to the group.')
    """
    Category for a patient who *belongs* to the tested group.
    """

    @staticmethod
    def get_categories() -> typing.Sequence[PatientCategory]:
        """
        The predicate bins a patient into :class:`BooleanPredicate.NO` or :class:`BooleanPredicate.YES` category.
        """
        return BooleanPredicate.NO, BooleanPredicate.YES
