import abc
import os
import typing

from collections import Counter

import hpotk
import hpotk.util

from gpsea.model import Patient

from .._partition import Partitioning


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
        return self._name

    def __eq__(self, other) -> bool:
        return isinstance(other, PatientCategory) \
            and self.cat_id == other.cat_id \
            and self.name == other.name \
            and self.description == other.description

    def __hash__(self) -> int:
        return hash((self.cat_id, self.name, self.description))


class Categorization:
    """
    `Categorization` represents one of discrete group a :class:`~gpsea.model.Patient` can be assigned into.
    """

    @staticmethod
    def from_raw_parts(
        cat_id: int,
        name: str,
        description: typing.Optional[str] = None,
    ):
        """
        Create `Categorization` from the `cat_id` identifier, `name`, and an optional `description`.
        """
        return Categorization(
            category=PatientCategory(
                cat_id=cat_id,
                name=name,
                description=description,
            )
        )

    def __init__(
        self,
        category: PatientCategory,
    ):
        self._category = hpotk.util.validate_instance(category, PatientCategory, 'category')

    @property
    def category(self) -> PatientCategory:
        return self._category

    def __eq__(self, other):
        return isinstance(other, Categorization) and self._category == other._category

    def __hash__(self):
        return hash((self._category,))

    def __repr__(self):
        return f'Categorization(category={self._category})'

    def __str__(self):
        return repr(self)


C = typing.TypeVar('C', bound=Categorization)
"""
A generic bound for types that extend :class:`Categorization`.
"""


class PolyPredicate(typing.Generic[C], Partitioning, metaclass=abc.ABCMeta):
    """
    `PolyPredicate` partitions a :class:`~gpsea.model.Patient` into one of several discrete groups
    represented by a :class:`~gpsea.analysis.predicate.Categorization`.

    The groups must be *exclusive* - the patient can be binned into one and only one group,
    and *exhaustive* - the groups must cover all possible scenarios.

    However, if the patient cannot be assigned into any meaningful category, `None` can be returned.
    As a rule of thumb, returning `None` will exclude the patient from the analysis.
    """

    @abc.abstractmethod
    def get_categorizations(self) -> typing.Sequence[C]:
        """
        Get a sequence of all categories which the predicate can produce.
        """
        pass

    def get_categories(self) -> typing.Iterator[PatientCategory]:
        """
        Get an iterator with :class:`PatientCategory` instances that the predicate can produce.
        """
        return (c.category for c in self.get_categorizations())

    @property
    def group_labels(self) -> typing.Collection[str]:
        """
        Get a collection with names of the :class:`PatientCategory` items that the predicate can produce.
        """
        return tuple(cat.name for cat in self.get_categories())

    def summarize_groups(self) -> str:
        cat_names = ', '.join(self.group_labels)
        return f"{self.variable_name}: {cat_names}"

    def summarize(
        self,
        out: typing.TextIO,
    ):
        """
        Summarize the predicate into the `out` handle.

        The summary includes the name, summary, and the groups the predicate can assign individuals into.
        """
        Partitioning.summarize(self, out)
        out.write(self.summarize_groups())
        out.write(os.linesep)

    def n_categorizations(self) -> int:
        """
        Get the number of categorizations the predicate can produce.
        """
        return len(self.get_categorizations())

    def get_category(
        self,
        cat_id: int,
    ) -> PatientCategory:
        """
        Get the category name for a :class:`PatientCategory.cat_id`.

        :param cat_id: an `int` with the id.
        :raises: ValueError if there is no such category was defined.
        """
        for ctg in self.get_categories():
            if ctg.cat_id == cat_id:
                return ctg
        raise ValueError(f'No category for {cat_id} was found')

    def get_category_name(
        self,
        cat_id: int,
    ) -> str:
        """
        Get the category name for a :class:`PatientCategory.cat_id`.

        :param cat_id: an `int` with the id.
        :raises: ValueError if there is no such category was defined.
        """
        return self.get_category(cat_id).name
    
    @staticmethod
    def _check_categorizations(
        categorizations: typing.Sequence[Categorization],
    ) -> typing.Sequence[str]:
        """
        Check that the categorizations meet the requirements.

        The requirements include:

        * the `cat_id` must be unique within the predicate
        """
        issues = []
        # There must be at least one category

        # `cat_id` must be unique for a predicate!
        cat_id_counts = Counter()
        for c in categorizations:
            cat_id_counts[c.category.cat_id] += 1

        for cat_id, count in cat_id_counts.items():
            if count > 1:
                issues.append(f'`cat_id` {cat_id} is present {count}>1 times')
        
        return issues

    @abc.abstractmethod
    def test(self, patient: Patient) -> typing.Optional[C]:
        """
        Assign a `patient` into a categorization.

        Return `None` if the patient cannot be assigned into any meaningful category.
        """
        pass
