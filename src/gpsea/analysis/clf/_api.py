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
    `PatientCategory` represents one of several exclusive discrete classes.

    Patient class has :attr:`cat_id`, a unique numeric identifier of the class,
    :attr:`name` with human-readable class name, and :attr:`description` with an optional verbose description.
    """

    def __init__(
        self,
        cat_id: int,
        name: str,
        description: typing.Optional[str] = None,
    ):
        self._cat_id = hpotk.util.validate_instance(cat_id, int, "cat_id")
        self._name = hpotk.util.validate_instance(name, str, "name")
        self._description = hpotk.util.validate_optional_instance(
            description, str, "description"
        )

    @property
    def cat_id(self) -> int:
        """
        Get an `int` with the unique numeric identifier of the class.
        """
        return self._cat_id

    @property
    def name(self) -> str:
        """
        Get a `str` with a human-readable name of the class.
        """
        return self._name

    @property
    def description(self) -> typing.Optional[str]:
        """
        Get a `str` with an optional detailed class description.
        """
        return self._description

    def __repr__(self) -> str:
        return (
            f"PatientCategory(cat_id={self.cat_id}, "
            f"name={self.name}, "
            f"description={self.description})"
        )

    def __str__(self) -> str:
        return self._name

    def __eq__(self, other) -> bool:
        return (
            isinstance(other, PatientCategory)
            and self.cat_id == other.cat_id
            and self.name == other.name
            and self.description == other.description
        )

    def __hash__(self) -> int:
        return hash((self.cat_id, self.name, self.description))


class Categorization:
    """
    `Categorization` represents one of discrete classes a :class:`~gpsea.model.Patient` can be assigned into.
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
        self._category = hpotk.util.validate_instance(
            category, PatientCategory, "category"
        )

    @property
    def category(self) -> PatientCategory:
        return self._category

    def __eq__(self, other):
        return isinstance(other, Categorization) and self._category == other._category

    def __hash__(self):
        return hash((self._category,))

    def __repr__(self):
        return f"Categorization(category={self._category})"

    def __str__(self):
        return repr(self)


C = typing.TypeVar("C", bound=Categorization)
"""
A generic bound for types that extend :class:`Categorization`.
"""


class Classifier(typing.Generic[C], Partitioning, metaclass=abc.ABCMeta):
    """
    `Classifier` partitions a :class:`~gpsea.model.Patient` into one of several discrete classes
    represented by a :class:`~gpsea.analysis.clf.Categorization`.

    The classes must be *exclusive* - the individual can be binned into one and only one class,
    and *exhaustive* - the classes must cover all possible scenarios.

    However, if the individual cannot be assigned into any meaningful class, `None` can be returned.
    As a rule of thumb, returning `None` will exclude the individual from the analysis.
    """

    @abc.abstractmethod
    def get_categorizations(self) -> typing.Sequence[C]:
        """
        Get a sequence of all categories which the classifier can produce.
        """
        pass

    def get_categories(self) -> typing.Iterator[PatientCategory]:
        """
        Get an iterator with :class:`PatientCategory` instances that the classifier can produce.
        """
        return (c.category for c in self.get_categorizations())

    @property
    def class_labels(self) -> typing.Collection[str]:
        """
        Get a collection with names of the :class:`PatientCategory` items that the classifier can produce.
        """
        return tuple(cat.name for cat in self.get_categories())

    def summarize_classes(self) -> str:
        cat_names = ", ".join(self.class_labels)
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
        out.write(self.summarize_classes())
        out.write(os.linesep)

    def n_categorizations(self) -> int:
        """
        Get the number of categorizations the classifier can produce.
        """
        return len(self.get_categorizations())

    def get_category(
        self,
        cat_id: int,
    ) -> PatientCategory:
        """
        Get the category name for a :attr:`PatientCategory.cat_id`.

        :param cat_id: an `int` with the id.
        :raises: ValueError if there is no such category was defined.
        """
        for ctg in self.get_categories():
            if ctg.cat_id == cat_id:
                return ctg
        raise ValueError(f"No category for {cat_id} was found")

    def get_category_name(
        self,
        cat_id: int,
    ) -> str:
        """
        Get the category name for a :attr:`PatientCategory.cat_id`.

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
                issues.append(f"`cat_id` {cat_id} is present {count}>1 times")

        return issues

    @abc.abstractmethod
    def test(self, individual: Patient) -> typing.Optional[C]:
        """
        Assign an `individual` into a class.

        Return `None` if the individual cannot be assigned into any meaningful class.
        """
        pass


class GenotypeClassifier(Classifier[Categorization], metaclass=abc.ABCMeta):
    """
    `GenotypeClassifier` is a base class for all types
    that assign an individual into a group based on the genotype.
    """

    pass


P = typing.TypeVar("P")
"""
Phenotype entity of interest, such as :class:`~hpotk.model.TermId`, representing an HPO term or an OMIM/MONDO term.

However, phenotype entity can be anything as long as it is :class:`~typing.Hashable` and comparable
(have `__eq__` and `__lt__` magic methods).
"""

YES = PatientCategory(1, 'Yes', 'The patient belongs to the group.')
"""
Category for a patient who *belongs* to the tested group.
"""

NO = PatientCategory(0, 'No', 'The patient does not belong to the group.')
"""
Category for a patient who does *not* belong to the tested group.
"""


class PhenotypeCategorization(typing.Generic[P], Categorization):
    """
    On top of the attributes of the `Categorization`, `PhenotypeCategorization` keeps track of the target phenotype `P`.
    """

    def __init__(
        self,
        category: PatientCategory,
        phenotype: P,
    ):
        super().__init__(category)
        self._phenotype = phenotype

    @property
    def phenotype(self) -> P:
        return self._phenotype

    def __repr__(self) -> str:
        return (
            "PhenotypeCategorization("
            f"category={self._category}, "
            f"phenotype={self._phenotype})"
        )

    def __str__(self) -> str:
        return repr(self)

    def __eq__(self, other) -> bool:
        return isinstance(other, PhenotypeCategorization) \
            and self._category == other._category \
            and self._phenotype == other._phenotype

    def __hash__(self) -> int:
        return hash((self._category, self._phenotype))


class PhenotypeClassifier(
    typing.Generic[P], Classifier[PhenotypeCategorization[P]], 
    metaclass=abc.ABCMeta,
):
    """
    Phenotype classifier assigns an individual into a class `P` based on the phenotype.

    The class `P` can be a :class:`~hpotk.model.TermId` representing an HPO term or an OMIM/MONDO term.

    Only one class can be investigated, and :attr:`phenotype` returns the investigated phenotype
    (e.g. *Arachnodactyly* `HP:0001166`).

    As another hallmark of this predicate, one of the categorizations must correspond to the group of patients
    who exibit the investigated phenotype. The categorization is provided
    via :attr:`present_phenotype_categorization` property.
    """

    @property
    @abc.abstractmethod
    def phenotype(self) -> P:
        """
        Get the phenotype entity of interest.
        """
        pass

    @property
    @abc.abstractmethod
    def present_phenotype_categorization(self) -> PhenotypeCategorization[P]:
        """
        Get the categorization which represents the group of the patients who exibit the investigated phenotype.
        """
        pass

    @property
    def present_phenotype_category(self) -> PatientCategory:
        """
        Get the patient category that correspond to the group of the patients who exibit the investigated phenotype.
        """
        return self.present_phenotype_categorization.category
