import abc
import typing

import hpotk

from gpsea.model import Patient

from .._api import PolyPredicate, PatientCategory, PatientCategories, Categorization

P = typing.TypeVar("P")
"""
Phenotype entity of interest, such as :class:`~hpotk.model.TermId`, representing an HPO term or an OMIM/MONDO term.

However, phenotype entity can be anything as long as it is :class:`~typing.Hashable` and comparable
(have `__eq__` and `__lt__` magic methods).
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


class PhenotypePolyPredicate(
    typing.Generic[P], PolyPredicate[PhenotypeCategorization[P]], metaclass=abc.ABCMeta
):
    """
    `PhenotypePolyPredicate` investigates a patient in context of one or more phenotype categories `P`.

    Each phenotype category `P` can be a :class:`~hpotk.model.TermId` representing an HPO term or an OMIM/MONDO term.

    Only one category can be investigated, and :attr:`phenotype` returns the investigated phenotype
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


class PropagatingPhenotypePredicate(PhenotypePolyPredicate[hpotk.TermId]):
    """
    `PropagatingPhenotypePredicate` tests if a patient is annotated with an HPO term.

    Note, `query` must be a term of the provided `hpo`!

    :param hpo: HPO object
    :param query: the HPO term to test
    :param missing_implies_phenotype_excluded: `True` if lack of an explicit annotation implies term's absence`.
    """

    def __init__(
        self,
        hpo: hpotk.MinimalOntology,
        query: hpotk.TermId,
        missing_implies_phenotype_excluded: bool = False,
    ):
        self._hpo = hpotk.util.validate_instance(hpo, hpotk.MinimalOntology, "hpo")
        self._query = hpotk.util.validate_instance(
            query, hpotk.TermId, "phenotypic_feature"
        )
        self._query_label = self._hpo.get_term_name(query)
        assert self._query_label is not None, f"Query {query} is in HPO"
        self._missing_implies_phenotype_excluded = hpotk.util.validate_instance(
            missing_implies_phenotype_excluded,
            bool,
            "missing_implies_phenotype_excluded",
        )
        self._phenotype_observed = PhenotypeCategorization(
            category=PatientCategories.YES,
            phenotype=self._query,
        )
        self._phenotype_excluded = PhenotypeCategorization(
            category=PatientCategories.NO,
            phenotype=self._query,
        )
        # Some tests depend on the order of `self._categorizations`.
        self._categorizations = (self._phenotype_observed, self._phenotype_excluded)

    def get_question_base(self) -> str:
        return f"Is {self._query_label} present in the patient"

    @property
    def phenotype(self) -> hpotk.TermId:
        return self._query

    @property
    def present_phenotype_categorization(self) -> PhenotypeCategorization[hpotk.TermId]:
        return self._phenotype_observed

    def get_categorizations(
        self,
    ) -> typing.Sequence[PhenotypeCategorization[hpotk.TermId]]:
        return self._categorizations

    def test(
        self, patient: Patient
    ) -> typing.Optional[PhenotypeCategorization[hpotk.TermId]]:
        """An HPO TermID is given when initializing the class.
        Given a Patient class, this function tests whether the patient has the
        given phenotype.

        Args:
            patient (Patient): A Patient class representing a patient.

        Returns:
            typing.Optional[PhenotypeCategorization[P]]: PhenotypeCategorization,
                                                        either "YES" if the phenotype
                                                        is listed and is not excluded, or
                                                        "NO" if the phenotype is listed and excluded,
                                                        otherwise will return None.
                                                        Unless _missing_implies_phenotype_excluded is True, then
                                                        will return "NO" if the phenotype is listed and excluded
                                                        or not listed.
        """
        self._check_patient(patient)

        if len(patient.phenotypes) == 0:
            return None

        for phenotype in patient.phenotypes:
            if phenotype.is_present:
                if any(
                    self._query == anc
                    for anc in self._hpo.graph.get_ancestors(
                        phenotype, include_source=True
                    )
                ):
                    return self._phenotype_observed
            else:
                if self._missing_implies_phenotype_excluded:
                    return self._phenotype_excluded
                else:
                    if any(
                        self._query == desc
                        for desc in self._hpo.graph.get_descendants(
                            phenotype, include_source=True
                        )
                    ):
                        return self._phenotype_excluded

        return None

    def __repr__(self):
        return f"PropagatingPhenotypeBooleanPredicate(query={self._query})"


class DiseasePresencePredicate(PhenotypePolyPredicate[hpotk.TermId]):
    """
    `DiseasePresencePredicate` tests if the patient was diagnosed with a disease.

    The predicate tests if the patient's diseases include a disease ID formatted as a :class:`~hpotk.model.TermId`.

    :param disease_id_query: the Disease ID to test
    """

    def __init__(self, disease_id_query: hpotk.TermId):
        self._query = hpotk.util.validate_instance(
            disease_id_query, hpotk.TermId, "disease_id_query"
        )

        self._diagnosis_present = PhenotypeCategorization(
            category=PatientCategories.YES,
            phenotype=disease_id_query,
        )
        self._diagnosis_excluded = PhenotypeCategorization(
            category=PatientCategories.NO,
            phenotype=disease_id_query,
        )

    def get_question_base(self) -> str:
        return f"Was {self._query} diagnosed in the patient"

    @property
    def phenotype(self) -> hpotk.TermId:
        return self._query

    @property
    def present_phenotype_categorization(self) -> PhenotypeCategorization[hpotk.TermId]:
        return self._diagnosis_present

    def get_categorizations(
        self,
    ) -> typing.Sequence[PhenotypeCategorization[hpotk.TermId]]:
        return self._diagnosis_present, self._diagnosis_excluded

    def test(
        self, patient: Patient
    ) -> typing.Optional[PhenotypeCategorization[hpotk.TermId]]:
        """
        Test if the `patient` was diagnosed with a disease.

        Args:
            patient (Patient): a `patient` instance to be tested.

        Returns:
            typing.Optional[PhenotypeCategorization[P]]: PhenotypeCategorization,
                                                        either "YES" if the phenotype
                                                        is listed and is not excluded, or
                                                        "NO" if the disease is not listed
                                                        or is excluded.
        """
        self._check_patient(patient)

        for dis in patient.diseases:
            if dis.is_present and dis.identifier == self._query:
                return self._diagnosis_present

        return self._diagnosis_excluded

    def __repr__(self):
        return f"DiseasePresencePredicate(query={self._query})"
