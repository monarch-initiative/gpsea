import abc
import typing

import hpotk

from genophenocorr.model import Patient

from .._api import PolyPredicate, PatientCategory, PatientCategories, Categorization

P = typing.TypeVar("P")
"""
Phenotype entity of interest, such as :class:`hpotk.model.TermId`, representing an HPO term or an OMIM/MONDO term.

However, phenotype entity can be anything as long as it is :class:`typing.Hashable` and comparable
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

    Each phenotype category `P` can be a :class:`hpotk.model.TermId` representing an HPO term or an OMIM/MONDO term.

    Most of the time, only one category is investigated, and :attr:`phenotype` will return a sequence with
    one element only (e.g. *Arachnodactyly* `HP:0001166`). However, multiple categories can be sought as well,
    such as when the predicate bins the patient into one of discrete diseases
    (e.g. Geleophysic dysplasia, Marfan syndrome, ...). Then the predicate should return a sequence of disease
    identifiers.
    """

    @property
    @abc.abstractmethod
    def phenotypes(self) -> typing.Sequence[P]:
        """
        Get the phenotype entities of interest.
        """
        pass


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
        self._query_label = self._hpo.get_term(query)
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

    def get_question(self) -> str:
        return f"Is {self._query_label} present in the patient?"

    @property
    def phenotypes(self) -> typing.Sequence[hpotk.TermId]:
        # We usually test just a single HPO term, so we return a tuple with a single member.
        return (self._query,)

    def get_categorizations(
        self,
    ) -> typing.Sequence[PhenotypeCategorization[hpotk.TermId]]:
        return self._phenotype_observed, self._phenotype_excluded

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

    The predicate tests if the patient's diseases include a disease ID formatted as a :class:`hpotk.model.TermId`.

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

    def get_question(self) -> str:
        return f"Was {self._query} diagnosed in the patient"

    @property
    def phenotypes(self) -> typing.Sequence[hpotk.TermId]:
        # We usually test just a single Disease, so we return a tuple with a single member.
        return (self._query,)

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


class CountingPhenotypeScorer:
    """
    `CountingPhenotypeScorer` assigns the patient with a phenotype score
    that is equivalent to the count of present phenotypes that are either
    an exact match to the `query` terms or their descendants.
    For instance, we may want to count whether an individual has brain, liver, kidney, and skin abormalities.
    In the case, the query would include the corresponding terms (e.g., Abnormal brain morphology HP:0012443).
    An individual can then have between 0 and 4 phenotype group abnormalities.
    This predicate is intended to be used with the Mann Whitney U test.

    Example
    ^^^^^^^
    We first need to load HPO using HPO toolkit:
    
        >>> import hpotk
        >>> store = hpotk.configure_ontology_store()
        >>> hpo = store.load_minimal_hpo(release='v2023-10-09')
        
    Now we can create `CountingPhenotypeScorer` to test for presence of brain, liver, kidney and skin abnormalities:

        >>> from genophenocorr.analysis.predicate.phenotype import CountingPhenotypeScorer
        >>> phenotype_groups = (
        ...     "HP:0012443",  # Abnormal brain morphology
        ...     "HP:0410042",  # Abnormal liver morphology
        ...     "HP:0012210",  # Abnormal renal morphology
        ...     "HP:0011121",  # Abnormal skin morphology
        ... )
        >>> scorer = CountingPhenotypeScorer.from_query_curies(
        ...     hpo=hpo,
        ...     query=phenotype_groups,
        ... )

    """

    @staticmethod
    def from_query_curies(
        hpo: hpotk.MinimalOntology,
        query: typing.Iterable[typing.Union[str, hpotk.TermId]],
    ):
        """
        Create a scorer to test for the number of phenotype terms that fall into the phenotype groups.

        :param hpo: HPO as represented by :class:`~hpotk.ontology.MinimalOntology` of HPO toolkit.
        :param query: an iterable of the top-level terms, either represented as CURIEs (`str`)
          or as term IDs.
        """
        query_term_ids = set()
        for q in query:
            # First check if the query items are Term IDs or curies.
            if isinstance(q, str):
                q = hpotk.TermId.from_curie(q)
            elif isinstance(q, hpotk.TermId):
                pass
            else:
                raise ValueError(
                    f"query argument must be iterable of hpotk TermId's or strings but we found {type(q)}"
                )

            # Now chack that the term IDs are HPO term IDs.
            if q not in hpo:
                raise ValueError(f"The query {q} was not found in the HPO")
            query_term_ids.add(q)

        if len(query_term_ids) == 0:
            raise ValueError("`query` must not be empty")

        # the query terms must not include a term and its ancestor
        for q in query_term_ids:
            for anc in hpo.graph.get_ancestors(q):
                if anc in query_term_ids:
                    raise ValueError(
                        f"Both {q} and its ancestor term {anc} were found in the query, "
                        + "but query terms must not include a term and its ancestor"
                    )

        return CountingPhenotypeScorer(
            hpo=hpo,
            query=query_term_ids,
        )

    def __init__(
        self,
        hpo: hpotk.MinimalOntology,
        query: typing.Iterable[hpotk.TermId],
    ):
        self._hpo = hpo
        self._query = set(query)

    def get_question(self) -> str:
        return "How many of the query HPO terms (or their descendants) does the individual display"

    def score(
        self,
        patient: Patient,
    ) -> float:
        """
        Get the count (number) of terms in the query set
        that have matching terms (exact matches or descendants) in the patient.
        Do not double count if the patient has two terms
        (e.g., two different descendants) of one of the query terms.
        """
        count = 0
        for q in self._query:
            for pf in patient.present_phenotypes():
                hpo_id = pf.identifier
                if hpo_id == q or any(
                    anc == q for anc in self._hpo.graph.get_ancestors(hpo_id)
                ):
                    count += 1
                    # We break the inner loop to continue the outer.
                    break

        # A sanity check - we cannot produce more counts than there are categories!
        assert 0 <= count <= len(self._query)

        return count

    def __call__(
        self,
        *args: typing.Any,
        **kwds: typing.Any,
    ) -> float:
        # TODO: move to `PhenotypeScorer` API.
        assert len(args) == 1 and isinstance(args[0], Patient), 'The first argument must be an instance of `Patient`'
        assert len(kwds) == 0, 'We do not take any key-word arguments'
        return self.score(args[0])
