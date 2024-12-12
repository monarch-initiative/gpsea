import typing

import hpotk

from gpsea.model import Patient

from ._api import PhenotypeClassifier, PhenotypeCategorization, YES, NO


class HpoClassifier(PhenotypeClassifier[hpotk.TermId]):
    """
    `HpoClassifier` tests if a patient is annotated with an HPO term.

    Note, `query` must be a term of the provided `hpo`!

    See :ref:`hpo-classifier` section for an example usage.

    :param hpo: HPO ontology
    :param query: the HPO term to test
    :param missing_implies_phenotype_excluded: `True` if lack of an explicit annotation implies term's absence`.
    """

    def __init__(
        self,
        hpo: hpotk.MinimalOntology,
        query: hpotk.TermId,
        missing_implies_phenotype_excluded: bool = False,
    ):
        assert isinstance(hpo, hpotk.MinimalOntology)
        self._hpo = hpo
        assert isinstance(query, hpotk.TermId)
        self._query = query
        self._query_label = self._hpo.get_term_name(query)
        assert self._query_label is not None, f"Query {query} is in HPO"
        assert isinstance(missing_implies_phenotype_excluded, bool)
        self._missing_implies_phenotype_excluded = missing_implies_phenotype_excluded
        self._phenotype_observed = PhenotypeCategorization(
            category=YES,
            phenotype=self._query,
        )
        self._phenotype_excluded = PhenotypeCategorization(
            category=NO,
            phenotype=self._query,
        )
        # Some tests depend on the order of `self._categorizations`.
        self._categorizations = (self._phenotype_observed, self._phenotype_excluded)

    @property
    def name(self) -> str:
        return "HPO Classifier"

    @property
    def description(self) -> str:
        return f"Test for presence of {self._query_label} [{self._query.value}]"

    @property
    def variable_name(self) -> str:
        return self._query.value

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
        self, patient: Patient,
    ) -> typing.Optional[PhenotypeCategorization[hpotk.TermId]]:
        self._check_patient(patient)

        if len(patient.phenotypes) == 0:
            return None

        for phenotype in patient.phenotypes:
            if phenotype.is_present:
                if self._query == phenotype.identifier or any(
                    self._query == anc
                    for anc in self._hpo.graph.get_ancestors(phenotype)
                ):
                    return self._phenotype_observed
            else:
                if self._missing_implies_phenotype_excluded:
                    return self._phenotype_excluded
                else:
                    if phenotype.identifier == self._query or any(
                        phenotype.identifier == anc
                        for anc in self._hpo.graph.get_ancestors(self._query)
                    ):
                        return self._phenotype_excluded

        return None

    def __eq__(self, value: object) -> bool:
        return isinstance(value, HpoClassifier) \
            and self._hpo.version == value._hpo.version \
            and self._query == value._query \
            and self._missing_implies_phenotype_excluded == value._missing_implies_phenotype_excluded
    
    def __hash__(self) -> int:
        return hash((self._hpo.version, self._query, self._missing_implies_phenotype_excluded))

    def __repr__(self):
        return f"HpoClassifier(query={self._query})"


class DiseasePresenceClassifier(PhenotypeClassifier[hpotk.TermId]):
    """
    `DiseasePresenceClassifier` tests if an individual was diagnosed with a disease.

    :param disease_id_query: a disease identifier formatted either as a CURIE `str` (e.g. ``OMIM:256000``)
      or as a :class:`~hpotk.TermId`.
    """

    def __init__(
        self,
        disease_id_query: typing.Union[str, hpotk.TermId],
    ):
        if isinstance(disease_id_query, str):
            self._query = hpotk.TermId.from_curie(disease_id_query)
        elif isinstance(disease_id_query, hpotk.TermId):
            self._query = disease_id_query
        else:
            raise AssertionError

        self._diagnosis_present = PhenotypeCategorization(
            category=YES,
            phenotype=disease_id_query,
        )
        self._diagnosis_excluded = PhenotypeCategorization(
            category=NO,
            phenotype=disease_id_query,
        )

    @property
    def name(self) -> str:
        return "Disease Classifier"

    @property
    def description(self) -> str:
        return f"Partition based on a diagnosis of {self._query.value}"

    @property
    def variable_name(self) -> str:
        return self._query.value

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
        self._check_patient(patient)

        for dis in patient.diseases:
            if dis.is_present and dis.identifier == self._query:
                return self._diagnosis_present

        return self._diagnosis_excluded

    def __eq__(self, value: object) -> bool:
        return isinstance(value, DiseasePresenceClassifier) \
            and self._query == value._query
    
    def __hash__(self) -> int:
        return hash((self._query,))

    def __repr__(self):
        return f"DiseasePresenceClassifier(query={self._query})"
