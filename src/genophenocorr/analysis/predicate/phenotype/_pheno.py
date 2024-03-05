import abc
import typing

import hpotk

from genophenocorr.model import Patient

from .._api import PolyPredicate, BooleanPredicate, PatientCategory


class PhenotypePolyPredicate(PolyPredicate, metaclass=abc.ABCMeta):
    pass


class BasePhenotypeBooleanPredicate(PhenotypePolyPredicate, BooleanPredicate, metaclass=abc.ABCMeta):
    pass


class PropagatingPhenotypeBooleanPredicate(BasePhenotypeBooleanPredicate):

    def __init__(self, hpo: hpotk.MinimalOntology,
                 phenotypic_feature: hpotk.TermId,
                 missing_implies_excluded: bool = False) -> None:
        self._hpo = hpotk.util.validate_instance(hpo, hpotk.MinimalOntology, 'hpo')
        self._query = hpotk.util.validate_instance(phenotypic_feature, hpotk.TermId, 'phenotypic_feature')
        self._missing_implies_excluded = hpotk.util.validate_instance(missing_implies_excluded, bool,
                                                                      'missing_implies_excluded')

    def get_question(self) -> str:
        query_label = self._hpo.get_term(self._query).name
        return f'Is \'{query_label}\' present in the patient?'

    def test(self, patient: Patient) -> typing.Optional[PatientCategory]:
        self._check_patient(patient)

        if len(patient.phenotypes) == 0:
            return None

        for phenotype in patient.phenotypes:
            if phenotype.is_observed:
                if any(self._query == anc for anc in self._hpo.graph.get_ancestors(phenotype, include_source=True)):
                    return BooleanPredicate.YES
            else:
                if self._missing_implies_excluded:
                    return BooleanPredicate.NO
                else:
                    if any(self._query == desc for desc in
                           self._hpo.graph.get_descendants(phenotype, include_source=True)):
                        return BooleanPredicate.NO

        return None


class PropagatingPhenotypePredicate(PhenotypePolyPredicate):
    """
    `PropagatingPhenotypePredicate` tests if the `patient` is annotated with a `query` HPO term.

    The predicate returns the following results:

    * :attr:`PropagatingPhenotypePredicate.PRESENT` if the patient is annotated with the `query` term or its descendant
    * :attr:`PropagatingPhenotypePredicate.EXCLUDED` presence of the `query` term or its ancestor was specifically
      excluded in the patient
    * :attr:`PropagatingPhenotypePredicate.NOT_MEASURED` if the patient is not annotated with the `query` and presence of `query`
      was *not* excluded
    """

    PRESENT = PatientCategory(cat_id=0,
                              name='Present',
                              description="""
                              The sample *is* annotated with the tested phenotype feature `q`.

                              This is either because the sample is annotated with `q` (exact match),
                              or because one of sample's annotations is a descendant `q` (annotation propagation).
                              For instance, we tested for a Seizure and the sample *had* a Clonic seizure
                              (a descendant of Seizure).
                              """)  #: :meta hide-value:

    EXCLUDED = PatientCategory(cat_id=1,
                               name='Excluded',
                               description="""
                               We are particular about the sample *not* having the tested feature `q`.

                               In other words, `q` was *excluded* in the sample or the sample is annotated with an excluded ancestor of `q`.

                               For instance, we tested for a Clonic seizure and the sample did *not* have any Seizure, which implies
                               *not* Clonic seizure.
                               """)  #: :meta hide-value:

    NOT_MEASURED = PatientCategory(cat_id=2,
                                   name='Not measured',
                                   description="""
                                   We do not know if the sample has or has not the tested feature.
                                   """)  #: :meta hide-value:

    def __init__(self, hpo: hpotk.MinimalOntology,
                 phenotypic_feature: hpotk.TermId) -> None:
        self._hpo = hpotk.util.validate_instance(hpo, hpotk.MinimalOntology, 'hpo')
        self._query = hpotk.util.validate_instance(phenotypic_feature, hpotk.TermId, 'phenotypic_feature')

    @staticmethod
    def get_categories() -> typing.Sequence[PatientCategory]:
        return PropagatingPhenotypePredicate.PRESENT, PropagatingPhenotypePredicate.EXCLUDED, PropagatingPhenotypePredicate.NOT_MEASURED

    def get_question(self) -> str:
        query_label = self._hpo.get_term(self._query).name
        return f'Is \'{query_label}\' present in the patient?'

    def test(self, patient: Patient) -> typing.Optional[PatientCategory]:
        self._check_patient(patient)

        if len(patient.phenotypes) == 0:
            return None

        for phenotype in patient.phenotypes:
            if phenotype.is_observed:
                if any(self._query == anc for anc in self._hpo.graph.get_ancestors(phenotype, include_source=True)):
                    return PropagatingPhenotypePredicate.PRESENT
            else:
                if any(self._query == desc for desc in self._hpo.graph.get_descendants(phenotype, include_source=True)):
                    return self.EXCLUDED

        return self.NOT_MEASURED

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return f'PropagatingPhenotypePredicate(query={self._query})'


class PhenotypePredicateFactory(metaclass=abc.ABCMeta):
    """
    `PhenotypePredicateFactory` creates a :class:`PhenotypePolyPredicate` to assess a phenotype query.
    """

    @staticmethod
    @abc.abstractmethod
    def get_categories() -> typing.Sequence[PatientCategory]:
        """
        Get a sequence of all categories which the `PhenotypePolyPredicate` provided by this factory can produce.
        """
        pass

    @abc.abstractmethod
    def get_predicate(self, query: hpotk.TermId) -> PhenotypePolyPredicate:
        """
        Get an instance of `PhenotypePolyPredicate` for the phenotype `query`.

        :param query: a term ID of the query phenotype.
        """
        pass


class PropagatingPhenotypeBooleanPredicateFactory(PhenotypePredicateFactory):

    def __init__(self, hpo: hpotk.MinimalOntology,
                 missing_implies_excluded: bool = False):
        self._hpo = hpotk.util.validate_instance(hpo, hpotk.MinimalOntology, 'hpo')
        self._missing_implies_excluded = missing_implies_excluded

    @staticmethod
    def get_categories() -> typing.Sequence[PatientCategory]:
        return PropagatingPhenotypeBooleanPredicate.get_categories()

    def get_predicate(self, query: hpotk.TermId) -> PhenotypePolyPredicate:
        return PropagatingPhenotypeBooleanPredicate(self._hpo, query, self._missing_implies_excluded)


class DiseaseDiagnosisPredicate(BooleanPredicate):
    """
    `DiseaseDiagnosisPredicate` tests if the `patient` has a certain disease

    The predicate returns the following results:

    * :attr:`BooleanPredicate.YES` if the patient has the disease
    * :attr:`BooleanPredicate.NO` if the patient is not annotated to have the disease
    """
    def __init__(self) -> None:
        ## TODO stuff disease diagnosis into INIT
        super().__init__()
        # store disease here

    def test(self, patient: Patient) -> typing.Optional[PatientCategory]:
        self._check_patient(patient)
        ## TODO check if patient has the same disease
        pass

