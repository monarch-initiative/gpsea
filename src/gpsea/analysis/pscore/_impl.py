import typing

import hpotk

from gpsea.model import Patient
from ._api import PhenotypeScorer


class CountingPhenotypeScorer(PhenotypeScorer):
    """
    `CountingPhenotypeScorer` assigns the patient with a phenotype score
    that is equivalent to the count of present phenotypes that are either
    an exact match to the `query` terms or their descendants.

    For instance, we may want to count whether an individual has brain, liver, kidney, and skin abormalities.
    In the case, the query would include the corresponding terms (e.g., Abnormal brain morphology HP:0012443).
    An individual can then have between 0 and 4 phenotype group abnormalities.
    This predicate is intended to be used with the Mann Whitney U test.
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

    # def __call__(
    #         self,
    #         *args: typing.Any,
    #         **kwds: typing.Any,
    # ) -> float:
    #     # TODO: move to `PhenotypeScorer` API.
    #     assert len(args) == 1 and isinstance(args[0], Patient), 'The first argument must be an instance of `Patient`'
    #     assert len(kwds) == 0, 'We do not take any key-word arguments'
    #     return self.score(args[0])


class DeVriesPhenotypeScorer(PhenotypeScorer):
    """
    `DeVriesPhenotypeScorer` computes "adapted De Vries Score"
    as described in `Feenstra et al <https://pubmed.ncbi.nlm.nih.gov/21712853>`_.
    """

    def __init__(
            self,
            hpo: hpotk.MinimalOntology,
    ):
        self._hpo = hpo

    def _developmental_delay_score(self, observed_term_ids: typing.List[str]) -> float:
        """
        calculate the dev delay component of the score
        Args:
            observed_term_ids: terms observed in patient

        Returns: a score between 0 and 2
        """
        gdd_tids = {'HP:0011344': 2, 'HP:0011344': 2,
                    'HP:0011342': 1, 'HP:0011343': 1, 'HP:0001263': 1}  # severe and profound GDD
        idd_tids = {'HP:0010864': 2, 'HP:0002187': 2,'HP:0001256': 1, 'HP:0002342': 1, 'HP:0001249': 1,
                          'HP:0006889': 0.5}  # mild, moderate, and unspecified GDD (borderline has 0.5)
        # check GDD terms with higher priority than ID terms
        for t in observed_term_ids:
            if t in gdd_tids:
                return gdd_tids.get(t)
        for t in observed_term_ids:
            if t in idd_tids:
                return idd_tids.get(t)
        return 0


    def _postnatal_growth_score(self, observed_term_ids: typing.List[str]) -> float:
        """
        calculate the postnatal growth component of the score
        Args:
            observed_term_ids:

        Returns:

        """
        microcephaly = 'HP:0000252'
        short_stature = 'HP:0004322'
        macrocephaly = 'HP: 0000256'
        tall_stature = 'HP:0000098'
        # to do -- implement me

    def _calculate_score(self, observed_term_ids: typing.List[str]) -> float:
        """
        calculate score based on list of strings with term identifiers or observed HPO terms.
        Args:
            observed_term_ids: list of strings with term identifiers or observed HPO terms

        Returns: de Vries score between 0 and 10

        """
        delay_score = self._developmental_delay_score(observed_term_ids)
        growth_score = self._postnatal_growth_score(observed_term_ids)
        ## todo -- complete
        return delay_score + growth_score

    def score(self, patient: Patient) -> float:
        """
        Compute the score for the `patient`.
        """
        # collect term identifiers as strings for all observed phenotypes
        observed_term_ids = [tid.identifier.value for tid in patient.present_phenotypes()]
        return self._calculate_score(observed_term_ids)
