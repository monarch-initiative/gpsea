import typing

import hpotk

from gpsea.model import Patient
from ._api import PhenotypeScorer


class CountingPhenotypeScorer(PhenotypeScorer):
    """
    `CountingPhenotypeScorer` assigns the patient with a phenotype score
    that is equivalent to the count of present phenotypes that are either
    an exact match to the `query` terms or their descendants.

    For instance, we may want to count whether an individual has brain, liver, kidney, and skin abnormalities.
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

            # Now check that the term IDs are HPO term IDs.
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


class DeVriesPhenotypeScorer(PhenotypeScorer):
    """
    `DeVriesPhenotypeScorer` computes "adapted De Vries Score"
    as described in `Feenstra et al. <https://pubmed.ncbi.nlm.nih.gov/21712853>`_.
    """

    def __init__(
        self,
        hpo: hpotk.MinimalOntology,
    ):
        self._hpo = hpo

        # severe and profound GDD
        self._gdd_tids = {
            'HP:0011344': 2, 'HP:0012736': 2,
            'HP:0011342': 1, 'HP:0011343': 1, 'HP:0001263': 1,
        }

        # mild, moderate, and unspecified GDD (borderline has 0.5)
        self._idd_tids = {
            'HP:0010864': 2, 'HP:0002187': 2,
            'HP:0001256': 1, 'HP:0002342': 1, 'HP:0001249': 1,
            'HP:0006889': 0.5,
        }

    def _developmental_delay_score(
        self,
        observed_term_ids: typing.Iterable[str],
    ) -> float:
        """
        Calculate the dev delay component of the score
        
        Args:
            observed_term_ids: terms observed in patient

        Returns: a score between 0 and 2
        """
        # Check GDD terms with higher priority than ID terms.
        # Global developmental delay
        for t in observed_term_ids:
            if t in self._gdd_tids:
                return self._gdd_tids[t]
        
        # Intellectual disability
        for t in observed_term_ids:
            if t in self._idd_tids:
                return self._idd_tids[t]
        
        return 0

    def _term_or_descendant(
        self,
        target_tid: str,
        observed_term_ids: typing.Iterable[str],
    ) -> int:
        """
        Args:
            target_tid: term of interest
            observed_term_ids: all terms observed in patient

        Returns:
            1 if the term or any descendant is present in the patient, otherwise 0
        """
        for term_id in observed_term_ids:
            if term_id == target_tid \
               or any(ancestor == target_tid for ancestor in self._hpo.graph.get_ancestors(term_id)):
                return 1
        
        return 0

    def _term_or_descendant_count(
        self,
        target_tid: str,
        observed_term_ids: typing.Iterable[str],
    ) -> int:
        """
        Args:
            target_tid: term of interest
            observed_term_ids: all terms observed in patient

        Returns:
            the total count of the terms equal to or descending from the target_tid
        """
        total_count = 0
        for term_id in observed_term_ids:
            for desc_tid in self._hpo.graph.get_ancestors(term_id, include_source=True):
                if desc_tid.value == target_tid:
                    total_count += 1
        return total_count

    def _postnatal_growth_score(
        self,
        observed_term_ids: typing.Iterable[str],
    ) -> int:
        """
        Calculate the postnatal growth component of the score.
        
        Args:
            observed_term_ids: terms observed in patient

        Returns: an `int` (between 0 and 2)
        """
        microcephaly = 'HP:0000252'
        short_stature = 'HP:0004322'
        macrocephaly = 'HP:0000256'
        tall_stature = 'HP:0000098'
        total_count = 0
        for tid in (microcephaly, short_stature, macrocephaly, tall_stature):
            total_count += self._term_or_descendant(tid, observed_term_ids)
        if total_count > 2:
            raise ValueError(f"Inconsistent annotations for postnatal growth score {total_count}:  {observed_term_ids}")
        return total_count

    def _facial_dysmorphism_score(
        self,
        observed_term_ids: typing.Collection[str],
    ) -> int:
        """
        This section assigns two points if two or more anomalies are identified in the following
        categories: hypertelorism, nasal anomalies and ear anomalies. Our implementation counts the total
        number of terms or descendants of the hypertelorism, Abnormal external nose morphology, and
        Abnormal pinna morphology.

        Args:
            observed_term_ids: terms observed in patient

        Returns: facial dysmorphism score (between 0 and 2)

        """
        hypertelorism = 'HP:0000316'
        external_nose = 'HP:0010938'
        pinna_morphology = 'HP:0000377'

        # No need to inspect descendants since Hypertelorism has none.
        total_count = 1 if hypertelorism in observed_term_ids else 0
        total_count += self._term_or_descendant_count(target_tid=external_nose, observed_term_ids=observed_term_ids)
        total_count += self._term_or_descendant_count(target_tid=pinna_morphology, observed_term_ids=observed_term_ids)
        if total_count > 1:
            return 2
        else:
            return 0

    def _congenital_score(
        self,
        observed_term_ids: typing.Iterable[str],
    ) -> int:
        """
        Non-facial dysmorphism and congenital abnormalities component.
        One point is assigned for either the corresponding HPO terms or any of their descendents up to a maximum of 2.
        
        Args:
            observed_term_ids:  terms observed in patient

        Returns:   Non-facial dysmorphism and congenital abnormalities score (between 0 and 2)

        """
        hypospadias = 'HP:0000047'
        abnormal_hand_morphology = 'HP:0005922'
        abnormal_heart_morphology = 'HP:0001627'
        # total_count = len([t for t in observed_term_ids if t == hypospadias])
        total_count = self._term_or_descendant_count(
            target_tid=hypospadias, observed_term_ids=observed_term_ids,
        )
        total_count += self._term_or_descendant_count(target_tid=abnormal_hand_morphology,
                                                      observed_term_ids=observed_term_ids)
        total_count += self._term_or_descendant_count(target_tid=abnormal_heart_morphology,
                                                      observed_term_ids=observed_term_ids)
        return min(2, total_count)

    def _prenatal_growth_score(
        self,
        observed_term_ids: typing.Iterable[str],
    ) -> int:
        """
        Two points are assigned if Prenatal-onset growth retardation is present.

        Args:
            observed_term_ids: list of strings with term identifiers or observed HPO terms

        Returns: score between 0 and 2

        """
        small_for_gestational_age = 'HP:0001518'
        intrauterine_growth_retardation = 'HP:0001511'
        targets = (small_for_gestational_age, intrauterine_growth_retardation)
        for tid in observed_term_ids:
            if tid in targets:
                return 2
        return 0

    def score(self, patient: Patient) -> float:
        """
        Calculate score based on list of strings with term identifiers or observed HPO terms.
        
        Args:
            patient: list of strings with term identifiers or observed HPO terms

        Returns: de Vries score between 0 and 10

        """
        observed_term_ids = tuple(tid.identifier.value for tid in patient.present_phenotypes())

        delay_score = self._developmental_delay_score(observed_term_ids)
        growth_score = self._postnatal_growth_score(observed_term_ids)
        facial_score = self._facial_dysmorphism_score(observed_term_ids)
        congen_score = self._congenital_score(observed_term_ids)
        prenatal_score = self._prenatal_growth_score(observed_term_ids)
        
        return delay_score + growth_score + facial_score + congen_score + prenatal_score
