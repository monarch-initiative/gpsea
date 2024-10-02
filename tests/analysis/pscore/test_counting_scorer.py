import typing

import hpotk
import pytest

from gpsea.model import Patient, Phenotype, SampleLabels, Sex
from gpsea.analysis.pscore import CountingPhenotypeScorer


class TestCountingPhenotypeScorer:

    @pytest.fixture
    def counting_scorer(
        self,
        hpo: hpotk.MinimalOntology,
    ) -> CountingPhenotypeScorer:
        return CountingPhenotypeScorer.from_query_curies(
            hpo=hpo,
            query=(
                "HP:0000478",  # Abnormality of the eye
                "HP:0000407",  # Sensorineural hearing impairment
            ),
        )

    @pytest.mark.parametrize(
        "curies, expected",
        [
            (
                # Exact term matches
                ["HP:0000478", "HP:0000407"],
                2,
            ),
            (
                # One match: Old-aged sensorineural hearing impairment
                [
                    "HP:0040113",
                ],
                1,
            ),
            (
                # Unrelated: Seizure and Arachnodactyly
                ["HP:0001250", "HP:0001166"],
                0,
            ),
        ],
    )
    def test_a_patient(
        self,
        curies: typing.Sequence[str],
        expected: int,
        counting_scorer: CountingPhenotypeScorer,
    ):
        patient = Patient.from_raw_parts(
            labels=SampleLabels("test"),
            sex=Sex.UNKNOWN_SEX,
            phenotypes=(
                Phenotype(
                    hpotk.TermId.from_curie(curie),
                    is_observed=True,
                )
                for curie in curies
            ),
            measurements=(),
            diseases=(),
            variants=(),
        )

        actual = counting_scorer.score(patient)

        assert actual == expected

    def test_get_question(
        self,
        counting_scorer: CountingPhenotypeScorer,
    ):
        actual = counting_scorer.get_question()
        assert actual == "How many of the query HPO terms (or their descendants) does the individual display"

    def test_creating_scorer_with_term_and_ancestor_fails(
        self,
        hpo: hpotk.MinimalOntology,
    ):
        with pytest.raises(ValueError) as e:
            query = (
                "HP:0001250",  # Seizure
                "HP:0012638",  # Abnormal nervous system physiology
            )
            CountingPhenotypeScorer.from_query_curies(hpo, query)

        assert e.value.args[0] == "Both HP:0001250 and its ancestor term HP:0012638 were found in the query, but query terms must not include a term and its ancestor"

    def test_creating_scorer_with_unknown_term_fails(
        self,
        hpo: hpotk.MinimalOntology,
    ):
        with pytest.raises(ValueError) as e:
            query = (
                "HP:0",  # Does not exist
                "HP:1",  # Does not exist
            )
            CountingPhenotypeScorer.from_query_curies(hpo, query)

        assert e.value.args[0] == "The query HP:0 was not found in the HPO"
