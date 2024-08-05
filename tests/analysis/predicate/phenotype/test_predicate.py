import typing

import hpotk
import pytest

from genophenocorr.model import Patient, Phenotype, SampleLabels
from genophenocorr.analysis.predicate import PatientCategory
from genophenocorr.analysis.predicate.phenotype import (
    CountingPhenotypePredicate,
    PhenotypeCategorization,
)


class TestCountingPhenotypePredicate:

    @pytest.fixture
    def counting_predicate(
            self,
            hpo: hpotk.MinimalOntology,
    ) -> CountingPhenotypePredicate:
        return CountingPhenotypePredicate.from_query_curies(
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
            counting_predicate: CountingPhenotypePredicate,
    ):
        patient = Patient(
            labels=SampleLabels("test"),
            phenotypes=(
                Phenotype(
                    hpotk.TermId.from_curie(curie),
                    name="Doesn't matter",
                    is_observed=True,
                )
                for curie in curies
            ),
            diseases=(),
            variants=(),
        )

        actual = counting_predicate.test(patient)

        assert actual.phenotype == expected

    def test_phenotypes(
            self,
            counting_predicate: CountingPhenotypePredicate,
    ):
        assert all(
            phenotype in counting_predicate.phenotypes for phenotype in (0, 1, 2)
        )

    def test_categories_and_categorizations(
            self,
            counting_predicate: CountingPhenotypePredicate,
    ):
        categories = list(counting_predicate.get_categories())
        assert categories == [
            PatientCategory(cat_id=0, name='0 phenotypes'), 
            PatientCategory(cat_id=1, name='1 phenotype'), 
            PatientCategory(cat_id=2, name='2 phenotypes'),
        ]

        assert counting_predicate.n_categorizations() == 3

        categorizations = list(counting_predicate.get_categorizations())
        examples = [
            PhenotypeCategorization(PatientCategory(cat_id=0, name='0 phenotypes'), 0),
            PhenotypeCategorization(PatientCategory(cat_id=1, name='1 phenotype'), 1),
        ]

        assert all(example in categorizations for example in examples)

    def test_get_question(
            self,
            counting_predicate: CountingPhenotypePredicate,
    ):
        actual = counting_predicate.get_question()
        assert actual == "How many of the target HPO terms (or their descendants) does the individual display?"

    def test_creating_predicate_with_term_and_ancestor_fails(
            self,
            hpo: hpotk.MinimalOntology,
    ):
        with pytest.raises(ValueError) as e:
            query = (
                "HP:0001250",  # Seizure
                "HP:0012638",  # Abnormal nervous system physiology
            )
            CountingPhenotypePredicate.from_query_curies(hpo, query)

        assert e.value.args[
                   0] == "Both HP:0001250 and its ancestor term HP:0012638 were found in the query, but query terms must not include a term and its ancestor"

    def test_creating_predicate_with_unknown_term_fails(
            self,
            hpo: hpotk.MinimalOntology,
    ):
        with pytest.raises(ValueError) as e:
            query = (
                "HP:0",  # Does not exist
                "HP:1",  # Does not exist
            )
            CountingPhenotypePredicate.from_query_curies(hpo, query)

        assert e.value.args[0] == "The query HP:0 was not found in the HPO"
