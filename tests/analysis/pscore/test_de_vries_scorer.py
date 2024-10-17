import typing

import hpotk
import pytest

from gpsea.analysis.pscore import DeVriesPhenotypeScorer
from gpsea.model import Patient, SampleLabels, Phenotype, Sex

intrauterine_growth_retardation = 'HP:0001511'
small_for_gestational_age = 'HP:0001518'
arachnodactyly = "HP:0001166"
seizure = "HP:0001250"
sensorineural_hearing_impairment = 'HP:0000407'
intellectual_disability_mild = 'HP:0001256'
intellectual_disability_profound = 'HP:0002187'
microcephaly = 'HP:0000252'
short_stature = 'HP:0004322'
hypertelorism = 'HP:0000316'
posteriorly_rotated_ears = 'HP:0000358'
underdeveloped_crus_of_the_helix = 'HP:0009898'  # external ear morphology
ventricular_septal_defect = 'HP:0001629'
metacarpal_synostosis = 'HP:0009701'  # hand morphology
hypospadias = 'HP:0000047'


class TestDeVriesScorer:

    @pytest.fixture
    def devries_scorer(
        self,
        hpo: hpotk.MinimalOntology,
    ) -> DeVriesPhenotypeScorer:
        return DeVriesPhenotypeScorer(hpo=hpo)

    @pytest.mark.parametrize(
        "term_set, expected",
        [
            ([intrauterine_growth_retardation], 2,),
            ([intrauterine_growth_retardation, small_for_gestational_age], 2,),  # superfluous, still should be 2
            ([sensorineural_hearing_impairment, ], 0,),  # Unrelated
            ([seizure, arachnodactyly], 0,),  # Unrelated
            ([intrauterine_growth_retardation, intellectual_disability_mild], 3,),
            ([intrauterine_growth_retardation, intellectual_disability_profound], 4,),
            ([intrauterine_growth_retardation, intellectual_disability_profound, microcephaly], 5,),
            ([intrauterine_growth_retardation, intellectual_disability_profound, microcephaly, short_stature], 6,),
            ([intrauterine_growth_retardation, intellectual_disability_profound, microcephaly, short_stature,
              hypertelorism], 6,),
            ([intrauterine_growth_retardation, intellectual_disability_profound, microcephaly, short_stature,
              hypertelorism, posteriorly_rotated_ears], 8,),
            ([intrauterine_growth_retardation, intellectual_disability_profound, microcephaly, short_stature,
              hypertelorism, posteriorly_rotated_ears, underdeveloped_crus_of_the_helix], 8,),
            ([intrauterine_growth_retardation, intellectual_disability_profound, microcephaly, short_stature,
              hypertelorism, posteriorly_rotated_ears, underdeveloped_crus_of_the_helix,
              ventricular_septal_defect], 9,),
            ([intrauterine_growth_retardation, intellectual_disability_profound, microcephaly, short_stature,
              hypertelorism, posteriorly_rotated_ears, underdeveloped_crus_of_the_helix,
              ventricular_septal_defect, hypospadias], 10,),
            ([intrauterine_growth_retardation, intellectual_disability_profound, microcephaly, short_stature,
              hypertelorism, posteriorly_rotated_ears, underdeveloped_crus_of_the_helix,
              ventricular_septal_defect, hypospadias, metacarpal_synostosis], 10,),

        ],
    )
    def test_a_patient(
            self,
            term_set: typing.Sequence[str],
            expected: int,
            devries_scorer: DeVriesPhenotypeScorer,
    ):
        patient = Patient.from_raw_parts(
            labels=SampleLabels("test"),
            sex=Sex.UNKNOWN_SEX,
            age_at_death=None,
            phenotypes=(
                Phenotype.from_raw_parts(
                    term_id=curie,
                    is_observed=True,
                )
                for curie in term_set
            ),
            measurements=(),
            diseases=(),
            variants=(),
        )

        actual = devries_scorer.score(patient)

        assert actual == expected
