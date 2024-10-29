import abc
import typing

import hpotk

from gpsea.model import Patient, Age, Timeline
from .._api import Endpoint
from .._base import Survival


class EndpointBase(Endpoint, metaclass=abc.ABCMeta):

    def __init__(
        self,
        timeline: Timeline,
    ):
        self._timeline = timeline

    def _compute_survival(
        self,
        age: typing.Optional[Age],
        is_censored: bool,
    ) -> typing.Optional[Survival]:
        if age is None or age.timeline != self._timeline:
            return None
        else:
            return Survival(
                value=age.days,
                is_censored=is_censored,
            )


class Death(EndpointBase):

    @property
    def name(self) -> str:
        return "Age of death"

    @property
    def description(self) -> str:
        return f"Compute time until {self._timeline.name.lower()} death"

    @property
    def variable_name(self) -> str:
        return "Age of death"

    def compute_survival(
        self,
        patient: Patient,
    ) -> typing.Optional[Survival]:
        # If the patient is alive we use the current age as `value` and `censored=True`
        # If the patient is deceased, we use the age at death as a `value` and `censored=False`
        if patient.vital_status is None:
            # Absence of the vital status prevents reasoning about death.
            return None

        if patient.vital_status.is_deceased:
            return self._compute_survival(
                age=patient.vital_status.age_of_death,
                is_censored=False,
            )
        else:
            # In absence of an explicit information regarding death,
            # we assume the individual was alive at the reported age.
            return self._compute_survival(
                age=patient.age,
                is_censored=True,
            )
        
    def question_base(self) -> str:
        return f"time until {self._timeline.name.lower()} death"

    def __eq__(self, value: object) -> bool:
        return isinstance(value, Death) and self._timeline == value._timeline

    def __hash__(self) -> int:
        return hash((self._timeline,))

    def __str__(self) -> str:
        return f"Death(timeline={self._timeline})"

    def __repr__(self) -> str:
        return str(self)


class PhenotypicFeatureOnset(EndpointBase):

    def __init__(
        self,
        timeline: Timeline,
        hpo: hpotk.MinimalOntology,
        term_id: hpotk.TermId,
    ):
        super().__init__(timeline)

        assert isinstance(hpo, hpotk.MinimalOntology)
        self._hpo = hpo

        assert isinstance(term_id, hpotk.TermId)
        self._term_id = term_id

        assert term_id in hpo, f"`term_id` {term_id.value} is not in HPO {hpo.version}"

    @property
    def name(self) -> str:
        return f"Onset of {self._hpo.get_term_name(self._term_id)}"

    @property
    def description(self) -> str:
        return f"Compute time until onset of {self._hpo.get_term_name(self._term_id)}"

    @property
    def variable_name(self) -> str:
        return f"Onset of {self._term_id.value}"

    def compute_survival(
        self,
        patient: Patient,
    ) -> typing.Optional[Survival]:
        # Search the present phenotypes. If the individual is annotated with multiple present descendants
        # (e.g. Focal seizure, Clonic seizure) of the target term (e.g. Seizure) then choose
        # the earliest onset, because the onset of the ancestor should be observable
        # since the onset of the first descendant.

        earliest_onset = None
        for present in patient.present_phenotypes():
            # Check if the onset is available ...
            if present.onset is not None and present.onset.timeline == self._timeline:
                # ... and if the individual is annotated with the target HPO or its descendant.
                if present.identifier == self._term_id or any(
                    anc == self._term_id
                    for anc in self._hpo.graph.get_ancestors(present)
                ):
                    if earliest_onset is None:
                        earliest_onset = present.onset
                    else:
                        earliest_onset = min(earliest_onset, present.onset)

        if earliest_onset is None:
            # Phenotype was not found, use the age of the individual and right-censor
            return self._compute_survival(
                age=patient.age,
                is_censored=True,
            )
        else:
            # Phenotype was found, use the earliest onset.
            return self._compute_survival(
                age=earliest_onset,
                is_censored=False,
            )

    def __eq__(self, value: object) -> bool:
        return (
            isinstance(value, PhenotypicFeatureOnset)
            and self._timeline == value._timeline
            and self._term_id == value._term_id
            and self._hpo.version == value._hpo.version
        )

    def __hash__(self) -> int:
        return hash((self._timeline, self._term_id, self._hpo.version))

    def __str__(self) -> str:
        return (
            "PhenotypicFeatureOnset("
            f"timeline={self._timeline}, "
            f"onset={self._term_id}, "
            f"hpo={self._hpo.version})"
        )

    def __repr__(self) -> str:
        return str(self)


class DiseaseOnset(EndpointBase):

    def __init__(
        self,
        timeline: Timeline,
        disease_id: hpotk.TermId,
    ):
        super().__init__(timeline)

        assert isinstance(disease_id, hpotk.TermId)
        self._disease_id = disease_id

    @property
    def name(self) -> str:
        return f"Onset of {self._disease_id.value}"

    @property
    def description(self) -> str:
        return f"Compute time until {self._disease_id.value} onset"

    @property
    def variable_name(self) -> str:
        return f"Onset of {self._disease_id.value}"

    def compute_survival(
        self,
        patient: Patient,
    ) -> typing.Optional[Survival]:
        for disease in patient.present_diseases():
            if disease.identifier == self._disease_id:
                return self._compute_survival(
                    age=disease.onset,
                    is_censored=False,
                )

        return self._compute_survival(
            age=patient.age,
            is_censored=True,
        )
        
    def question_base(self) -> str:
        return f"time until {self._timeline.name.lower()} diagnosis of {self._disease_id.value}"

    def __eq__(self, value: object) -> bool:
        return (
            isinstance(value, DiseaseOnset)
            and self._timeline == value._timeline
            and self._disease_id == value._disease_id
        )

    def __hash__(self) -> int:
        return hash((self._timeline, self._disease_id))

    def __str__(self) -> str:
        return "DiseaseOnset(" \
            f"timeline={self._timeline}, " \
            f"disease_id={self._disease_id})"

    def __repr__(self) -> str:
        return str(self)


GESTATIONAL_DEATH = Death(timeline=Timeline.GESTATIONAL)
POSTNATAL_DEATH = Death(timeline=Timeline.POSTNATAL)


def death(
    timeline: typing.Literal["gestational", "postnatal"] = "postnatal",
) -> Endpoint:
    """
    Get :class:`~gpsea.analysis.temporal.Endpoint` for computing time
    until death of an individual or until the individual is lost from the study
    without knowing about the time of death.

    The time of death is computed from individual's vital status with the following rules:
     
    * If the individual is labeled as :attr:`~gpsea.model.Status.DECEASED`,
      we compute the survival from the age of death.
    * If the individual is :attr:`~gpsea.model.Status.ALIVE` or the status is missing,
      we use the age at last encounter as the censored survival.
    * If the age at last encounter is missing or if the age does not match the target timeline
      (e.g. `timeline==postnatal` but the individual has `gestational` age) then we cannot compute the survival
      and the endpoint returns `None`.
    """
    age_timeline = _decode_timeline(timeline)
    if age_timeline == Timeline.GESTATIONAL:
        return GESTATIONAL_DEATH
    elif age_timeline == Timeline.POSTNATAL:
        return POSTNATAL_DEATH
    else:
        raise ValueError(f"Unsupported timeline {timeline}")


def disease_onset(
    disease_id: typing.Union[str, hpotk.TermId],
    timeline: typing.Literal["gestational", "postnatal"] = "postnatal",
) -> Endpoint:
    """
    Get :class:`~gpsea.analysis.temporal.Endpoint` to compute time
    until onset of a disease or until the individual is lost from the study.

    The onset of diagnosis is computed from the onset field of
    a :class:`~gpsea.model.Disease` with the following rules:
     
    * If the individual is diagnosed with the target disease and its onset is known,
      then the survival is computed from the disease onset.
    * If the individual is *not* diagnosed with the disease and the age at last encounter is known,
      this age is used as censored survival.
    * If the age at last encounter is missing or if the age does not match the target timeline
      (e.g. `timeline==postnatal` but the individual's age is on gestational timeline) then we cannot compute
      the time until disease onset and the endpoint returns `None`.
    """
    age_timeline = _decode_timeline(timeline)
    disease_id = _validate_term_id(disease_id)

    return DiseaseOnset(
        timeline=age_timeline,
        disease_id=disease_id,
    )


def hpo_onset(
    hpo: hpotk.MinimalOntology,
    term_id: typing.Union[str, hpotk.TermId],
    timeline: typing.Literal["gestational", "postnatal"] = "postnatal",
) -> Endpoint:
    """
    Get :class:`~gpsea.analysis.temporal.Endpoint` to compute time
    until onset of an HPO term or until the individual is lost from the study.

    The HPO term onset is computed from the onset field of
    a :class:`~gpsea.model.Phenotype` with the following rules:
     
    * If the individual is annotated with the target HPO term and its onset is known,
      then the survival is computed from the term's onset.
    * If the individual is *not* diagnosed with the term and the age at last encounter is known,
      this age is used as censored survival.
    * If the age at last encounter is missing or if the age does not match the target timeline
      (e.g. `timeline==postnatal` but the individual has `gestational` age) then we cannot compute
      time until phenotype onset and the endpoint returns `None`.
    """
    age_timeline = _decode_timeline(timeline)
    term_id = _validate_term_id(term_id)

    return PhenotypicFeatureOnset(
        timeline=age_timeline,
        hpo=hpo,
        term_id=term_id,
    )


def _decode_timeline(
    timeline: typing.Literal["gestational", "postnatal"] = "postnatal",
) -> Timeline:
    if timeline == "gestational":
        return Timeline.GESTATIONAL
    elif timeline == "postnatal":
        return Timeline.POSTNATAL
    else:
        raise ValueError(f"Unsupported timeline {timeline}")


def _validate_term_id(
    term_id: typing.Union[str, hpotk.TermId],
) -> hpotk.TermId:
    if isinstance(term_id, str):
        return hpotk.TermId.from_curie(term_id)
    elif isinstance(term_id, hpotk.TermId):
        return term_id
    else:
        raise ValueError(
            f"`term_id` must be a `str` or `hpotk.TermId` but was {term_id}"
        )
