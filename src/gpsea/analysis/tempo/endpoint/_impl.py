import abc
import typing

import hpotk

from gpsea.model import Patient, Age, AgeKind
from .._api import Endpoint
from .._base import Survival


class EndpointBase(Endpoint, metaclass=abc.ABCMeta):

    def __init__(
        self,
        kind: AgeKind,
    ):
        self._kind = kind

    def _compute_survival(
        self,
        age: typing.Optional[Age],
        is_censored: bool,
    ) -> typing.Optional[Survival]:
        if age is None or age.kind != self._kind:
            return None
        else:
            return Survival(
                value=age.days,
                is_censored=is_censored,
            )


class Death(EndpointBase):

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

    def __eq__(self, value: object) -> bool:
        return isinstance(value, Death) and self._kind == value._kind

    def __hash__(self) -> int:
        return hash((self._kind,))

    def __str__(self) -> str:
        return f"Death(kind={self._kind})"

    def __repr__(self) -> str:
        return str(self)


class PhenotypicFeatureOnset(EndpointBase):

    def __init__(
        self,
        kind: AgeKind,
        hpo: hpotk.MinimalOntology,
        term_id: hpotk.TermId,
    ):
        super().__init__(kind)

        assert isinstance(hpo, hpotk.MinimalOntology)
        self._hpo = hpo

        assert isinstance(term_id, hpotk.TermId)
        self._term_id = term_id

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
            if present.onset is not None and present.onset.kind == self._kind:
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
            and self._kind == value._kind
            and self._term_id == value._term_id
            and self._hpo.version == value._hpo.version
        )

    def __hash__(self) -> int:
        return hash((self._kind, self._term_id, self._hpo.version))

    def __str__(self) -> str:
        return (
            "PhenotypicFeatureOnset("
            f"kind={self._kind}, "
            f"onset={self._term_id}, "
            f"hpo={self._hpo.version})"
        )

    def __repr__(self) -> str:
        return str(self)


class DiseaseOnset(EndpointBase):

    def __init__(
        self,
        kind: AgeKind,
        disease_id: hpotk.TermId,
    ):
        super().__init__(kind)

        assert isinstance(disease_id, hpotk.TermId)
        self._disease_id = disease_id

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

    def __eq__(self, value: object) -> bool:
        return (
            isinstance(value, DiseaseOnset)
            and self._kind == value._kind
            and self._disease_id == value._disease_id
        )

    def __hash__(self) -> int:
        return hash((self._kind, self._disease_id))

    def __str__(self) -> str:
        return "DiseaseOnset(" \
            f"kind={self._kind}, " \
            f"disease_id={self._disease_id})"

    def __repr__(self) -> str:
        return str(self)


GESTATIONAL_DEATH = Death(kind=AgeKind.GESTATIONAL)
POSTNATAL_DEATH = Death(kind=AgeKind.POSTNATAL)


def death(
    kind: typing.Literal["gestational", "postnatal"] = "postnatal",
) -> Endpoint:
    age_kind = _decode_kind(kind)
    if age_kind == AgeKind.GESTATIONAL:
        return GESTATIONAL_DEATH
    elif age_kind == AgeKind.POSTNATAL:
        return POSTNATAL_DEATH
    else:
        raise ValueError(f"Unsupported kind {kind}")


def disease_onset(
    disease_id: typing.Union[str, hpotk.TermId],
    kind: typing.Literal["gestational", "postnatal"] = "postnatal",
) -> Endpoint:
    age_kind = _decode_kind(kind)
    disease_id = _validate_term_id(disease_id)

    return DiseaseOnset(
        kind=age_kind,
        disease_id=disease_id,
    )


def hpo_onset(
    hpo: hpotk.MinimalOntology,
    term_id: typing.Union[str, hpotk.TermId],
    kind: typing.Literal["gestational", "postnatal"] = "postnatal",
) -> Endpoint:
    age_kind = _decode_kind(kind)
    term_id = _validate_term_id(term_id)

    return PhenotypicFeatureOnset(
        kind=age_kind,
        hpo=hpo,
        term_id=term_id,
    )


def _decode_kind(
    kind: typing.Literal["gestational", "postnatal"] = "postnatal",
) -> AgeKind:
    if kind == "gestational":
        return AgeKind.GESTATIONAL
    elif kind == "postnatal":
        return AgeKind.POSTNATAL
    else:
        raise ValueError(f"Unsupported kind {kind}")


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
