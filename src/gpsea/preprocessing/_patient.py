import abc

import typing

from stairval.notepad import Notepad

from gpsea.model import Patient, Cohort

T = typing.TypeVar('T')
"""
The input for `PatientCreator`.

It can be any object that contains the patient data (e.g. a phenopacket).
"""


class PatientCreator(typing.Generic[T], metaclass=abc.ABCMeta):
    """
    `PatientCreator` can create a `Patient` from some input `T`.
    """
    
    @abc.abstractmethod
    def process(
        self,
        item: T,
        notepad: Notepad,
    ) -> typing.Optional[Patient]:
        pass


class CohortCreator(typing.Generic[T]):
    """
    `CohortCreator` creates a cohort from an iterable of some `T` where `T` represents a cohort member.
    """

    def __init__(
        self,
        patient_creator: PatientCreator[T],
    ):
        # Check that we're getting a `PatientCreator`.
        # Unfortunately, we cannot check that `T`s of `PatientCreator` and `CohortCreator` actually match
        # due to Python's loosey-goosey nature.
        assert isinstance(patient_creator, PatientCreator)
        self._pc = patient_creator

    def process(
        self,
        inputs: typing.Iterable[T],
        notepad: Notepad,
    ) -> Cohort:
        patients = []
        patient_labels = set()
        duplicate_pat_labels = set()

        for i, pp in enumerate(inputs):
            sub = notepad.add_subsection(f'patient #{i}')
            patient = self._pc.process(pp, sub)
            if patient is not None:
                if patient.labels in patient_labels:
                    duplicate_pat_labels.add(patient.labels)
                patient_labels.add(patient.labels)
                patients.append(patient)

        # What happens if a sample has
        if len(duplicate_pat_labels) > 0:
            label_summaries = [d.label_summary() for d in duplicate_pat_labels]
            label_summaries.sort()
            notepad.add_error(
                f"Patient ID/s {', '.join(label_summaries)} have a duplicate",
                "Please verify every patient has an unique ID.",
            )

        # We should have >1 patients in the cohort, right?
        if len(patients) <= 1:
            notepad.add_error(
                f'Cohort must include {len(patients)}>1 members',
                'Fix issues in patients to enable the analysis',
            )

        return Cohort.from_patients(patients)
