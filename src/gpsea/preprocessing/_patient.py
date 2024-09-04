import abc

import typing

import hpotk

from gpsea.model import Patient, Cohort

from ._audit import Auditor, Notepad

T = typing.TypeVar('T')
"""
The input for `PatientCreator`.

It can be any object that contains the patient data (e.g. a phenopacket).
"""


class PatientCreator(typing.Generic[T], Auditor[T, Patient], metaclass=abc.ABCMeta):
    """
    `PatientCreator` can create a `Patient` from some input `T`.

    `PatientCreator` is an `Auditor`, hence the input is sanitized and any errors are reported to the caller.
    """
    pass


class CohortCreator(typing.Generic[T], Auditor[typing.Iterable[T], Cohort]):
    """
    `CohortCreator` creates a cohort from an iterable of some `T` where `T` represents a cohort member.
    """

    def __init__(self, patient_creator: PatientCreator[T]):
        # Check that we're getting a `PatientCreator`.
        # Unfortunately, we cannot check that `T`s of `PatientCreator` and `CohortCreator` actually match
        # due to Python's loosey-goosey nature.
        self._pc = hpotk.util.validate_instance(patient_creator, PatientCreator, 'patient_creator')

    def process(self, inputs: typing.Iterable[T], notepad: Notepad) -> Cohort:
        patients = []
        patient_labels = set()
        duplicate_pat_labels = set()

        for i, pp in enumerate(inputs):
            sub = notepad.add_subsection(f'patient #{i}')
            patient = self._pc.process(pp, sub)
            if patient.labels in patient_labels:
                duplicate_pat_labels.add(patient.labels)
            patient_labels.add(patient.labels)
            patients.append(patient)

        # What happens if a sample has
        if len(duplicate_pat_labels) > 0:
            label_summaries = [d.label_summary() for d in duplicate_pat_labels]
            label_summaries.sort()
            notepad.add_error(f"Patient ID/s {', '.join(label_summaries)} have a duplicate",
                              "Please verify every patient has an unique ID.")

        # We should have >1 patients in the cohort, right?
        if len(patients) <= 1:
            notepad.add_warning(f'Cohort must include {len(patients)}>1 members',
                                'Fix issues in patients to enable the analysis')

        return Cohort.from_patients(patients)
