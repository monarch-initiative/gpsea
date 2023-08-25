import os
import typing

import hpotk
from google.protobuf.json_format import Parse
from phenopackets import Phenopacket

from genophenocorr.cohort import Cohort
from genophenocorr.patient import PhenopacketPatientCreator


def load_phenopacket_folder(pp_directory: str,
                            patient_creator: PhenopacketPatientCreator) -> Cohort:
    """
    Creates a Patient object for each phenopacket formatted JSON file in the given directory `pp_directory`.

    :param pp_directory: path to a folder with phenopacket JSON files. An error is raised if the path does not point to
      a directory with at least one phenopacket.
    :param patient_creator: patient creator for turning a phenopacket into a :class:`genophenocorr.Patient`
    :return: a cohort made of the phenopackets
    """
    if not os.path.isdir(pp_directory):
        raise ValueError("Could not find directory of Phenopackets.")
    hpotk.util.validate_instance(patient_creator, PhenopacketPatientCreator, 'patient_creator')

    # load Phenopackets
    pps = _load_phenopacket_dir(pp_directory)
    if len(pps) == 0:
        raise ValueError(f"No JSON Phenopackets were found in {pp_directory}")

    # turn phenopackets into patients using patient creator
    patients = [patient_creator.create_patient(pp) for pp in pps]

    # create cohort from patients
    return Cohort.from_patients(patients)


def _load_phenopacket_dir(pp_dir: str) -> typing.Sequence[Phenopacket]:
    patients = []
    for patient_file in os.listdir(pp_dir):
        if patient_file.endswith('.json'):
            phenopacket_path = os.path.join(pp_dir, patient_file)
            pp = _load_phenopacket(phenopacket_path)
            patients.append(pp)
    return patients


def _load_phenopacket(phenopacket_path: str) -> Phenopacket:
    with open(phenopacket_path) as f:
        return Parse(f.read(), Phenopacket())
