import abc
import logging

import hpotk
import typing
# pyright: reportGeneralTypeIssues=false

from genophenocorr.model import ProteinMetadata, Phenotype, Patient, Variant
from ._phenotype import PhenotypeCreator
from ._api import FunctionalAnnotator


T = typing.TypeVar('T')


class PatientCreator(typing.Generic[T], metaclass=abc.ABCMeta):
    """A metaclass that can be used to establish a class that creates a Patient object 

    Methods:
        create_patient(item:Generic): Creates a Patient from the data in a given item
    """

    @abc.abstractmethod
    def create_patient(self, item: T) -> Patient:
        """Creates a Patient from the data in a given item

        Args:
            item (Generic[T]): An object with subject data
        Returns:
            Patient: A Patient object
        """
        pass