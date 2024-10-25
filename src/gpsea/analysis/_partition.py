import abc
import os
import typing

from gpsea.model import Patient

from ._util import Summarizable


class Partitioning(Summarizable, metaclass=abc.ABCMeta):
    """
    `Partitioning` is a superclass of all classes that assign a group,
    compute a score or survival for an individual.
    """
    
    @property
    @abc.abstractmethod
    def name(self) -> str:
        """
        Get the name of the partitioning.
        """
        pass

    @property
    @abc.abstractmethod
    def summary(self) -> str:
        """
        Get a description of the partitioning.
        """
        pass

    def summarize(
        self,
        other: typing.Optional[str] = None,
    ):
        """
        Summarize the item while also considering `other` (default `None`).
        """
        if other is None:
            return os.linesep.join((self.name, self.summary))
        else:
            return os.linesep.join((self.name, self.summary, other))

    @staticmethod
    def _check_patient(patient: Patient):
        """
        Check if the `patient` meets the partitioning requirements.
        """
        if not isinstance(patient, Patient):
            raise ValueError(f"patient must be type Patient but was type {type(patient)}")


class ContinuousPartitioning(Partitioning, metaclass=abc.ABCMeta):
    """
    `ContinuousPartitioning` computes a score that is a real number.

    The class is just a marker class at this time.
    """
    pass
