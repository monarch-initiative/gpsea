import abc
import typing

from genophenocorr.model import Variant, VariantCoordinates, ProteinMetadata, TranscriptInfoAware, TranscriptCoordinates

T = typing.TypeVar('T')


class VariantCoordinateFinder(typing.Generic[T], metaclass=abc.ABCMeta):

    @abc.abstractmethod
    def find_coordinates(self, item: T) -> VariantCoordinates:
        """
        Determine :class:`VariantCoordinates` from an `item` of some sort.
        """
        pass


class FunctionalAnnotator(metaclass=abc.ABCMeta):

    @abc.abstractmethod
    def annotate(self, variant_coordinates: VariantCoordinates) -> Variant:
        pass


class TranscriptCoordinateService(metaclass=abc.ABCMeta):

    @abc.abstractmethod
    def fetch(self, tx: TranscriptInfoAware) -> TranscriptCoordinates:
        pass


class ProteinMetadataService(metaclass=abc.ABCMeta):
    """A metaclass that can be used to establish a class that creates ProteinMetadata objects

    Methods:
        annotate(protein_id:str): Gets metadata and creates ProteinMetadata for given protein ID
    """

    @abc.abstractmethod
    def annotate(self, protein_id: str) -> typing.Sequence[ProteinMetadata]:
        """Get metadata for given protein ID

        Args:
            protein_id (string): A protein ID
        Returns:
            Sequence[ProteinMetadata]: A sequence of ProteinMetadata objects, or an empty sequence if no data was found
        """
        pass
