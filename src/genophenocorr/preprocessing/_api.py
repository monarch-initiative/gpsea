import abc
import typing

from genophenocorr.model import VariantCoordinates, Genotype, ProteinMetadata, TranscriptInfoAware, TranscriptCoordinates, TranscriptAnnotation

T = typing.TypeVar('T')


class VariantCoordinateFinder(typing.Generic[T], metaclass=abc.ABCMeta):

    @abc.abstractmethod
    def find_coordinates(self, item: T) -> typing.Tuple[VariantCoordinates, Genotype]:
        """
        Determine :class:`VariantCoordinates` and :class:`Genotype` from an `item` of some sort.

        Raises:
            ValueError: if there is an error of any kind.
        """
        pass


class FunctionalAnnotator(metaclass=abc.ABCMeta):

    @abc.abstractmethod
    def annotate(self, variant_coordinates: VariantCoordinates) -> typing.Sequence[TranscriptAnnotation]:
        """
        Compute functional annotations for the variant coordinates. The annotations can be empty.

        Returns: a sequence of transcript annotations
        Raises:
            ValueError if the annotation cannot proceed due to the remote resource being offline, etc.
        """
        pass


class TranscriptCoordinateService(metaclass=abc.ABCMeta):
    """
    `TranscriptCoordinateService` gets transcript (tx) coordinates for a given transcript ID.
    """

    @abc.abstractmethod
    def fetch(self, tx: typing.Union[str, TranscriptInfoAware]) -> TranscriptCoordinates:
        """
        Get tx coordinates for a tx ID or an entity that knows about the tx ID.

        The method will raise an exception in case of an issue.

        Args:
            tx: a `str` with tx ID (e.g. `NM_002834.5`) or an entity that knows about the transcript ID
            (e.g. :class:`genophenocorr.model.TranscriptAnnotation`).

        Returns: the transcript coordinates.
        """
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
