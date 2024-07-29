import abc
import typing

from genophenocorr.model import VariantCoordinates, ProteinMetadata, TranscriptInfoAware, TranscriptCoordinates, TranscriptAnnotation

T = typing.TypeVar('T')


class VariantCoordinateFinder(typing.Generic[T], metaclass=abc.ABCMeta):

    @abc.abstractmethod
    def find_coordinates(self, item: T) -> typing.Optional[VariantCoordinates]:
        """
        Try to find :class:`VariantCoordinates` from an `item` of some sort.

        The variant coordinates may not be available all the time, 
        and `None` may be returned.

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
    """
    A service for obtaining annotations for a given protein accession ID.

    The annotations include elements of the :class:`ProteinMetadata` class.
    """

    @abc.abstractmethod
    def annotate(self, protein_id: str) -> ProteinMetadata:
        """
        Prepare `ProteinMetadata` for a protein with given `protein_id` accession ID.

        Args:
            protein_id (string): A accession ID `str` (e.g. `NP_001027558.1`)
        Returns:
            ProteinMetadata: a `ProteinMetadata` container with the protein metadata
        Raises:
            `ValueError` if the metadata cannot be generated for *any* reason.
        """
        pass
