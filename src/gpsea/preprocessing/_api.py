import abc
import io
import os
import sys
import typing

from gpsea.model import (
    VariantCoordinates,
    ProteinMetadata,
    TranscriptInfoAware,
    TranscriptCoordinates,
    TranscriptAnnotation,
    ImpreciseSvInfo,
    FeatureInfo,
    FeatureType,
    Region,
    SimpleProteinFeature
)

from ._audit import NotepadTree

T = typing.TypeVar("T")


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
    def annotate(
        self, variant_coordinates: VariantCoordinates
    ) -> typing.Sequence[TranscriptAnnotation]:
        """
        Compute functional annotations for the variant coordinates. The annotations can be empty.

        Returns: a sequence of transcript annotations
        Raises:
            ValueError if the annotation cannot proceed due to the remote resource being offline, etc.
        """
        pass


class ImpreciseSvFunctionalAnnotator(metaclass=abc.ABCMeta):
    """
    Annotator for large SVs that lack the exact breakpoint coordinates.
    """

    @abc.abstractmethod
    def annotate(self, item: ImpreciseSvInfo) -> typing.Sequence[TranscriptAnnotation]:
        """
        Compute functional annotations for a large SV.

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
    def fetch(
        self,
        tx: typing.Union[str, TranscriptInfoAware],
    ) -> TranscriptCoordinates:
        """
        Get tx coordinates for a tx ID or an entity that knows about the tx ID.

        The method will raise an exception in case of an issue.

        Args:
            tx: a `str` with tx ID (e.g. `NM_002834.5`) or an entity that knows about the transcript ID
              (e.g. :class:`~gpsea.model.TranscriptAnnotation`).

        Returns: the transcript coordinates.
        """
        pass


class GeneCoordinateService(metaclass=abc.ABCMeta):
    """
    `GeneCoordinateService` gets transcript (Tx) coordinates for a gene ID.
    """

    @abc.abstractmethod
    def fetch_for_gene(self, gene: str) -> typing.Sequence[TranscriptCoordinates]:
        """
        Get Tx coordinates for a gene ID.

        The method will raise an exception in case of an issue.

        Args:
            gene: a `str` with tx ID (e.g. `HGNC:3603`)

        Returns:
            typing.Sequence[TranscriptCoordinates]: a sequence of transcript coordinates for the gene.
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


class DataFrameProteinMetadataService(ProteinMetadataService):
    """
    A class for obtaining annotations about a protein from a user-supplied pandas DataFrame.
    We expect to obtain the gene symbol, protein identifier, and regions
    The structure of the DataFrame should include the following columns and column headers
    +---------+----------------+------------------+--------+-------+-------+
    | symbol  | identifier     | region           | start  | end   | length|  
    +---------+----------------+------------------+--------+-------+-------+
    | ITPR1   | NP_001365381.1 | Suppresor domain | 1      | 223   |2758   |
    +---------+----------------+------------------+--------+-------+-------+
    | ITPR1   | NP_001365381.1 | IP3 binding      | 224    | 578   |2758   |
    +---------+----------------+------------------+--------+-------+-------+

    TODO ALSO ALLOW REGION TYPES -- FOR NOW JUST REGION

    """
    import pandas as pd
    def __init__(self, df: pd.DataFrame) -> None:
        super().__init__()
        self._df = df

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
        expected_headers = ["symbol", "identifier", "region", "start", "end", "length"]
        if expected_headers != list(self._df):
            raise ValueError(f"Got unexpected DataFrame headers: {list(self._df)}")
        symbols = self._df['symbol'].unique()
        if len(symbols) != 1:
            raise ValueError(f"Expected a single gene symbol but got {self._df['symbol'].unique()}")
        gene_symbol = symbols[0]
        identifiers = self._df['identifier'].unique()
        if len(identifiers) != 1:
            raise ValueError(f"Expected a single identifier but got {self._df['identifier'].unique()}")
        identifier = identifiers[0]
        lengths = self._df['length'].unique()
        if len(lengths) != 1:
            raise ValueError(f"Expected a single length but got {self._df['length'].unique()}")
        length = lengths[0]
        region_list = list()
        for _, row in self._df.iterrows():
            region_name = row["region"]
            region_start = row["start"]
            region_end = row["end"]
            finfo = FeatureInfo(name=region_name, region=Region(region_start=region_start, region_name=region_name))
            ## TODO allow all RegionTypes, extend Pandas
            pfeature = SimpleProteinFeature(info=finfo, feature_type=FeatureType.REGION)
            region_list.append(pfeature)
        return ProteinMetadata(protein_id=identifier, label=gene_symbol, protein_features=region_list, protein_length=length)


class PreprocessingValidationResult:
    """
    The result of input validation of patient data.
    """

    def __init__(
        self,
        policy: str,
        notepad: NotepadTree,
    ):
        self._policy = policy
        self._notepad = notepad

    @property
    def policy(self) -> str:
        """
        Get the used validation policy.
        """
        return self._policy

    def is_ok(self) -> bool:
        """
        Test if the result is OK considering the used validation policy.

        If `False` is returned, use :func:`summarize` method to print out the results.

        :returns: `True` if the analysis can proceed or `False` if errors/warnings were found.
        """
        if self._policy == 'none':
            # No validation
            return True
        elif self._policy == 'lenient':
            return not self._notepad.has_errors(include_subsections=True)
        elif self._policy == 'strict':
            return not self._notepad.has_errors_or_warnings(include_subsections=True)
        else:
            # Bug, please report to the developers.
            raise ValueError(f'Unexpected policy {self._policy}')

    def summarize(
        self,
        file: io.TextIOBase = sys.stdout,
        indent: int = 2,
    ):
        """
        Summarize the validation results into the provided `file`.

        :param file: where to write the validation summary (e.g. :class:`io.StringIO`, )
        :param indent: a non-negative `int` for the  to indent the output
        """
        assert isinstance(indent, int) and indent >= 0

        file.write(f"Validated under {self._policy} policy")
        file.write(os.linesep)

        n_errors = sum(node.error_count() for node in self._notepad.iterate_nodes())
        n_warnings = sum(node.warning_count() for node in self._notepad.iterate_nodes())
        if n_errors > 0 or n_warnings > 0:
            file.write("Showing errors and warnings")
            file.write(os.linesep)

            for node in self._notepad.iterate_nodes():
                if node.has_errors_or_warnings(include_subsections=True):
                    # We must report the node label even if there are no issues with the node.
                    l_pad = " " * (node.level * indent)
                    file.write(l_pad + node.label)
                    file.write(os.linesep)

                    if node.has_errors():
                        file.write(l_pad + " errors:")
                        file.write(os.linesep)
                        for error in node.errors():
                            file.write(
                                l_pad + " " + error.message + f". {error.solution}"
                                if error.solution
                                else ""
                            )
                            file.write(os.linesep)
                    if node.has_warnings():
                        file.write(l_pad + " warnings:")
                        file.write(os.linesep)
                        for warning in node.warnings():
                            file.write(
                                l_pad + " ·" + warning.message + f". {warning.solution}"
                                if warning.solution
                                else ""
                            )
                            file.write(os.linesep)
        else:
            file.write("No errors or warnings were found")
            file.write(os.linesep)
