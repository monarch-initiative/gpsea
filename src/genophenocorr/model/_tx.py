import typing

import hpotk

from .genome import GenomicRegion


class TranscriptCoordinates:
    """
    `TranscriptCoordinates` knows about genomic region of the transcript, exonic/intronic regions, as well as
    the coding and non-coding regions.
    """

    def __init__(self, identifier: str,
                 region: GenomicRegion,
                 exons: typing.Iterable[GenomicRegion],
                 cds_start: typing.Optional[int],
                 cds_end: typing.Optional[int]):
        self._id = hpotk.util.validate_instance(identifier, str, 'identifier')
        self._region = hpotk.util.validate_instance(region, GenomicRegion, 'region')
        self._exons = tuple(exons)

        if not isinstance(cds_start, int) or not isinstance(cds_end, int):
            raise ValueError(f'CDS coordinates must be ints but were {cds_start}, {cds_end}')
        if not self._region.contains_pos(cds_start + 1) or not self._region.contains_pos(cds_end):
            # `+1` to convert the 0-based start into a 1-based coordinate for a moment
            raise ValueError(f'CDS coordinates {cds_start:,}, {cds_end:,} must be in the tx region: ({self._region.end}, {self._region.end}]')
        self._cds_start = cds_start
        self._cds_end = cds_end

    @property
    def identifier(self) -> str:
        """
        Get transcript's identifier (e.g. `ENST00000123456.7`, `NM_123456.7`).
        """
        return self._id

    @property
    def region(self) -> GenomicRegion:
        """
        Get the genomic region spanned by the transcript, corresponding to 5'UTR, exonic, intronic, and 3'UTR regions.
        """
        return self._region

    @property
    def exons(self) -> typing.Sequence[GenomicRegion]:
        """
        Get the exon regions.
        """
        return self._exons

    @property
    def cds_start(self) -> typing.Optional[int]:
        """
        Get the 0-based (excluded) start coordinate of the first base of the start codon of the transcript or `None`
        if the transcript is not coding.
        """
        return self._cds_start

    @property
    def cds_end(self) -> typing.Optional[int]:
        """
        Get the 0-based (included) end coordinate of the last base of the termination codon of the transcript or `None`
        if the transcript is not coding.
        """
        return self._cds_end

    def compute_n_codons(self) -> typing.Optional[int]:
        """

        Returns:

        """

    def __eq__(self, other):
        return (isinstance(other, TranscriptCoordinates)
                and self.identifier == other.identifier
                and self.region == other.region
                and self.exons == other.exons
                and self.cds_start == other.cds_start
                and self.cds_end == other.cds_end)

    def __hash__(self):
        return hash((self._id, self._region, self._exons, self._cds_start, self._cds_end))

    def __str__(self):
        return f'TranscriptCoordinates(identifier={self.identifier}, region={self.region})'

    def __repr__(self):
        return str(self)
