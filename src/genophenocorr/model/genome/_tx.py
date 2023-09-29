import typing

import hpotk

from ._genome import GenomicRegion


class TranscriptCoordinates:
    """
    `TranscriptCoordinates` knows about genomic region of the transcript, exonic/intronic regions, as well as
    the coding and non-coding regions.
    """

    def __init__(self, identifier: str,
                 region: GenomicRegion,
                 exons: typing.Iterable[GenomicRegion],
                 cds_start: int,
                 cds_end: int):
        self._id = hpotk.util.validate_instance(identifier, str, 'identifier')
        self._region = hpotk.util.validate_instance(region, GenomicRegion, 'region')
        self._exons = tuple(exons)
        # TODO - check
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
    def cds_start(self) -> int:
        return self._cds_start

    @property
    def cds_end(self) -> int:
        return self._cds_end

    def __eq__(self, other):
        return (isinstance(other, TranscriptCoordinates)
                and self.identifier == other.identifier
                and self.region == other.region
                and self.exons == other.exons
                and self.cds_start == other.cds_start
                and self.cds_end == other.cds_end)

    def __str__(self):
        return f'TranscriptCoordinates(identifier={self.identifier}, region={self.region})'

    def __repr__(self):
        return str(self)
    