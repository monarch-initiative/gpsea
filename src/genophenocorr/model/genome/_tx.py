import typing

import hpotk

from ._genome import GenomicRegion


class Transcript:

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

    def __str__(self):
        return f'Transcript(identifier={self.identifier}, region={self.region})'

    def __repr__(self):
        return str(self)
    