import typing

import hpotk

from .genome import GenomicRegion


class TranscriptCoordinates:
    """
    `TranscriptCoordinates` knows about genomic region of the transcript, exonic/intronic regions, as well as
    the coding and non-coding regions.

    If both CDS coordinates are `None`, then the transcript coordinates are assumed to represent a non-coding transcript.
    """

    def __init__(
        self, 
        identifier: str,
        region: GenomicRegion,
        exons: typing.Iterable[GenomicRegion],
        cds_start: typing.Optional[int],
        cds_end: typing.Optional[int],
        is_preferred: typing.Optional[bool] = None
    ):
        self._id = hpotk.util.validate_instance(identifier, str, 'identifier')
        self._region = hpotk.util.validate_instance(region, GenomicRegion, 'region')
        self._exons = tuple(exons)

        if cds_start is None and cds_end is None:
            # non-coding transcript, no validation necessary
            pass
        else:
            # coding transcript
            if cds_start is None or cds_end is None:
                raise ValueError(f'Both cds_start={cds_start} and cds_end={cds_end} must not be `None`')
            if not isinstance(cds_start, int) or not isinstance(cds_end, int):
                raise ValueError(f'CDS coordinates must be `int`s but were cds_start={cds_start}, cds_end={cds_end}')
            if not self._region.contains_pos(cds_start + 1) or not self._region.contains_pos(cds_end):
                # `+1` to convert the 0-based start into a 1-based coordinate for a moment
                raise ValueError(f'CDS coordinates cds_start={cds_start:,}, cds_end={cds_end:,} must be '
                                 f'in the tx region: ({self._region.end}, {self._region.end}]')
        self._cds_start = cds_start
        self._cds_end = cds_end
        self._is_preferred = is_preferred

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

    def is_coding(self) -> bool:
        """
        Return `True` if the transcript is coding/translated.
        """
        return self._cds_start is not None and self._cds_end is not None

    def get_coding_base_count(self) -> typing.Optional[int]:
        """
        Calculate the number of coding bases present in the transcript. Note, the count does *NOT* include
        the termination codon since it does not code for an aminoacid.
        Returns: an `int` with the coding base count or `None` if the transcript is non-coding.
        """
        if not self.is_coding():
            return None

        n_bases = 0
        for exon in self.exons:
            start = max(self._cds_start, exon.start)
            end = min(self._cds_end, exon.end)
            n_bases += max(end - start, 0)

        return n_bases - 3  # minus stop codon

    def get_codon_count(self) -> typing.Optional[int]:
        """
        Calculate the count of codons present in the transcript. Note, the count does *NOT* include the termination codon!

        Returns: the number of codons of the transcript or `None` if the transcript is non-coding.
        """
        n_coding_bases = self.get_coding_base_count()
        assert n_coding_bases % 3 == 0, (f"Transcript {self._id} has {n_coding_bases:,} "
                                         f"coding bases that is not divisible by 3!")
        return int(n_coding_bases / 3)

    def get_five_prime_utrs(self) -> typing.Sequence[GenomicRegion]:
        """
        Get a sequence of genomic regions that correspond to 5' untranslated regions of the transcript.

        Returns: a sequence of genomic regions, an empty sequence if the transcript is non-coding.
        """
        if not self.is_coding():
            return ()

        utrs = []

        for exon in self.exons:
            if exon.start >= self._cds_start:
                break
            if exon.end <= self._cds_start:
                # An entire exon is UTR.
                utrs.append(exon)
            else:
                # Just a part of an exon is UTR.
                utrs.append(GenomicRegion(exon.contig, exon.start, self._cds_start, exon.strand))

        return tuple(utrs)

    def get_three_prime_utrs(self) -> typing.Sequence[GenomicRegion]:
        """
        Get a sequence of genomic regions that correspond to 3' untranslated regions of the transcript.

        Note, the termination codon is *NOT* included in the regions!

        Returns: a sequence of genomic regions, an empty sequence if the transcript is non-coding.
        """
        if not self.is_coding():
            return ()

        utrs = []

        for exon in self.exons:
            if exon.end > self._cds_end:
                if self._cds_end <= exon.start:
                    # An entire exon is UTR.
                    utrs.append(exon)
                else:
                    # Just a part of an exon is UTR.
                    utrs.append(GenomicRegion(exon.contig, self._cds_end, exon.end, exon.strand))

        return tuple(utrs)

    def get_cds_regions(self) -> typing.Sequence[GenomicRegion]:
        """
        Get a sequence of genomic regions that correspond to coding regions of the transcript, including
        BOTH the initiation and termination codons.

        Returns: a sequence of genomic regions, an empty sequence if the transcript is non-coding.
        """
        if not self.is_coding():
            return ()

        coding_regions = []
        for exon in self.exons:
            if self._cds_start >= exon.end or self._cds_end <= exon.start:
                # An entire exon is UTR.
                continue
            else:
                cds_start_is_at_or_before_start = self._cds_start <= exon.start
                cds_end_is_at_or_after_end = exon.end <= self._cds_end
                if cds_start_is_at_or_before_start or cds_end_is_at_or_after_end:
                    # At least one base of the exon must be coding.
                    if cds_start_is_at_or_before_start and cds_end_is_at_or_after_end:
                        # The entire exon is coding
                        coding_regions.append(exon)
                    else:
                        # Part of the exon is coding, another part is UTR
                        cds = GenomicRegion(exon.contig,
                                            max(self._cds_start, exon.start),
                                            min(exon.end, self._cds_end),
                                            exon.strand)
                        coding_regions.append(cds)

        return tuple(coding_regions)
    
    @property
    def is_preferred(self) -> typing.Optional[bool]:
        """
        Check if the transcript belongs among the preferred transcripts of the gene. 
        This usually means that the transcript was chosen by the MANE project or similar.
        
        Returns `None` if the info is not available.
        """
        return self._is_preferred

    def __eq__(self, other):
        return (isinstance(other, TranscriptCoordinates)
                and self._id == other._id
                and self._region == other._region
                and self._exons == other._exons
                and self._cds_start == other._cds_start
                and self._cds_end == other._cds_end
                and self._is_preferred == other._is_preferred)

    def __hash__(self):
        return hash((self._id, self._region, self._exons, self._cds_start, self._cds_end, self._is_preferred))

    def __str__(self):
        return f'TranscriptCoordinates(identifier={self.identifier}, region={self.region})'

    def __repr__(self):
        return str(self)
