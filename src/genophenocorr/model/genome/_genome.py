import abc
import typing

import hpotk


class Contig(typing.Sized):

    def __init__(self, name: str, gb_acc: str, refseq_name: str, ucsc_name: str, length: int):
        self._name = hpotk.util.validate_instance(name, str, 'name')

        self._gb_acc = hpotk.util.validate_instance(gb_acc, str, 'gb_acc')
        self._refseq = hpotk.util.validate_instance(refseq_name, str, 'refseq_name')
        self._ucsc = hpotk.util.validate_instance(ucsc_name, str, 'ucsc_name')

        self._len = hpotk.util.validate_instance(length, int, 'length')
        if self._len < 0:
            raise ValueError(f'Length must not be negative but got {self._len}')

    @property
    def name(self) -> str:
        return self._name

    @property
    def genbank_acc(self) -> str:
        return self._gb_acc

    @property
    def refseq_name(self) -> str:
        return self._refseq

    @property
    def ucsc_name(self) -> str:
        return self._ucsc

    def __len__(self) -> int:
        return self._len

    def __str__(self):
        return f"Contig(name={self.name}, refseq_name={self.refseq_name}, ucsc_name={self.ucsc_name}, len={self._len})"

    def __repr__(self):
        return str(self)

    def __eq__(self, other):
        return (isinstance(other, Contig)
                and self.name == other.name
                and self.refseq_name == other.refseq_name
                and self.ucsc_name == other.ucsc_name
                and len(self) == len(other))

    def __hash__(self):
        return hash((self.name, self.refseq_name, self.ucsc_name, len(self)))


class GenomeBuild:

    def __init__(self, identifier: str, contigs: typing.Iterable[Contig]):
        self._id = identifier
        self._contigs = tuple(contigs)
        self._contig_by_name = {}
        for contig in self._contigs:
            self._contig_by_name[contig.name] = contig
            self._contig_by_name[contig.genbank_acc] = contig
            self._contig_by_name[contig.refseq_name] = contig
            self._contig_by_name[contig.ucsc_name] = contig

    @property
    def identifier(self) -> str:
        return self._id

    @property
    def contigs(self) -> typing.Sequence[Contig]:
        return self._contigs

    def contig_by_name(self, name: str) -> typing.Optional[Contig]:
        try:
            return self._contig_by_name[name]
        except KeyError:
            return None

    def __str__(self):
        return f"GenomeBuild(identifier={self.identifier}, n_contigs={len(self.contigs)})"

    def __repr__(self):
        return f"GenomeBuild(identifier={self.identifier}, contigs={self.contigs})"


class Region(typing.Sized):
    """
    `Region` represents a contiguous region/slice of a biological sequence, such as DNA, RNA, or protein.

    The :attr:`start` and :attr:`end` represent 0-based coordinates of the region.
    The region has length that corresponds to the number of spanned bases/aminoacids.

    :param start: 0-based (excluded) start coordinate of the region.
    :param end: 0-based (included) end coordinate of the region.
    """

    def __init__(self, start: int, end: int):
        if not isinstance(start, int) or not isinstance(end, int):
            raise ValueError(f'`start` and `end` must be ints but were `{type(start)}`, `{type(end)}`')
        if start < 0 or end < 0:
            raise ValueError(f'`start` and `end` must be positive but were `{start}`, `{end}`')
        if start > end:
            raise ValueError(f'`start` {start} must be at or before `end` {end}')
        self._start = start
        self._end = end

    @property
    def start(self) -> int:
        """
        Get 0-based (excluded) start coordinate of the region.
        """
        return self._start

    @property
    def end(self) -> int:
        """
        Get 0-based (included) end coordinate of the region.
        """
        return self._end

    def __len__(self) -> int:
        return self._end - self._start

    def __eq__(self, other):
        return (isinstance(other, Region)
                and self.start == other.start
                and self.end == other.end)

    def __str__(self):
        return f'Region(start={self.start}, end={self.end})'

    def __repr__(self):
        return str(self)


class Stranded(metaclass=abc.ABCMeta):
    """
    Mixin for classes that are on double-stranded sequence.
    """

    @property
    @abc.abstractmethod
    def strand(self) -> bool:
        pass

    def is_forward_strand(self) -> bool:
        return self.strand

    def is_reverse_strand(self) -> bool:
        return not self.strand


class GenomicRegion(Stranded, Region):
    """
    `GenomicRegion` represents a region located on strand of a DNA contig.

    :param contig: name of the contig, e.g. `15`, `chrX`.
    :param start: 0-based (excluded) start coordinate of the region.
    :param end: 0-based (included) end coordinate of the region.
    :param strand: the strand of the genomic region, `True` for forward strand or `False` for reverse.
    """

    def __init__(self, contig: Contig, start: int, end: int, strand: bool):
        super().__init__(start, end)
        self._contig = contig
        self._strand = strand

    @property
    def contig(self) -> Contig:
        return self._contig

    def start_on_strand(self, strand: bool) -> int:
        if self.strand == strand:
            return self.start
        else:
            return len(self.contig) - self.end

    def end_on_strand(self, strand: bool) -> int:
        if self.strand == strand:
            return self.end
        else:
            return len(self.contig) - self.start

    @property
    def strand(self) -> bool:
        return self._strand

    def __str__(self):
        return f'GenomicRegion(contig={self.contig.name}, start={self.start}, end={self.end}, strand={self.strand})'

    def __repr__(self):
        return str(self)
