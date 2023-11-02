import abc
import enum
import typing

import hpotk


class Contig(typing.Sized):
    """
    `Contig` represents identifiers and length of a contiguous sequence of genome assembly.

    The identifiers include:

    * :attr:`name` e.g. `1`
    * :attr:`genbank_acc` e.g. `CM000663.2`
    * :attr:`refseq_name` e.g. `NC_000001.11`
    * :attr:`ucsc_name` e.g. `chr1`

    The length of a `Contig` represents the number of bases of the contig sequence.

    You should not try to create a `Contig` on your own, but always get it from a :class:`GenomeBuild`.
    """

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


class GenomeBuildIdentifier:
    """
    Identifier of the genome build consisting of :attr:`major_assembly` (e.g. GRCh38) and :attr:`patch` (e.g. p13).

    :param major_assembly: major assembly `str`
    :param patch: assembly patch `str`
    """

    def __init__(self, major_assembly: str, patch: str):
        self._major_assembly = major_assembly
        self._patch = patch

    @property
    def major_assembly(self) -> str:
        """
        Get major assembly, e.g. `GRCh38`.
        """
        return self._major_assembly

    @property
    def patch(self) -> str:
        """
        Get assembly patch , e.g. `p13`.
        """
        return self._patch

    @property
    def identifier(self):
        """
        Get genome build identifier consisting of major assembly + patch, e.g. `GRCh38.p13`
        """
        return self._major_assembly + '.' + self._patch

    def __eq__(self, other):
        return (isinstance(other, GenomeBuildIdentifier)
                and self.major_assembly == other.major_assembly
                and self.patch == other.patch)

    def __hash__(self):
        return hash((self.major_assembly, self.patch))

    def __str__(self):
        return f"GenomeBuildIdentifier(prefix={self._major_assembly}, patch={self._patch})"

    def __repr__(self):
        return f"GenomeBuildIdentifier(prefix={self._major_assembly}, patch={self._patch})"


class GenomeBuild:
    """
    `GenomeBuild` is a container for the :attr:`genome_build_id` and the :attr:`contigs` of the build.

    The build supports retrieving contig by various identifiers:

    .. doctest:: genome-build

    >>> from genophenocorr.model.genome import GRCh38

    >>> chr1 = GRCh38.contig_by_name('1')  # by sequence name

    >>> assert chr1 == GRCh38.contig_by_name('CM000663.2')    # by GenBank identifier
    >>> assert chr1 == GRCh38.contig_by_name('NC_000001.11')  # by RefSeq accession
    >>> assert chr1 == GRCh38.contig_by_name('chr1')    # by UCSC name
    """

    def __init__(self, identifier: GenomeBuildIdentifier, contigs: typing.Iterable[Contig]):
        self._id = hpotk.util.validate_instance(identifier, GenomeBuildIdentifier, 'identifier')
        self._contigs = tuple(contigs)
        self._contig_by_name = {}
        for contig in self._contigs:
            self._contig_by_name[contig.name] = contig
            self._contig_by_name[contig.genbank_acc] = contig
            self._contig_by_name[contig.refseq_name] = contig
            self._contig_by_name[contig.ucsc_name] = contig

    @property
    def identifier(self) -> str:
        return self.genome_build_id.identifier

    @property
    def genome_build_id(self) -> GenomeBuildIdentifier:
        return self._id

    @property
    def contigs(self) -> typing.Sequence[Contig]:
        return self._contigs

    def contig_by_name(self, name: str) -> typing.Optional[Contig]:
        """
        Get a contig with `name`.

        The name can come in various formats:

        * sequence name, e.g. `1`
        * GenBank accession, e.g. `CM000663.2`
        * RefSeq accession, e.g. `NC_000001.11`
        * UCSC name, e.g. `chr1`

        :param name: a `str` with contig name.
        """
        try:
            return self._contig_by_name[name]
        except KeyError:
            return None

    def __str__(self):
        return f"GenomeBuild(identifier={self._id.identifier}, n_contigs={len(self.contigs)})"

    def __repr__(self):
        return f"GenomeBuild(identifier={self._id.identifier}, contigs={self.contigs})"


def _a_contains_b(a_start: int, a_end: int, b_start: int, b_end: int) -> bool:
    return a_start <= b_start and b_end <= a_end


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

    def overlaps_with_region(self, other) -> bool:
        """
        Test if this `Region` overlaps with the `other`.

        :param other: another :class:`Region`
        """
        if isinstance(other, Region):
            raise ValueError(f'`other` is not instance of `Region`: {type(other)}')

        return self.overlaps_with(other.start, other.end)

    def overlaps_with(self, start: int, end: int) -> bool:
        """
        Test if this `Region` overlaps with `start` and `end` coordinates of another region.

        .. warning::

          `start` must be at or before `end`. Otherwise, the results are UNDEFINED.

        :param start: 0-based start coordinate of the other region
        :param end: 0-based end coordinate of the other region
        """
        if self.is_empty():
            return _a_contains_b(start, end, self._start, self._end)
        if end - start == 0:
            return _a_contains_b(self._start, self._end, start, end)

        return self.start < end and start < self.end

    def contains_region(self, other) -> bool:
        """
        Test if this `Region` contains the `other` region.

        :param other: another :class:`Region`
        """
        if isinstance(other, Region):
            raise ValueError(f'`other` is not instance of `Region`: {type(other)}')

        return self.contains(other.start, other.end)

    def contains(self, start: int, end: int) -> bool:
        """
        Test if this `Region` contains the another region denoted by `start` and `end` coordinates.

        .. warning::

          `start` must be at or before `end`. Otherwise, the results are UNDEFINED.

        :param start: 0-based start coordinate of the other region
        :param end: 0-based end coordinate of the other region
        """
        return _a_contains_b(self._start, self._end, start, end)

    def is_empty(self) -> bool:
        """
        Return `True` if the region is empty, i.e. it spans 0 units/bases/aminoacids...
        """
        return self._end - self._start == 0

    def __len__(self) -> int:
        return self._end - self._start

    def __eq__(self, other):
        return (isinstance(other, Region)
                and self.start == other.start
                and self.end == other.end)

    def __hash__(self):
        return hash((self._start, self._end))

    def __str__(self):
        return f'Region(start={self.start}, end={self.end})'

    def __repr__(self):
        return str(self)


class Strand(enum.Enum):
    """
    `Strand` is an enum to model positive and negative strands of double-stranded sequences, such as DNA.
    """

    POSITIVE = ('+',)
    """
    The positive strand of a double stranded sequence.
    """

    NEGATIVE = ('-',)
    """
    The negative strand of a double stranded sequence.
    """

    def __init__(self, symbol: str):
        self._symbol = symbol

    @property
    def symbol(self) -> str:
        """
        Get a `str` with the strand's sign.
        """
        return self._symbol

    def is_positive(self):
        """
        Return `True` if this is the *positive* strand.
        """
        return self == Strand.POSITIVE

    def is_negative(self):
        """
        Return `True` if this is the *negative* strand.
        """
        return self == Strand.NEGATIVE

    def opposite(self):
        """
        Get the opposite strand of the current strand.
        """
        return Strand.POSITIVE if self == Strand.NEGATIVE else Strand.NEGATIVE

    def __str__(self):
        return self._symbol


class Stranded(metaclass=abc.ABCMeta):
    """
    Mixin for classes that are on double-stranded sequences.
    """

    @property
    @abc.abstractmethod
    def strand(self) -> Strand:
        pass


class Transposable(Stranded, metaclass=abc.ABCMeta):
    """
    `Transposable` elements know how to flip themselves to arbitrary :class:`Strand` of a sequence.
    """

    @abc.abstractmethod
    def with_strand(self, other: Strand):
        pass

    def to_opposite_strand(self):
        return self.with_strand(self.strand.opposite())

    def to_positive_strand(self):
        return self.with_strand(Strand.POSITIVE)

    def to_negative_strand(self):
        return self.with_strand(Strand.NEGATIVE)


class GenomicRegion(Transposable, Region):
    """
    `GenomicRegion` represents a region located on strand of a DNA contig.

    :param contig: name of the contig, e.g. `15`, `chrX`.
    :param start: 0-based (excluded) start coordinate of the region.
    :param end: 0-based (included) end coordinate of the region.
    :param strand: the strand of the genomic region, `True` for forward strand or `False` for reverse.
    """

    def __init__(self, contig: Contig, start: int, end: int, strand: Strand):
        super().__init__(start, end)
        self._contig = contig
        self._strand = strand
        if end > len(self._contig):
            raise ValueError(f'Genomic region end {end} must not extend '
                             f'beyond contig {self._contig.name} bounds [0,{len(self._contig)}]')

    @property
    def contig(self) -> Contig:
        return self._contig

    def start_on_strand(self, other: Strand) -> int:
        if self.strand == other:
            return self.start
        else:
            return len(self.contig) - self.end

    def end_on_strand(self, other: Strand) -> int:
        if self.strand == other:
            return self.end
        else:
            return len(self.contig) - self.start

    @property
    def strand(self) -> Strand:
        return self._strand

    def with_strand(self, other: Strand):
        if self.strand == other:
            return self
        else:
            return GenomicRegion(self.contig, self.start_on_strand(other), self.end_on_strand(other), other)

    def __eq__(self, other):
        return (isinstance(other, GenomicRegion)
                and self.contig == other.contig
                and self.start == other.strand
                and self.end == other.end
                and self.strand == other.strand)

    def __hash__(self):
        return hash((self.contig, self.start, self.end, self.strand))

    def __str__(self):
        return f'GenomicRegion(contig={self.contig.name}, start={self.start}, end={self.end}, strand={self.strand})'

    def __repr__(self):
        return str(self)
