import abc
import enum
import typing

import hpotk

from .genome import Region, GenomicRegion, Contig, Strand
from ._gt import Genotyped, Genotypes
from ._variant_effects import VariantEffect


class TranscriptInfoAware(metaclass=abc.ABCMeta):
    """
    The implementors know about basic gene/transcript identifiers.
    """

    @property
    @abc.abstractmethod
    def gene_id(self) -> str:
        """
        Returns:
            string: The gene symbol (e.g. SURF1)
        """
        pass

    @property
    @abc.abstractmethod
    def transcript_id(self) -> str:
        """
        Returns:
            string: The transcript RefSeq identifier (e.g. NM_123456.7)
        """
        pass


class TranscriptAnnotation(TranscriptInfoAware):
    """
    `TranscriptAnnotation` represent a result of the functional annotation of a variant
    with respect to single transcript of a gene.
    """

    def __init__(
        self,
        gene_id: str,
        tx_id: str,
        hgvs_cdna: typing.Optional[str],
        is_preferred: bool,
        variant_effects: typing.Iterable[VariantEffect],
        affected_exons: typing.Optional[typing.Iterable[int]],
        protein_id: typing.Optional[str],
        hgvsp: typing.Optional[str],
        protein_effect_coordinates: typing.Optional[Region],
    ):
        self._gene_id = hpotk.util.validate_instance(gene_id, str, "gene_id")
        self._tx_id = hpotk.util.validate_instance(tx_id, str, "tx_id")
        self._hgvs_cdna = hpotk.util.validate_optional_instance(
            hgvs_cdna, str, "hgvs_cdna"
        )
        self._is_preferred = hpotk.util.validate_instance(
            is_preferred, bool, "is_preferred"
        )
        self._variant_effects = tuple(variant_effects)
        if affected_exons is not None:
            self._affected_exons = tuple(affected_exons)
        else:
            self._affected_exons = None
        if protein_id is not None:
            self._protein_id = hpotk.util.validate_instance(
                protein_id, str, "protein_id"
            )
        else:
            self._protein_id = None
        if hgvsp is not None:
            self._hgvsp = hpotk.util.validate_instance(hgvsp, str, "hgvsp")
        else:
            self._hgvsp = None
        self._protein_effect_location = hpotk.util.validate_optional_instance(
            protein_effect_coordinates,
            Region,
            "protein_effect_coordinates",
        )

    @property
    def gene_id(self) -> str:
        """
        Get the HGVS symbol of the affected gene (e.g. *SURF1*).
        """
        return self._gene_id

    @property
    def transcript_id(self) -> str:
        """
        Get the RefSeq identifier of the transcript (e.g. `NM_123456.7`).
        """
        return self._tx_id

    @property
    def is_preferred(self) -> bool:
        """
        Return `True` if the transcript is the preferred transcript of a gene,
        such as MANE transcript, canonical Ensembl transcript.
        """
        return self._is_preferred

    @property
    def hgvs_cdna(self) -> typing.Optional[str]:
        """
        Get the HGVS description of the sequence variant (e.g. ``NM_123456.7:c.9876G>T``)
        or `None` if not available.
        """
        return self._hgvs_cdna

    @property
    def variant_effects(self) -> typing.Sequence[VariantEffect]:
        """
        Get a sequence of the predicted functional variant effects.
        """
        return self._variant_effects

    @property
    def overlapping_exons(self) -> typing.Optional[typing.Sequence[int]]:
        """
        Get a sequence of 1-based exon indices (the index of the 1st exon is `1`)
        that overlap with the variant.
        """
        return self._affected_exons

    @property
    def protein_id(self) -> typing.Optional[str]:
        """
        Get the ID of the protein encoded by the :attr:`transcript_id`
        or `None` if the transcript is not protein-coding.
        """
        return self._protein_id

    @property
    def hgvsp(self) -> typing.Optional[str]:
        """
        Get the HGVS description of the protein sequence variant (e.g. ``NP_001027559.1:p.Gln421Ter``)
        or `None` if not available.
        """
        return self._hgvsp

    @property
    def protein_effect_location(self) -> typing.Optional[Region]:
        """
        Get the :class:`~gpsea.model.genome.Region` with start and end amino-acid coordinates
        affected by the variant.
        """
        return self._protein_effect_location

    def __str__(self) -> str:
        return (
            f"TranscriptAnnotation(gene_id:{self._gene_id},"
            f"transcript_id:{self._tx_id},"
            f"hgvs_cdna:{self._hgvs_cdna},"
            f"is_preferred:{self._is_preferred},"
            f"variant_effects:{self._variant_effects},"
            f"overlapping_exons:{self._affected_exons},"
            f"protein_id:{self._protein_id},"
            f"hgvsp:{self._hgvsp},"
            f"protein_effect_location:{self._protein_effect_location})"
        )

    def __eq__(self, other) -> bool:
        return (
            isinstance(other, TranscriptAnnotation)
            and self._gene_id == other._gene_id
            and self._tx_id == other._tx_id
            and self._hgvs_cdna == other._hgvs_cdna
            and self._is_preferred == other._is_preferred
            and self._variant_effects == other._variant_effects
            and self._affected_exons == other._affected_exons
            and self._protein_id == other._protein_id
            and self._hgvsp == other._hgvsp
            and self._protein_effect_location == other._protein_effect_location
        )

    def __repr__(self) -> str:
        return str(self)

    def __hash__(self) -> int:
        return hash(
            (
                self._gene_id,
                self._tx_id,
                self._hgvs_cdna,
                self._is_preferred,
                self._variant_effects,
                self._affected_exons,
                self._protein_id,
                self._hgvsp,
                self._protein_effect_location,
            )
        )


class VariantClass(enum.Enum):
    """
    `VariantClass` represents a high-level variant category
    which mostly corresponds to the structural variant categories
    of the Variant Call Format specification,
    but includes type for single nucleotide variants (SNV) and multi-nucleotide variant (MNV).
    """

    DEL = 0
    """
    A deletion - a variant with a net loss of sequence from the alternative allele
    regardless of its size.
    
    Both a deletion of 1 base pair and a deletion of 1000 base pairs are acceptable.
    """

    DUP = 1
    """
    Duplication (tandem or interspersed).
    """

    INS = 2
    """
    Insertion of a novel sequence.
    """

    INV = 3
    """
    Inversion of a chromosome segment.
    """

    MNV = 4
    """
    Multi-nucleotide variant (e.g. `AG>CT`) that is not a duplication, deletion, or insertion.
    May be called INDEL.
    """

    SNV = 5
    """
    Single nucleotide variant.
    """

    BND = 6
    """
    A breakend.
    """


class VariantCoordinates:
    """
    A representation of coordinates of sequence and symbolic variants.

    Note, the breakend variants are not currently supported.
    """

    @staticmethod
    def from_vcf_literal(
        contig: Contig,
        pos: int,
        ref: str,
        alt: str,
    ):
        """
        Create `VariantCoordinates` from a variant in VCF literal notation.

        Note, this function must *not* be used to create a VCF symbolic variant
        (e.g. `<DEL>` or translocation).
        Use :func:`from_vcf_symbolic` instead.

        **Example**

        Create a variant from a VCF line:
        ```
        #CHROM  POS     ID  REF ALT ...
        chr1    1001    .   C   T
        ```

        We first must decide on the genome build. Most of the time, we should use GRCh38:

        >>> from gpsea.model.genome import GRCh38
        >>> build = GRCh38

        Then, we access the contig for ``'chr1'``:

        >>> chr1 = build.contig_by_name('chr1')

        Last, we create the variant coordinates:

        >>> from gpsea.model import VariantCoordinates
        >>> vc = VariantCoordinates.from_vcf_literal(
        ...     contig=chr1, pos=1001, ref='C', alt='T',
        ... )

        Now can test the properties:

        >>> vc.start, vc.end, vc.ref, vc.alt, len(vc), vc.change_length
        (1000, 1001, 'C', 'T', 1, 0)

        Args:
            contig: a :class:`Contig` for the chromosome
            pos: a 1-based coordinate of the first base of the reference allele, as described in VCF standard
            ref: a `str` with the REF allele. Should meet the requirements of the VCF standard.
            alt: a `str` with the ALT allele. Should meet the requirements of the VCF standard.
        """
        # TODO: test that we're getting proper alleles
        region = GenomicRegion(
            contig=contig,
            start=pos - 1,
            end=pos + len(ref) - 1,
            strand=Strand.POSITIVE,
        )

        change_length = len(ref) - len(alt)

        return VariantCoordinates(
            region=region,
            ref=ref,
            alt=alt,
            change_length=change_length,
        )

    @staticmethod
    def from_vcf_symbolic(
        contig: Contig,
        pos: int,
        end: int,
        ref: str,
        alt: str,
        svlen: int,
    ):
        """
        Create `VariantCoordinates` from a variant in VCF symbolic notation.

        Note, this function must *not* be used to create a VCF sequence/literal variant.
        Use :func:`from_vcf_literal` instead.

        **Example**

        Let's create a symbolic variant from the line:

        ```
        #CHROM   POS      ID   REF   ALT     QUAL   FILTER   INFO
        2        321682   .    T     <DEL>   6      PASS     SVTYPE=DEL;END=321887;SVLEN=-205
        ```

        We first must decide on the genome build. Most of the time, we should use GRCh38:

        >>> from gpsea.model.genome import GRCh38
        >>> contig = GRCh38.contig_by_name('2')

        Now, we create the coordinates as:

        >>> vc = VariantCoordinates.from_vcf_symbolic(
        ...     contig=contig, pos=321682, end=321887,
        ...     ref='T', alt='<DEL>', svlen=-205,
        ... )

        Now can test the properties:

        >>> vc.start, vc.end, vc.ref, vc.alt, len(vc), vc.change_length
        (321681, 321887, 'T', '<DEL>', 206, -205)

        Args:
            contig: a :class:`Contig` for the chromosome
            pos: a 1-based coordinate of the first base of the affected reference allele region
            end: a 1-based coordinate of the last base of the affected reference allele region
            ref: a `str` with the REF allele. Most of the time, it is one of `{'N', 'A', 'C', 'G', 'T'}`
            alt: a `str` with the ALT allele, e.g. one of `{'<DEL>', '<DUP>', '<INS>', '<INV>'}`
            svlen: an `int` with change length (the difference between ref and alt allele lengths)
        """
        # TODO: test that we're getting proper alleles
        region = GenomicRegion(
            contig=contig,
            start=pos - 1,  # convert to 0-based coordinates,
            end=end,
            strand=Strand.POSITIVE,
        )

        return VariantCoordinates(
            region=region,
            ref=ref,
            alt=alt,
            change_length=svlen,
        )

    def __init__(self, region: GenomicRegion, ref: str, alt: str, change_length: int):
        self._region = hpotk.util.validate_instance(region, GenomicRegion, "region")
        self._ref = hpotk.util.validate_instance(ref, str, "ref")
        self._alt = hpotk.util.validate_instance(alt, str, "alt")
        self._change_length = hpotk.util.validate_instance(
            change_length, int, "change_length"
        )

    @property
    def chrom(self) -> str:
        """
        Get the label of the chromosome/contig where the variant is located.
        """
        return self._region.contig.name

    @property
    def start(self) -> int:
        """
        Get the 0-based start coordinate (excluded) of the first base of the :attr:`ref` allele.
        """
        return self._region.start

    @property
    def end(self) -> int:
        """
        Get the 0-based end coordinate (included) of the last base of the :attr:`ref` allele.
        """
        return self._region.end

    @property
    def region(self) -> GenomicRegion:
        """
        Get the genomic region spanned by the :attr:`ref` allele.
        """
        return self._region

    @property
    def ref(self) -> str:
        """
        Get the reference allele (e.g. "A", "CCT", "N"). The allele may be an empty string.
        """
        return self._ref

    @property
    def alt(self) -> str:
        """
        Get the alternate allele (e.g. "A", "GG", "<DEL>").

        The allele may be an empty string for sequence variants.
        The symbolic alternate allele follow the VCF notation and use the `<` and `>` characters
        (e.g. "<DEL>", "<INS:ME:SINE>").
        """
        return self._alt

    @property
    def change_length(self) -> int:
        """
        Get the change of length between the `ref` and `alt` alleles due to the variant presence.

        See :ref:`change-length-of-an-allele` for more info.
        """
        return self._change_length

    @property
    def variant_key(self) -> str:
        """
        Get a readable representation of the variant's coordinates.

        For instance, ``X_12345_12345_C_G`` for a sequence variant or ``22_10001_20000_INV`` for a symbolic variant.
        If the key is larger than 50 characters, the 'ref' and/or 'alt' (if over 10 bps)
        are changed to just show number of bps.
        Example: ``X_1000001_1000027_TAAAAAAAAAAAAAAAAAAAAAAAAAA_T`` -> ``X_1000001_1000027_--27bp--_T``

        .. note::

          Both *start* and *end* coordinates use 1-based (included) coordinate system.
        """
        if self.is_structural():
            return f"{self.chrom}_{self.start + 1}_{self.end}_{self.alt[1:-1]}"
        else:
            key = f"{self.chrom}_{self.start + 1}_{self.end}_{self.ref}_{self.alt}"
            if len(key) > 50:
                if len(self.ref) > 10:
                    ref = f"--{len(self.ref)}bp--"
                else:
                    ref = self.ref
                if len(self.alt) > 10:
                    alt = f"--{len(self.alt)}bp--"
                else:
                    alt = self.alt
                return f"{self.chrom}_{self.start + 1}_{self.end}_{ref}_{alt}"
            else:
                return key

    @property
    def variant_class(self) -> VariantClass:
        """
        Get a :class:`VariantClass` category.
        """
        if self.is_structural():
            # Expecting a `str` like <DEL>, <INS>, <DUP>, <INV>, ...
            return VariantClass[self.alt[1:-1]]
        else:
            if len(self.ref) > len(self.alt):
                return VariantClass.DEL
            elif len(self.ref) < len(self.alt):
                # may also be a duplication, but it's hard to say from this
                return VariantClass.INS
            else:
                if len(self.ref) == 1:
                    return VariantClass.SNV
                else:
                    return VariantClass.MNV

    def is_structural(self) -> bool:
        """
        Checks if the variant coordinates use structural variant notation as described by Variant Call Format
        (`VCF <https://en.wikipedia.org/wiki/Variant_Call_Format>`_).

        Ane example of *structural* variant notation::

          chr5  101 . N <DEL> .  .  SVTYPE=DEL;END=120;SVLEN=-10


        as opposed to the *sequence* (literal) notation::

          chr5  101 . NACGTACGTAC N

        :return: `True` if the variant coordinates use structural variant notation.
        """
        return (
            len(self._alt) != 0
            and self._alt.startswith("<")
            and self._alt.endswith(">")
        )

    def __len__(self):
        """
        Get the number of bases on the ref allele that are affected by the variant.
        """
        return len(self._region)

    def __eq__(self, other) -> bool:
        return (
            isinstance(other, VariantCoordinates)
            and self.region == other.region
            and self.ref == other.ref
            and self.alt == other.alt
            and self.change_length == other.change_length
        )

    def __hash__(self) -> int:
        return hash((self._region, self._ref, self._alt, self._change_length))

    def __str__(self) -> str:
        return (
            f"VariantCoordinates(region={self.region}, "
            f"ref={self.ref}, alt={self.alt}, "
            f"change_length={self.change_length})"
        )

    def __repr__(self) -> str:
        return str(self)


class ImpreciseSvInfo:
    """
    Data regarding a structural variant (SV) with imprecise breakpoint coordinates.
    """

    def __init__(
        self,
        structural_type: hpotk.TermId,
        variant_class: VariantClass,
        gene_id: str,
        gene_symbol: str,
    ):
        self._structural_type = hpotk.util.validate_instance(
            structural_type, hpotk.TermId, "structural_type"
        )
        self._variant_class = hpotk.util.validate_instance(
            variant_class, VariantClass, "variant_class"
        )
        self._gene_id = gene_id
        self._gene_symbol = gene_symbol

    @property
    def structural_type(self) -> hpotk.TermId:
        """
        Get term ID of the structural type (e.g. ``SO:1000029`` for chromosomal deletion).
        """
        return self._structural_type

    @property
    def variant_class(self) -> VariantClass:
        """
        Get a :class:`VariantClass` category.
        """
        return self._variant_class

    @property
    def gene_id(self) -> str:
        """
        Get a `str` with gene identifier CURIE (e.g. ``HGNC:3603``) or `None` if the identifier is not available.
        """
        return self._gene_id

    @property
    def gene_symbol(self) -> str:
        """
        Get a `str` with HGVS gene symbol (e.g. *FBN1*) or `None` if the symbol is not available.
        """
        return self._gene_symbol

    @property
    def variant_key(self) -> str:
        """
        Get a readable representation of the variant.
        """
        return f"{self._structural_type}_{self._gene_id}_{self._gene_symbol}"

    def __eq__(self, value: object) -> bool:
        return (
            isinstance(value, ImpreciseSvInfo)
            and self._structural_type == value._structural_type
            and self._variant_class == value._variant_class
            and self._gene_id == value._gene_id
            and self._gene_symbol == value._gene_symbol
        )

    def __hash__(self) -> int:
        return hash(
            (
                self._structural_type,
                self._variant_class,
                self._gene_id,
                self._gene_symbol,
            )
        )

    def __str__(self) -> str:
        return (
            f"ImpreciseSv("
            f"structural_type={self._structural_type}, "
            f"variant_class={self._variant_class}, "
            f"gene_id={self._gene_id}, "
            f"gene_symbol={self._gene_symbol})"
        )

    def __repr__(self) -> str:
        return str(self)


class VariantInfo:
    """
    `VariantInfo` consists of either variant coordinates or imprecise SV data.

    The class is conceptually similar to Rust enum - only one of the fields can be set at any point in time.
    """

    def __init__(
        self,
        variant_coordinates: typing.Optional[VariantCoordinates] = None,
        sv_info: typing.Optional[ImpreciseSvInfo] = None,
    ):
        if (variant_coordinates is None) != (sv_info is None):
            # At most one field can be set.
            self._variant_coordinates = variant_coordinates
            self._sv_info = sv_info
        else:
            raise ValueError(
                f"Only one field can be set: variant_coordinates={variant_coordinates}, sv_info={sv_info}"
            )

    @property
    def variant_coordinates(self) -> typing.Optional[VariantCoordinates]:
        """
        Get variant coordinates if available.
        """
        return self._variant_coordinates

    def has_variant_coordinates(self) -> bool:
        """
        Returns `True` if the variant coordinates are available.
        """
        return self.variant_coordinates is not None

    @property
    def sv_info(self) -> typing.Optional[ImpreciseSvInfo]:
        """
        Get information about large imprecise SV.
        """
        return self._sv_info

    def has_sv_info(self) -> bool:
        """
        Returns `True` if the variant is a large imprecise SV and the exact coordinates are thus unavailable.
        """
        return self.sv_info is not None

    @property
    def variant_key(self) -> str:
        """
        Get a readable representation of the variant's coordinates or the large SV info.
        """
        if self.has_variant_coordinates():
            return self.variant_coordinates.variant_key
        elif self.has_sv_info():
            return self.sv_info.variant_key
        else:
            self._handle_missing_state()

    @property
    def variant_class(self) -> VariantClass:
        """
        Get a :class:`VariantClass` category.
        """
        if self.has_variant_coordinates():
            return self.variant_coordinates.variant_class
        elif self.has_sv_info():
            return self.sv_info.variant_class
        else:
            self._handle_missing_state()

    def is_structural(self) -> bool:
        """
        Test if the variant is a structural variant.

        This can either be because the variant coordinates
        use the structural variant notation (see :meth:`VariantCoordinates.is_structural`)
        or if the variant is large imprecise SV.
        """
        if self.has_variant_coordinates():
            return self.variant_coordinates.is_structural()
        elif self.has_sv_info():
            return True
        else:
            self._handle_missing_state()

    def _handle_missing_state(self):
        raise ValueError(
            "VariantInfo should have either variant coordinates or SV info!"
        )

    def __eq__(self, value: object) -> bool:
        return (
            isinstance(value, VariantInfo)
            and self._variant_coordinates == value._variant_coordinates
            and self._sv_info == value._sv_info
        )

    def __hash__(self) -> int:
        return hash((self._variant_coordinates, self._sv_info))

    def __str__(self) -> str:
        return (
            f"VariantInfo("
            f"variant_coordinates={self._variant_coordinates}, "
            f"sv_info={self._sv_info})"
        )

    def __repr__(self) -> str:
        return str(self)


class VariantInfoAware(metaclass=abc.ABCMeta):
    """
    An entity where :class:`VariantInfo` is available.
    """

    @property
    @abc.abstractmethod
    def variant_info(self) -> VariantInfo:
        """
        Get the variant data with coordinates or other info available for large imprecise SVs.
        """
        pass


class FunctionalAnnotationAware(metaclass=abc.ABCMeta):

    @property
    @abc.abstractmethod
    def tx_annotations(self) -> typing.Sequence[TranscriptAnnotation]:
        pass

    def get_tx_anno_by_tx_id(
        self, transcript_id: str
    ) -> typing.Optional[TranscriptAnnotation]:
        """Given a transcript ID, this will return the `TranscriptAnnotation` associated with that
        variant and transcript.

        Args:
            transcript_id (str): A transcript ID - i.e. 'NM_170707.4'

        Returns:
            typing.Optional[TranscriptAnnotation]: The Transcript Annotation if available.
        """
        for tx_ann in self.tx_annotations:
            if tx_ann.transcript_id == transcript_id:
                return tx_ann
        return None

    def get_hgvs_cdna_by_tx_id(self, transcript_id: str) -> typing.Optional[str]:
        """Given a transcript ID, will return the hgvs cdna string associated with that variant and transcript.

        Args:
            transcript_id (str): A transcript ID - i.e. 'NM_170707.4'

        Returns:
            str or None: The hgvs cdna if available - i.e. 'NM_170707.4:c.1824C>T'
        """
        for tx_ann in self.tx_annotations:
            if tx_ann.transcript_id == transcript_id:
                return tx_ann.hgvs_cdna
        return None

    def get_preferred_tx_annotation(self) -> typing.Optional[TranscriptAnnotation]:
        """Get the `TranscriptAnnotation` that represents the result of the functional annotation
        with respect to the preferred transcript of a gene.

        Returns `None` if transcript annotations is no preferred transcript found.

        Returns:
            typing.Optional[TranscriptAnnotation]: The `TranscriptAnnotation` with respect
                                                   to the preferred transcript
                                                   or `None` if the preferred transcript info is not available.
        """
        for tx in self.tx_annotations:
            if tx.is_preferred:
                return tx
        return None


class Variant(VariantInfoAware, FunctionalAnnotationAware, Genotyped):
    """
    Variant includes three lines of information:
     * the variant data with coordinates or other info available for large imprecise SVs,
     * results of the functional annotation with respect to relevant transcripts, and
     * the genotypes for the known samples
    """

    @staticmethod
    def create_variant_from_scratch(
        variant_coordinates: VariantCoordinates,
        gene_name: str,
        trans_id: str,
        hgvs_cdna: str,
        is_preferred: bool,
        consequences: typing.Iterable[VariantEffect],
        exons_effected: typing.Sequence[int],
        protein_id: typing.Optional[str],
        hgvsp: typing.Optional[str],
        protein_effect_start: typing.Optional[int],
        protein_effect_end: typing.Optional[int],
        genotypes: Genotypes,
    ):
        variant_info = VariantInfo(variant_coordinates=variant_coordinates)
        if protein_effect_start is not None and protein_effect_end is not None:
            protein_effect = Region(protein_effect_start, protein_effect_end)
        else:
            protein_effect = None
        transcript = TranscriptAnnotation(
            gene_name,
            trans_id,
            hgvs_cdna,
            is_preferred,
            consequences,
            exons_effected,
            protein_id,
            hgvsp,
            protein_effect,
        )
        return Variant(variant_info, (transcript,), genotypes)

    def __init__(
        self,
        variant_info: VariantInfo,
        tx_annotations: typing.Iterable[TranscriptAnnotation],
        genotypes: Genotypes,
    ):
        self._variant_info = variant_info
        self._tx_annotations = tuple(tx_annotations)
        self._gts = genotypes

    @property
    def variant_info(self) -> VariantInfo:
        """
        Get the representation of the variant data for sequence and symbolic variants, 
        as well as for large imprecise SVs.
        """
        return self._variant_info

    @property
    def tx_annotations(self) -> typing.Sequence[TranscriptAnnotation]:
        """A collection of TranscriptAnnotations that each represent results of the functional annotation
        of a variant with respect to single transcript of a gene.

        Returns:
            Sequence[TranscriptAnnotation]: A sequence of TranscriptAnnotation objects
        """
        return self._tx_annotations

    @property
    def genotypes(self) -> Genotypes:
        return self._gts

    def __eq__(self, other) -> bool:
        return (
            isinstance(other, Variant)
            and self._variant_info == other._variant_info
            and self._tx_annotations == other._tx_annotations
            and self._gts == other._gts
        )

    def __hash__(self) -> int:
        return hash((self._variant_info, self._tx_annotations, self._gts))

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return (
            f"Variant(variant_info={self._variant_info}, "
            f"tx_annotations={self._tx_annotations}, "
            f"genotypes={self._gts})"
        )
