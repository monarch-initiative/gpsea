import abc
import typing
import warnings

import hpotk

from .genome import Region, GenomicRegion
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
    """Class that represents results of the functional annotation of a variant with respect to single transcript of a gene.

    Args:
        gene_id (string): The gene symbol associated with the transcript
        tx_id (string): The transcript ID
        hgvsc (string): The HGVS "coding-DNA" ID if available, else None
        is_preferred (bool): The transcript is a MANE transcript, canonical Ensembl transcript, etc.
        variant_effects (Iterable[string]): An iterable of predicted effects given by VEP
        affected_exons (Iterable[integer]): An iterable of exons affected by the variant. Returns None if none are affected.
        protein_id (string): The protein ID for the protein encoded by the transcript.
        protein_effect_coordinates (Region, optional): An optional :class:`Region` with start and end coordinates
         of the effect on the protein sequence.
    """

    def __init__(self, gene_id: str,
                 tx_id: str,
                 hgvs_cdna: typing.Optional[str],
                 is_preferred: bool,
                 variant_effects: typing.Iterable[VariantEffect],
                 affected_exons: typing.Optional[typing.Iterable[int]],
                 protein_id: typing.Optional[str],
                 protein_effect_coordinates: typing.Optional[Region]):
        self._gene_id = hpotk.util.validate_instance(gene_id, str, 'gene_id')
        self._tx_id = hpotk.util.validate_instance(tx_id, str, 'tx_id')
        self._hgvs_cdna = hpotk.util.validate_optional_instance(hgvs_cdna, str, 'hgvs_cdna')
        self._is_preferred = hpotk.util.validate_instance(is_preferred, bool, 'is_preferred')
        self._variant_effects = tuple(variant_effects)
        if affected_exons is not None:
            self._affected_exons = tuple(affected_exons)
        else:
            self._affected_exons = None
        if protein_id is not None:
            self._protein_id = hpotk.util.validate_instance(protein_id, str, 'protein_id')
        else:
            self._protein_id = None
        self._protein_effect_location = hpotk.util.validate_optional_instance(protein_effect_coordinates, Region,
                                                                             'protein_effect_coordinates')

    @property
    def gene_id(self) -> str:
        """
        Returns:
            string: The gene symbol (e.g. SURF1)
        """
        return self._gene_id

    @property
    def transcript_id(self) -> str:
        """
        Returns:
            string: The transcript RefSeq identifier (e.g. NM_123456.7)
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
        Returns:
            string: The HGVS "coding-DNA" representation of the variant (e.g. NM_123456.7:c.9876G>T)
            or `None` if not available.
        """
        return self._hgvs_cdna

    @property
    def variant_effects(self) -> typing.Sequence[VariantEffect]:
        """
        Returns:
            Sequence[string]: A sequence of variant effects.
                    Definitions of these can be found at: http://www.sequenceontology.org/
        """
        return self._variant_effects

    @property
    def overlapping_exons(self) -> typing.Sequence[int]:
        """
        Returns:
            Sequence[integer]: A sequence of 1-based exon indices (the index of the 1st exon is `1`)
            that overlap with the variant.
        """
        return self._affected_exons

    @property
    def protein_id(self) -> typing.Optional[str]:
        """
        Returns:
            Optional[str]: The protein accession ID for the protein relevant to the variant
        """
        return self._protein_id

    @property
    def protein_effect_location(self) -> typing.Optional[Region]:
        """
        Returns:
            Region: a :class:`genophenocorr.model.genome.Region` with start and end position on the protein sequence
            that the variant affects.
        """
        return self._protein_effect_location

    def __str__(self) -> str:
        return f"TranscriptAnnotation(gene_id:{self.gene_id}," \
               f"transcript_id:{self.transcript_id}," \
               f"hgvs_cdna:{self.hgvs_cdna}," \
               f"is_preferred:{self.is_preferred}," \
               f"variant_effects:{self.variant_effects}," \
               f"overlapping_exons:{self.overlapping_exons}," \
               f"protein_id:{self.protein_id}," \
               f"protein_effect_location:{self.protein_effect_location})"

    def __eq__(self, other) -> bool:
        return isinstance(other, TranscriptAnnotation) \
            and self.gene_id == other.gene_id \
            and self.hgvs_cdna == other.hgvs_cdna \
            and self.is_preferred == other.is_preferred \
            and self.transcript_id == other.transcript_id \
            and self.variant_effects == other.variant_effects \
            and self.overlapping_exons == other.overlapping_exons \
            and self.protein_id == other.protein_id \
            and self.protein_effect_location == other.protein_effect_location

    def __repr__(self) -> str:
        return str(self)

    def __hash__(self) -> int:
        return hash((self.gene_id, self.hgvs_cdna, self.is_preferred, self.transcript_id, self.overlapping_exons, self.variant_effects,
                     self.protein_id, self.protein_effect_location))


class VariantCoordinates:
    """A representation of coordinates of sequence and symbolic variants.
    The breakend variants are not supported.

    Attributes:
        region (GenomicRegion): The region spanned by the variant reference allele
        ref (string): The reference allele
        alt (string): The alternate allele
        change_length (integer): The change between the ref and alt alleles due to the variant presence
    """

    def __init__(self, region: GenomicRegion, ref: str, alt: str, change_length: int):
        self._region = hpotk.util.validate_instance(region, GenomicRegion, 'region')
        self._ref = hpotk.util.validate_instance(ref, str, 'ref')
        self._alt = hpotk.util.validate_instance(alt, str, 'alt')
        self._change_length = hpotk.util.validate_instance(change_length, int, 'change_length')

    @property
    def chrom(self) -> str:
        """
        Returns:
            string: The label of the chromosome/contig where the variant is located.
        """
        return self._region.contig.name

    @property
    def start(self) -> int:
        """
        Returns:
            integer: The 0-based start coordinate (excluded) of the ref allele.
        """
        return self._region.start

    @property
    def end(self) -> int:
        """
        Returns:
            integer: The 0-based end coordinate (included) of the ref allele.
        """
        return self._region.end

    @property
    def region(self) -> GenomicRegion:
        """
        Returns:
            GenomicRegion: The genomic region spanned by the ref allele.
        """
        return self._region

    @property
    def ref(self) -> str:
        """
        Returns:
            string: The reference allele (e.g. "A", "N"). The allele may be an empty string.
        """
        return self._ref

    @property
    def alt(self) -> str:
        """
        Returns:
            string: The alternate allele (e.g. "A", "GG", "<DEL>"). The allele may be an empty string for sequence variants.
            The symbolic alternate allele follow the VCF notation and use the `<` and `>` characters
            (e.g. "<DEL>", "<INS:ME:SINE>").
        """
        return self._alt

    @property
    def change_length(self) -> int:
        """
        Get the change of length between the `ref` and `alt` alleles due to the variant presence.

        SNVs lead to change length of zero, deletions and insertions/duplications lead to negative
        and positive change lengths, respectively.
        """
        return self._change_length

    @property
    def variant_key(self) -> str:
        """
        Get a readable representation of the variant's coordinates.

        For instance, ``X_12345_12345_C_G`` for a sequence variant or ``22_10001_20000_INV`` for a symbolic variant.

        .. note::

          Both *start* and *end* coordinates use 1-based (included) coordinate system.
        """
        if self.is_structural():
            return f'{self.chrom}_{self.start + 1}_{self.end}_{self.alt[1:-1]}'
        else:
            key = f'{self.chrom}_{self.start + 1}_{self.end}_{self.ref}_{self.alt}'
            if len(key) > 50:
                ref = None
                alt = None
                if len(self.ref) > 10:
                    ref = f"--{len(self.ref)}bp--"
                if len(self.alt) > 10:
                    alt = f"--{len(self.alt)}bp--"
                return f"{self.chrom}_{self.start + 1}_{self.end}_{ref if not None else self.ref}_{alt if not None else self.alt}"
            else:
                return key

    @property
    def variant_class(self) -> str:
        """
        Get a `str` with the variant class. (e.g. `DUP`, `SNV`, `INS`, `MNV`, `INV`, ...).
        """

        if self.is_structural():
            # Expecting a `str` like <DEL>, <INS>, <DUP>, <INV>, ...
            return self.alt[1:-1]
        else:
            if len(self.ref) > len(self.alt):
                return 'DEL'
            elif len(self.ref) < len(self.alt):
                # may also be a duplication, but it's hard to say from this
                return 'INS'
            else:
                if len(self.ref) == 1:
                    return 'SNV'
                else:
                    return 'MNV'

    def is_structural(self) -> bool:
        """Checks if the variant coordinates use structural variant notation as described by Variant Call Format
        (`VCF <https://en.wikipedia.org/wiki/Variant_Call_Format>`_).

        Ane example of *structural* variant notation::

          chr5  101 . N <DEL> .  .  SVTYPE=DEL;END=120;SVLEN=-10


        as opposed to the *sequence* (literal) notation::

          chr5  101 . NACGTACGTAC N

        :return: `True` if the variant coordinates use structural variant notation.
        """
        return len(self._alt) != 0 and self._alt.startswith('<') and self._alt.endswith('>')

    def __len__(self):
        """
        Get the number of bases on the ref allele that are affected by the variant.
        """
        return len(self._region)

    def __eq__(self, other) -> bool:
        return isinstance(other, VariantCoordinates) \
            and self.region == other.region \
            and self.ref == other.ref \
            and self.alt == other.alt \
            and self.change_length == other.change_length

    def __hash__(self) -> int:
        return hash((self._region, self._ref, self._alt, self._change_length))

    def __str__(self) -> str:
        return f"VariantCoordinates(region={self.region}, " \
               f"ref={self.ref}, alt={self.alt}, " \
               f"change_length={self.change_length})"

    def __repr__(self) -> str:
        return str(self)


class VariantCoordinateAware(metaclass=abc.ABCMeta):

    @property
    @abc.abstractmethod
    def variant_coordinates(self) -> VariantCoordinates:
        pass

    @property
    def variant_string(self) -> str:
        warnings.warn('variant_string` was deprecated and will be removed in v0.2.0. '
                      'Use `variant_coordinates.variant_key` instead', DeprecationWarning, stacklevel=2)
        # TODO[0.2.0] - remove
        return self.variant_coordinates.variant_key

    @property
    def variant_class(self) -> str:
        """
        Returns:
            string: The variant class. (e.g. `DUP`, `SNV`, `INS`, `MNV`, `INV`, ...)
        """
        return self.variant_coordinates.variant_class


class FunctionalAnnotationAware(metaclass=abc.ABCMeta):

    @property
    @abc.abstractmethod
    def tx_annotations(self) -> typing.Sequence[TranscriptAnnotation]:
        pass
    
    def get_tx_anno_by_tx_id(self, transcript_id:str) -> typing.Optional[TranscriptAnnotation]:
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
    
    def get_hgvs_cdna_by_tx_id(self, transcript_id:str) -> typing.Optional[str]:
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
        """A TranscriptAnnotation that represents the result of the functional annotation
        of a variant with respect to the preferred transcript of a gene. Will return None if there is 
        no preferred transcript found. 

        Returns:
            typing.Optional[TranscriptAnnotation]: A TranscriptAnnotation object representing the 
                                                    preferred transcript if one can be found
        """
        for tx in self.tx_annotations:
            if tx.is_preferred:
                return tx
        return None


class Variant(VariantCoordinateAware, FunctionalAnnotationAware, Genotyped):
    """Class that represents results of the functional annotation of a variant with all included transcripts.

    :param var_coordinates: the coordinates of the variant.
    :param tx_annotations: an iterable of functional annotations.
    :param genotypes: the genotypes
    Attributes:
        variant_coordinates (VariantCoordinates):
        variant_string (string): A readable representation of the variant coordinates
        tx_annotations (Sequence[TranscriptAnnotation], Optional): A sequence of TranscriptAnnotation objects representing transcripts affected by this variant
        variant_class (string): The variant class (e.g. Duplication, SNV, etc.)
    """

    @staticmethod
    def create_variant_from_scratch(variant_coordinates: VariantCoordinates,
                                    gene_name: str,
                                    trans_id: str,
                                    hgvs_cdna: str,
                                    is_preferred: bool,
                                    consequences: typing.Iterable[VariantEffect],
                                    exons_effected: typing.Sequence[int],
                                    protein_id: typing.Optional[str],
                                    protein_effect_start: int,
                                    protein_effect_end: int,
                                    genotypes: Genotypes):
        protein_effect = Region(protein_effect_start, protein_effect_end)
        transcript = TranscriptAnnotation(gene_name, trans_id, hgvs_cdna, is_preferred, consequences, exons_effected,
                                          protein_id, protein_effect)
        return Variant(variant_coordinates, (transcript,), genotypes)

    def __init__(self, var_coordinates: VariantCoordinates,
                 tx_annotations: typing.Iterable[TranscriptAnnotation],
                 genotypes: Genotypes):
        """Constructs all necessary attributes for a Variant object

        Args:
            var_coordinates (VariantCoordinates): A VariantCoordinates object with coordinates for this Variant
            tx_annotations (typing.Sequence[TranscriptAnnotation]): A sequence of TranscriptAnnotation objects representing transcripts affected by this variant
            genotypes (Genotypes): genotypes container
        """
        self._var_coordinates = var_coordinates
        self._tx_annotations = tuple(tx_annotations)
        self._gts = genotypes

    @property
    def variant_coordinates(self) -> VariantCoordinates:
        """
        Returns:
            VariantCoordinates: A representation of coordinates of a sequence and symbolic variant.
        """
        return self._var_coordinates

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
        return isinstance(other, Variant) \
            and self.variant_coordinates == other.variant_coordinates \
            and self.tx_annotations == other.tx_annotations \
            and self.genotypes == other.genotypes

    def __hash__(self) -> int:
        return hash(
            (self.variant_coordinates, self.tx_annotations, self.genotypes))

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return (f"Variant(variant_coordinates:{str(self.variant_coordinates)}, "
                f"tx_annotations:{self.tx_annotations}, "
                f"genotypes:{self.genotypes})")

