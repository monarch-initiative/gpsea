import abc
import typing
import warnings

import hpotk

from .genome import GenomicRegion
from ._gt import Genotyped, Genotypes
from ._protein import ProteinMetadata


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

    Attributes:
        gene_id (string): The gene symbol associated with the transcript
        transcript_id (string): The transcript ID
        hgvsc_id (string): The HGVS "coding-DNA" ID if available, else None
        is_preferred (bool): The transcript is a MANE transcript, canonical Ensembl transcript, etc.
        variant_effects (Sequence[string]): A sequence of predicted effects given by VEP
        overlapping_exons (Sequence[integer]): A sequence of exons affected by the variant. Returns None if none are affected.
        protein_affected (ProteinMetadata): A ProteinMetadata object representing the protein affected by this transcript
        protein_effect_location (Tuple(integer, integer)): The start and end coordinates of the effect on the protein sequence.
    """

    def __init__(self, gene_id: str,
                 tx_id: str,
                 hgvsc: typing.Optional[str],
                 is_preferred: bool,
                 variant_effects,
                 affected_exons: typing.Optional[typing.Sequence[int]],
                 affected_protein: typing.Sequence[ProteinMetadata],
                 protein_effect_start: typing.Optional[int],
                 protein_effect_end: typing.Optional[int]):
        """Constructs all necessary attributes for a TranscriptAnnotation object

        Args:
            gene_id (string): The gene symbol associated with the transcript
            tx_id (string): The transcript ID
            hgvsc (string, Optional): The HGVS "coding-DNA" ID if available, else None
            variant_effects (Sequence[string]): A sequence of predicted effects given by VEP
            affected_exons (Sequence[integer], Optional): A sequence of exons affected by the variant. Returns None if none are affected.
            affected_protein (Sequence[ProteinMetadata]): A ProteinMetadata object representing the protein affected by this transcript
            protein_effect_start (integer, Optional): The start coordinate of the effect on the protein sequence.
            protein_effect_end (integer, Optional): The end coordinate of the effect on the protein sequence.
        """
        self._gene_id = gene_id
        self._tx_id = tx_id
        self._hgvsc_id = hgvsc
        self._is_preferred = is_preferred
        self._variant_effects = tuple(variant_effects)
        if affected_exons is not None:
            self._affected_exons = tuple(affected_exons)
        else:
            self._affected_exons = None
        self._affected_protein = tuple(affected_protein)
        self._protein_effect_location = (protein_effect_start, protein_effect_end)

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
    def hgvsc_id(self) -> typing.Optional[str]:
        """
        Returns:
            string: The HGVS "coding-DNA" representation of the variant (e.g. NM_123456.7:c.9876G>T)
            or `None` if not available.
        """
        return self._hgvsc_id

    @property
    def variant_effects(self) -> typing.Sequence[str]:
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
    def protein_affected(self) -> typing.Sequence[ProteinMetadata]:
        """
        Returns:
            Sequence[ProteinMetadata]: The ProteinMetadata object representing the protein that is affected by the alteration of this transcript
        """
        return self._affected_protein

    @property
    def protein_effect_location(self) -> typing.Tuple[int, int]:
        """
        Returns:
            Tuple(integer, integer): The start and end position on the protein sequence that the variant effects. (e.g. (1234, 1235))
        """
        return self._protein_effect_location

    def __str__(self) -> str:
        return f"TranscriptAnnotation(gene_id:{self.gene_id}," \
               f"transcript_id:{self.transcript_id}," \
               f"hgvsc_id:{self.hgvsc_id}," \
               f"is_preferred:{self.is_preferred}," \
               f"variant_effects:{self.variant_effects}," \
               f"overlapping_exons:{self.overlapping_exons}," \
               f"protein_affected:{self.protein_affected}," \
               f"protein_effect_location:{self.protein_effect_location})"

    def __eq__(self, other) -> bool:
        return isinstance(other, TranscriptAnnotation) \
            and self.gene_id == other.gene_id \
            and self.hgvsc_id == other.hgvsc_id \
            and self.is_preferred == other.is_preferred \
            and self.transcript_id == other.transcript_id \
            and self.variant_effects == other.variant_effects \
            and self.overlapping_exons == other.overlapping_exons \
            and self.protein_affected == other.protein_affected \
            and self.protein_effect_location == other.protein_effect_location

    def __repr__(self) -> str:
        return str(self)

    def __hash__(self) -> int:
        return hash((self.gene_id, self.hgvsc_id, self.is_preferred, self.transcript_id, self.overlapping_exons, self.variant_effects,
                     self.protein_affected, self.protein_effect_location))



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
        Returns:
            integer: The change between the ref and alt alleles due to the variant presence. SNVs lead to change length of zero,
            deletions and insertions/duplications lead to negative and positive change lengths, respectively.
        """
        return self._change_length

    @property
    def variant_key(self) -> str:
        """
        Get a readable representation of the variant's coordinates.

        For instance, `X_12345_12345_C_G` for sequence variant or `22_10001_20000_INV`
        Note that both start and end coordinates use 1-based (included) coordinate system.
        """
        if self.is_structural():
            return f'{self.chrom}_{self.start + 1}_{self.end}_{self.alt[1:-1]}'
        else:
            return f'{self.chrom}_{self.start + 1}_{self.end}_{self.ref}_{self.alt}'

    @property
    def variant_class(self) -> str:
        """
        Returns:
            string: The variant class. (e.g. `DUP`, `SNV`, `INS`, `MNV`, `INV`, ...)
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
        """Checks if the variant coordinates use structural variant notation
        (e.g. `chr5  101 . N <DEL> .  .  SVTYPE=DEL;END=120;SVLEN=-10`)
        as opposed to the sequence/literal notation (`chr5  101 . NACGTACGTAC N`).

        Returns:
            boolean: True if the variant coordinates use structural variant notation
        """
        return len(self._alt) != 0 and self._alt.startswith('<') and self._alt.endswith('>')

    def __len__(self):
        """
        Get the number of bases on the ref allele that are affected by the variant.
        """
        return len(self._region)

    def __eq__(self, other) -> bool:
        return isinstance(other, VariantCoordinates) \
            and self.alt == other.alt \
            and self.ref == other.ref \
            and self.chrom == other.chrom \
            and self.start == other.start \
            and self.end == other.end \
            and self.change_length == other.change_length

    def __str__(self) -> str:
        return f"VariantCoordinates(chrom={self.chrom}, " \
               f"start={self.start}, end={self.end}, " \
               f"ref={self.ref}, alt={self.alt}, " \
               f"change_length={self.change_length})"

    def __repr__(self) -> str:
        return str(self)

    def __hash__(self) -> int:
        return hash((self._region, self._ref, self._alt, self._change_length))


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
                                    hgvsc_id: str,
                                    is_preferred: bool,
                                    consequences: typing.Sequence[str],
                                    exons_effected: typing.Sequence[int],
                                    protein: typing.Sequence[ProteinMetadata],
                                    protein_effect_start: int,
                                    protein_effect_end: int,
                                    genotypes: Genotypes):
        transcript = TranscriptAnnotation(gene_name, trans_id, hgvsc_id, is_preferred, consequences, exons_effected, protein,
                                          protein_effect_start, protein_effect_end)
        return Variant(variant_coordinates, [transcript], genotypes)

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
