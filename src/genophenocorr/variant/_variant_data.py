import typing
from genophenocorr.protein import ProteinMetadata


class VariantCoordinates:
    """A representation of coordinates of sequence and symbolic variants.
    The breakend variants are not supported.

    Attributes:
        chrom (string): Chromosome variant coordinates are located on
        start (integer): The 0-based starting coordinate of the ref allele
        end (integer): The 0-based ending coordinate of the ref allele
        ref (string): The reference allele
        alt (string): The alternate allele
        change_length (integer): The change between the ref and alt alleles due to the variant presence
        genotype (string): The genotype of the variant (e.g. Heterozygous, Homozygous)
    Methods:
        is_structural: Returns True if the variant coordinates use structural variant notation
        as_string: Returns a readable representation of the variant coordinate
    """

    def __init__(self, chrom: str, start: int, end: int, ref: str, alt: str, change_length: int, genotype: str):
        """Constructs all necessary attributes for a VariantCoordinates object

        Args:
            chrom (string): Chromosome variant coordinates are located on
            start (integer): The 0-based starting coordinate of the ref allele
            end (integer): The 0-based ending coordinate of the ref allele
            ref (string): The reference allele
            alt (string): The alternate allele
            change_length (integer): The change between the ref and alt alleles due to the variant presence
            genotype (string): The genotype of the variant (e.g. Heterozygous, Homozygous)
        """
        # TODO(lnrekerle) - instance/type check
        # TODO - id?
        self._chrom = chrom
        self._start = start
        self._end = end
        self._ref = ref
        self._alt = alt
        self._change_length = change_length
        self._genotype = genotype

    @property
    def chrom(self) -> str:
        """
        Returns:
            string: The label of the chromosome/contig where the variant is located.
        """
        return self._chrom

    @property
    def start(self) -> int:
        """
        Returns:
            integer: The 0-based start coordinate (excluded) of the ref allele.
        """
        return self._start

    @property
    def end(self) -> int:
        """
        Returns:
            integer: The 0-based end coordinate (included) of the ref allele.
        """
        return self._end

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
    def genotype(self) -> str:
        # TODO - add a doc string with an example. Do we return `0/1`, or `HET`, or GENO_0000135
        #  (http://purl.obolibrary.org/obo/GENO_0000135)? All these are strings.
        return self._genotype

    def is_structural(self) -> bool:
        """Checks if the variant coordinates use structural variant notation
        (e.g. `chr5  101 . N <DEL> .  .  SVTYPE=DEL;END=120;SVLEN=-10`)
        as opposed to the sequence/literal notation (`chr5  101 . NACGTACGTAC N`).
        
        Returns:
            boolean: True if the variant coordinates use structural variant notation
        """
        return len(self._alt) != 0 and self._alt.startswith('<') and self._alt.endswith('>')

    def as_string(self) -> str:
        """
        Returns:
            string: A readable representation of the variant coordinates
        """
        return f"{self.chrom}_{self.start}_{self.end}_{self.ref}_{self.alt}_{self.genotype}"

    def __len__(self):
        """
        Get the number of bases on the ref allele that are affected by the variant.
        """
        return self._end - self._start

    def __eq__(self, other) -> bool:
        return isinstance(other, VariantCoordinates) \
            and self.alt == other.alt \
            and self.ref == other.ref \
            and self.chrom == other.chrom \
            and self.start == other.start \
            and self.end == other.end \
            and self.change_length == other.change_length \
            and self.genotype == other.genotype

    def __str__(self) -> str:
        return f"VariantCoordinates(chrom={self.chrom}, " \
            f"start={self.start}, end={self.end}, " \
            f"ref={self.ref}, alt={self.alt}, " \
            f"change_length={self.change_length}, " \
            f"genotype={self.genotype})"

    def __repr__(self) -> str:
        return str(self)

    def __hash__(self) -> int:
        return hash((self._chrom, self._start, self._end, self._ref, self._alt, self._change_length, self._genotype))


class TranscriptAnnotation:
    """Class that represents results of the functional annotation of a variant with respect to single transcript of a gene.

    Attributes:
        gene_id (string): The gene symbol associated with the transcript
        transcript_id (string): The transcript ID
        hgvsc_id (string): The HGVS "coding-DNA" ID if available, else None
        variant_effects (Sequence[string]): A sequence of predicted effects given by VEP
        overlapping_exons (Sequence[integer]): A sequence of exons affected by the variant. Returns None if none are affected.
        protein_affected (ProteinMetadata): A ProteinMetadata object representing the protein affected by this transcript
        protein_effect_location (Tuple(integer, integer)): The start and end coordinates of the effect on the protein sequence.
    """
    def __init__(self, gene_id: str,
                 tx_id: str,
                 hgvsc: typing.Optional[str],
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
    def hgvsc_id(self):
        """
        Returns:
            string: The HGVS "coding-DNA" representation of the variant (e.g. NM_123456.7:c.9876G>T)
        """
        return self._hgvsc_id

    @property
    def variant_effects(self):
        """
        Returns:
            Sequence[string]: A sequence of variant effects. 
                    Definitions of these can be found at: http://www.sequenceontology.org/
        """
        return self._variant_effects

    @property
    def overlapping_exons(self):
        """
        Returns:
            Sequence[integer]: A sequence of IDs of the exons that overlap with the variant.
        """
        return self._affected_exons

    @property
    def protein_affected(self) -> ProteinMetadata:
        """
        Returns:
            ProteinMetadata: The ProteinMetadata object representing the protein that is affected by the alteration of this transcript
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
            f"variant_effects:{self.variant_effects}," \
            f"overlapping_exons:{self.overlapping_exons}," \
            f"protein_affected:{self.protein_affected}," \
            f"protein_effect_location:{self.protein_effect_location})"

    def __eq__(self, other) -> bool:
        return isinstance(other, TranscriptAnnotation) \
            and self.gene_id == other.gene_id \
            and self.hgvsc_id == other.hgvsc_id \
            and self.transcript_id == other.transcript_id \
            and self.variant_effects == other.variant_effects \
            and self.overlapping_exons == other.overlapping_exons \
            and self.protein_affected == other.protein_affected \
            and self.protein_effect_location == other.protein_effect_location

    def __repr__(self) -> str:
        return str(self)

    def __hash__(self) -> int:
        return hash((self.gene_id, self.hgvsc_id, self.transcript_id, self.overlapping_exons, self.variant_effects, self.protein_affected, self.protein_effect_location))


class Variant:
    """Class that represents results of the functional annotation of a variant with all included transcripts.

    Attributes:
        variant_coordinates (VariantCoordinates): A VariantCoordinates object with coordinates for this Variant
        variant_string (string): A readable representation of the variant coordinates
        genotype (string, Optional): The genotype of the variant (e.g. Homozygous, Heterozygous)
        tx_annotations (Sequence[TranscriptAnnotation], Optional): A sequence of TranscriptAnnotation objects representing transcripts affected by this variant
        variant_class (string): The variant class (e.g. Duplication, SNV, etc.)
    """

    def __init__(self, var_id: str,
                 var_class: str,
                 var_coordinates: VariantCoordinates,
                 tx_annotations: typing.Optional[typing.Sequence[TranscriptAnnotation]],
                 genotype: typing.Optional[str]):
        """Constructs all necessary attributes for a Variant object

        Args:
            var_coordinates (VariantCoordinates): A VariantCoordinates object with coordinates for this Variant
            var_id (string): A readable representation of the variant coordinates
            genotype (string, Optional): The genotype of the variant (e.g. Homozygous, Heterozygous)
            tx_annotations (Sequence[TranscriptAnnotation], Optional): A sequence of TranscriptAnnotation objects representing transcripts affected by this variant
            var_class (string): The variant class (e.g. Duplication, SNV, etc.)
        """
        self._id = var_id
        self._var_coordinates = var_coordinates
        self._var_class = var_class
        if tx_annotations is None:
            self._tx_annotations = None
        else:
            self._tx_annotations = tuple(tx_annotations)
        self._genotype = genotype

    @property
    def variant_coordinates(self) -> VariantCoordinates:
        """
        Returns:
            VariantCoordinates: A representation of coordinates of a sequence and symbolic variant.
        """
        return self._var_coordinates

    @property
    def variant_string(self) -> str:
        """
        Returns:
            string: A readable representation of the variant's coordinates. 
                Format - "Chromosome_Start_Reference/Alternative" or
                "Chromosome_Start_StructuralType"
        """
        return self._id

    @property
    def genotype(self) -> str:
        """Optional parameter. Required for recessive tests. 
        Possible values: Heterozygous, Homozygous, Hemizygous
        
        Returns:
            string: Genotype of the variant
        """
        return self._genotype

    @property
    def tx_annotations(self) -> typing.Sequence[TranscriptAnnotation]:
        """A collection of TranscriptAnnotations that each represent results of the functional annotation 
        of a variant with respect to single transcript of a gene.

        Returns:
            Sequence[TranscriptAnnotation]: A sequence of TranscriptAnnotation objects
        """
        return self._tx_annotations

    @property
    def variant_class(self) -> str:
        """
        Returns:
            string: The variant class. (e.g. Duplication, SNV, Deletion, etc.)
        """
        return self._var_class

    def __eq__(self, other) -> bool:
        return isinstance(other, Variant) \
            and self.variant_string == other.variant_string \
            and self.genotype == other.genotype \
            and self.variant_class == other.variant_class \
            and self.variant_coordinates == other.variant_coordinates \
            and self.tx_annotations == other.tx_annotations

    def __hash__(self) -> int:
        return hash((self.variant_coordinates, self.variant_string, self.variant_class, self.genotype, self.tx_annotations))

    def __repr__(self) -> str:
        return str(self)

    def __str__(self) -> str:
        return f"Variant(variant_coordinates:{str(self.variant_coordinates)}," \
            f"variant_string:{self.variant_string}," \
            f"genotype:{self.genotype}," \
            f"tx_annotations:{self.tx_annotations}," \
            f"variant_class:{self.variant_class})"
