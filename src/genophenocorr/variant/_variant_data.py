import typing


class VariantCoordinates:
    """A representation of coordinates of sequence and symbolic variants.

    The breakend variants are not supported.
    """

    def __init__(self, chrom: str, start: int, end: int, ref: str, alt: str, change_length: int, genotype: str):
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
        Get label of the chromosome/contig where the variant is located.
        """
        return self._chrom

    @property
    def start(self) -> int:
        """
        Get 0-based start coordinate (excluded) of the ref allele.
        """
        return self._start

    @property
    def end(self) -> int:
        """
        Get 0-based end coordinate (included) of the ref allele.
        """
        return self._end

    @property
    def ref(self) -> str:
        """
        Get reference allele (e.g. "A", "N"). The allele may be an empty string.
        """
        return self._ref

    @property
    def alt(self) -> str:
        """
        Get alternate allele (e.g. "A", "GG", "<DEL>"). The allele may be an empty string for sequence variants.
        The symbolic alternate allele follow the VCF notation and use the `<` and `>` characters
        (e.g. "<DEL>", "<INS:ME:SINE>").
        """
        return self._alt

    @property
    def change_length(self) -> int:
        """
        Get the change between the ref and alt alleles due to the variant presence. SNVs lead to change length of zero,
        deletions and insertions/duplications lead to negative and positive change lengths, respectively.
        """
        return self._change_length

    @property
    def genotype(self) -> str:
        # TODO - add a doc string with an example. Do we return `0/1`, or `HET`, or GENO_0000135
        #  (http://purl.obolibrary.org/obo/GENO_0000135)? All these are strings.
        return self._genotype

    def is_structural(self) -> bool:
        """
        Return `True` if the variant coordinates use structural variant notation
        (e.g. `chr5  101 . N <DEL> .  .  SVTYPE=DEL;END=120;SVLEN=-10`)
        as opposed to the sequence/literal notation (`chr5  101 . NACGTACGTAC N`).
        """
        return len(self._alt) != 0 and self._alt.startswith('<') and self._alt.endswith('>')

    def as_string(self) -> str:
        return f"{self.chrom}_{self.start}_{self.end}_{self.ref}_{self.alt}"

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
    """
    Class that represents results of the functional annotation of a variant with respect to single transcript of a gene.
    """

    def __init__(self, gene_id: str,
                 tx_id: str,
                 hgvsc: typing.Optional[str],
                 variant_effects,
                 affected_exons: typing.Optional[typing.Sequence[int]],
                 affected_protein: str,
                 protein_effect_start: typing.Optional[int],
                 protein_effect_end: typing.Optional[int]):
        self._gene_id = gene_id
        self._tx_id = tx_id
        self._hgvsc_id = hgvsc
        self._variant_effects = tuple(variant_effects)
        if affected_exons is not None:
            self._affected_exons = tuple(affected_exons)
        else:
            self._affected_exons = None
        self._affected_protein = affected_protein
        self._protein_effect_location = (protein_effect_start, protein_effect_end)

    @property
    def gene_id(self) -> str:
        """
        Get the gene symbol (e.g. SURF1)
        """
        return self._gene_id

    @property
    def transcript_id(self) -> str:
        """
        Get the transcript identifier (e.g. NM_123456.7)
        """
        return self._tx_id

    @property
    def hgvsc_id(self):
        """
        Get the HGVS "coding-DNA" representation of the variant (e.g. NM_123456.7:c.9876G>T)
        """
        return self._hgvsc_id

    @property
    def variant_effects(self):
        """
        Get a sequence of variant effects. 
        Definitions of these can be found at: http://www.sequenceontology.org/
        """
        return self._variant_effects

    @property
    def overlapping_exons(self):
        """
        Get a sequence of IDs of the exons that overlap with the variant.
        """
        return self._affected_exons

    @property
    def protein_affected(self) -> str:
        """
        Get the protein ID that is affected by the alteration of this transcript (e.g. NP_037407.4)
        """
        return self._affected_protein

    @property
    def protein_effect_location(self) -> typing.Tuple[int, int]:
        """
        Get the start and end position on the protein sequence that the variant effects. (e.g. [1234, 1235])
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
    """
    Class that represents results of the functional annotation of a variant with all included transcripts.
    """

    def __init__(self, var_id: str,
                 var_class: str,
                 var_coordinates: VariantCoordinates,
                 current_tx: str,
                 tx_annotations: typing.Optional[typing.Sequence[TranscriptAnnotation]],
                 genotype: typing.Optional[str]):
        # TODO - revert
        self._id = var_id
        self._var_coordinates = var_coordinates
        self._var_class = var_class
        self._current_tx = current_tx
        if tx_annotations is not None:
            self._tx_annotations = tuple(tx_annotations)
        else:
            self._tx_annotations = None
        self._genotype = genotype

    @property
    def variant_coordinates(self) -> VariantCoordinates:
        """
        A representation of coordinates of a sequence and symbolic variant.
        """
        return self._var_coordinates

    @property
    def variant_string(self) -> str:
        """
        A readable representation of the variant's coordinates. 
        Format - "Chromosome_Start_Reference/Alternative" or
                 "Chromosome_Start_StructuralType"
        """
        return self._id

    @property
    def genotype(self) -> str:
        """
        Optional parameter. Required for recessive tests. 
        Possible values: Heterozygous, Homozygous, Hemizygous
        """
        return self._genotype

    @property
    def selected_transcript(self) -> str:
        """
        The ID of either the canonical transcript given by VEP or the transcript chosen by the user.
        """
        return self._current_tx

    @property
    def tx_annotations(self) -> typing.Sequence[TranscriptAnnotation]:
        """
        A collection of TranscriptAnnotations that each represent results of the functional annotation 
        of a variant with respect to single transcript of a gene.
        """
        return self._tx_annotations

    @property
    def variant_class(self) -> str:
        """
        The variant class. (e.g. Duplication, SNV, Deletion, etc.)
        """
        return self._var_class

    def __eq__(self, other) -> bool:
        # TODO - this seems a bit odd. Investigate.
        return isinstance(other, Variant) \
            and self.variant_string == other.variant_string \
            and self.genotype == other.genotype \
            and self.variant_class == other.variant_class \
            and self.variant_coordinates == other.variant_coordinates 
            #and self.tx_annotations == other.tx_annotations

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
