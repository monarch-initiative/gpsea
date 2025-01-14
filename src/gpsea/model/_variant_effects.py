import enum
import typing

import hpotk


class VariantEffect(enum.Enum):
    """
    `VariantEffect` represents consequences of a variant on transcript that are supported by GPSEA.

    .. doctest::

      >>> from gpsea.model import VariantEffect
      >>> missense = VariantEffect.MISSENSE_VARIANT
      >>> print(missense)
      missense_variant

    The `VariantEffect` has a :attr:`curie` attribute that represents the ontology class from
    `Sequence Ontology <http://www.sequenceontology.org/>`_.

    .. doctest::

      >>> missense.curie
      'SO:0001583'
    """

    TRANSCRIPT_AMPLIFICATION = "SO:0001889"
    TRANSCRIPT_ABLATION = "SO:0001893"
    TRANSCRIPT_TRANSLOCATION = "SO:0001883"
    SPLICE_ACCEPTOR_VARIANT = "SO:0001574"
    SPLICE_DONOR_VARIANT = "SO:0001575"
    STOP_GAINED = "SO:0001587"
    FRAMESHIFT_VARIANT = "SO:0001589"
    STOP_LOST = "SO:0001578"
    START_LOST = "SO:0002012"
    INFRAME_INSERTION = "SO:0001821"
    INFRAME_DELETION = "SO:0001822"
    MISSENSE_VARIANT = "SO:0001583"
    PROTEIN_ALTERING_VARIANT = "SO:0001818"
    SPLICE_REGION_VARIANT = "SO:0001630"
    SPLICE_DONOR_5TH_BASE_VARIANT = "SO:0001787"
    SPLICE_DONOR_REGION_VARIANT = "SO:0002170"
    SPLICE_POLYPYRIMIDINE_TRACT_VARIANT = "SO:0002169"
    INCOMPLETE_TERMINAL_CODON_VARIANT = "SO:0001626"
    START_RETAINED_VARIANT = "SO:0002019"
    STOP_RETAINED_VARIANT = "SO:0001567"
    SYNONYMOUS_VARIANT = "SO:0001819"
    CODING_SEQUENCE_VARIANT = "SO:0001580"
    MATURE_MIRNA_VARIANT = "SO:0001620"
    FIVE_PRIME_UTR_VARIANT = "SO:0001623"
    THREE_PRIME_UTR_VARIANT = "SO:0001624"
    NON_CODING_TRANSCRIPT_EXON_VARIANT = "SO:0001792"
    INTRON_VARIANT = "SO:0001627"
    NMD_TRANSCRIPT_VARIANT = "SO:0001621"
    NON_CODING_TRANSCRIPT_VARIANT = "SO:0001619"
    UPSTREAM_GENE_VARIANT = "SO:0001631"
    DOWNSTREAM_GENE_VARIANT = "SO:0001632"
    TFBS_ABLATION = "SO:0001895"
    TFBS_AMPLIFICATION = "SO:0001892"
    TF_BINDING_SITE_VARIANT = "SO:0001782"
    REGULATORY_REGION_ABLATION = "SO:0001894"
    REGULATORY_REGION_AMPLIFICATION = "SO:0001891"
    FEATURE_ELONGATION = "SO:0001907"
    REGULATORY_REGION_VARIANT = "SO:0001566"
    FEATURE_TRUNCATION = "SO:0001906"
    INTERGENIC_VARIANT = "SO:0001628"
    SEQUENCE_VARIANT = "SO:0001060"

    def to_display(self) -> str:
        """
        Get a concise name of the variant effect that is suitable for showing to humans.

        Example
        ^^^^^^^

        >>> from gpsea.model import VariantEffect
        >>> VariantEffect.MISSENSE_VARIANT.to_display()
        'missense'
        >>> VariantEffect.SPLICE_DONOR_5TH_BASE_VARIANT.to_display()
        'splice donor 5th base'

        :returns: a `str` with the name or `'n/a'` if the variant effect was not assigned a concise name.
        """
        return effect_to_display.get(self, "n/a")

    @staticmethod
    def structural_so_id_to_display(so_term: typing.Union[hpotk.TermId, str]) -> str:
        """
        Get a `str` with a concise name for representing a Sequence Ontology (SO) term identifier.

        Example
        ^^^^^^^

        >>> from gpsea.model import VariantEffect
        >>> VariantEffect.structural_so_id_to_display('SO:1000029')
        'chromosomal deletion'

        :param so_term: a CURIE `str` or a :class:`~hpotk.TermId` with the query SO term.
        :returns: a `str` with the concise name for the SO term or `'n/a'` if a name has not been assigned yet.
        """
        if isinstance(so_term, hpotk.TermId):
            so_term = so_term.value

        return so_id_to_display.get(so_term, "n/a")

    def __init__(self, curie: str):
        self._curie = curie

    @property
    def curie(self) -> str:
        """
        Get a compact URI (CURIE) of the variant effect
        (e.g. `SO:0001583` for a missense variant).
        """
        return self._curie

    def __str__(self) -> str:
        return self.name.lower()


effect_to_display = {
    VariantEffect.TRANSCRIPT_AMPLIFICATION: "transcript amplification",
    VariantEffect.TRANSCRIPT_ABLATION: "transcript ablation",
    VariantEffect.TRANSCRIPT_TRANSLOCATION: "transcript translocation",
    VariantEffect.SPLICE_ACCEPTOR_VARIANT: "splice acceptor",
    VariantEffect.SPLICE_DONOR_VARIANT: "splice donor",
    VariantEffect.STOP_GAINED: "stop gained",
    VariantEffect.FRAMESHIFT_VARIANT: "frameshift",
    VariantEffect.STOP_LOST: "stop lost",
    VariantEffect.START_LOST: "start lost",
    VariantEffect.INFRAME_INSERTION: "inframe insertion",
    VariantEffect.INFRAME_DELETION: "inframe deletion",
    VariantEffect.MISSENSE_VARIANT: "missense",
    VariantEffect.PROTEIN_ALTERING_VARIANT: "protein altering",
    VariantEffect.SPLICE_REGION_VARIANT: "splice region",
    VariantEffect.SPLICE_DONOR_5TH_BASE_VARIANT: "splice donor 5th base",
    VariantEffect.SPLICE_DONOR_REGION_VARIANT: "splice donor",
    VariantEffect.SPLICE_POLYPYRIMIDINE_TRACT_VARIANT: "splice polypyrimidine",
    VariantEffect.INCOMPLETE_TERMINAL_CODON_VARIANT: "incomplete terminal codon",
    VariantEffect.START_RETAINED_VARIANT: "start retained",
    VariantEffect.STOP_RETAINED_VARIANT: "stop retainined",
    VariantEffect.SYNONYMOUS_VARIANT: "synonymous",
    VariantEffect.CODING_SEQUENCE_VARIANT: "coding sequence",
    VariantEffect.MATURE_MIRNA_VARIANT: "mature miRNA",
    VariantEffect.FIVE_PRIME_UTR_VARIANT: "5UTR",
    VariantEffect.THREE_PRIME_UTR_VARIANT: "3UTR",
    VariantEffect.NON_CODING_TRANSCRIPT_EXON_VARIANT: "non-coding transcript exon",
    VariantEffect.INTRON_VARIANT: "intronic",
    VariantEffect.NMD_TRANSCRIPT_VARIANT: "NMD transcript",
    VariantEffect.NON_CODING_TRANSCRIPT_VARIANT: "non-coding transcript",
    VariantEffect.UPSTREAM_GENE_VARIANT: "upstream of gene",
    VariantEffect.DOWNSTREAM_GENE_VARIANT: "downstream of gene",
    VariantEffect.TFBS_ABLATION: "TFBS ablation",
    VariantEffect.TFBS_AMPLIFICATION: "TFBS amplification",
    VariantEffect.TF_BINDING_SITE_VARIANT: "TFBS binding site",
    VariantEffect.REGULATORY_REGION_ABLATION: "regulatory region ablation",
    VariantEffect.REGULATORY_REGION_AMPLIFICATION: "regulatory region amplification",
    VariantEffect.FEATURE_ELONGATION: "feature elongation",
    VariantEffect.REGULATORY_REGION_VARIANT: "regulatory region",
    VariantEffect.FEATURE_TRUNCATION: "feature truncation",
    VariantEffect.INTERGENIC_VARIANT: "intergenic",
    VariantEffect.SEQUENCE_VARIANT: "sequence variant",
}

so_id_to_display = {
    "SO:1000029": "chromosomal deletion",
    "SO:1000037": "chromosomal duplication",
    "SO:1000030": "chromosomal_inversion",
    "SO:1000044": "chromosomal_translocation",
}
