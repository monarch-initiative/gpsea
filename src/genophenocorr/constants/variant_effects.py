import abc
import os
import json


class VariantEffect():

    def __init__(self, term, accession) -> None:
        self._id = term
        self._so_id = accession

    @property
    def effect_id(self):
        return self._id

    @property
    def effect_name(self):
        return self.effect_id.replace('_', ' ')

    @property
    def SO_id(self):
        return self._so_id

TRANSCRIPT_ABLATION = VariantEffect('transcript_ablation', "SO:0001893")
SPLICE_ACCEPTOR_VARIANT = VariantEffect('splice_acceptor_variant', "SO:0001574")
SPLICE_DONOR_VARIANT =  VariantEffect('splice_donor_variant', "SO:0001575")
STOP_GAINED =  VariantEffect('stop_gained', "SO:0001587")
FRAMESHIFT_VARIANT =  VariantEffect('frameshift_variant', "SO:0001589")
STOP_LOST =  VariantEffect('stop_lost', "SO:0001578")
START_LOST =  VariantEffect('start_lost', "SO:0002012")
TRANSCRIPT_AMPLIFICATION =  VariantEffect('transcript_amplification', "SO:0001889")
INFRAME_INSERTION =  VariantEffect('inframe_insertion',  "SO:0001821")
INFRAME_DELETION =  VariantEffect('inframe_deletion', "SO:0001822")
MISSENSE_VARIANT =  VariantEffect('missense_variant', "SO:0001583")
PROTEIN_ALTERING_VARIANT =  VariantEffect('protein_altering_variant', "SO:0001818")
SPLICE_REGION_VARIANT =  VariantEffect('splice_region_variant', "SO:0001630")
SPLICE_DONOR_5TH_BASE_VARIANT =  VariantEffect('splice_donor_5th_base_variant', "SO:0001787")
SPLICE_DONOR_REGION_VARIANT =  VariantEffect('splice_donor_region_variant', "SO:0002170")
SPLICE_POLYPYRIMIDINE_TRACT_VARIANT =  VariantEffect('splice_polypyrimidine_tract_variant', "SO:0002169")
INCOMPLETE_TERMINAL_CODON_VARIANT =  VariantEffect('incomplete_terminal_codon_variant', "SO:0001626")
START_RETAINED_VARIANT =  VariantEffect('start_retained_variant', "SO:0002019")
STOP_RETAINED_VARIANT =  VariantEffect('stop_retained_variant', "SO:0001567")
SYNONYMOUS_VARIANT =  VariantEffect('synonymous_variant', "SO:0001819")
CODING_SEQUENCE_VARIANT =  VariantEffect('coding_sequence_variant', "SO:0001580")
MATURE_MIRNA_VARIANT =  VariantEffect('mature_miRNA_variant', "SO:0001620")
FIVE_PRIME_UTR_VARIANT =  VariantEffect('5_prime_UTR_variant', "SO:0001623")
THREE_PRIME_UTR_VARIANT =  VariantEffect('3_prime_UTR_variant', "SO:0001624")
NON_CODING_TRANSCRIPT_EXON_VARIANT =  VariantEffect('non_coding_transcript_exon_variant', "SO:0001792")
INTRON_VARIANT =  VariantEffect('intron_variant', "SO:0001627")
NMD_TRANSCRIPT_VARIANT =  VariantEffect('NMD_transcript_variant', "SO:0001621",)
NON_CODING_TRANSCRIPT_VARIANT =  VariantEffect('non_coding_transcript_variant', "SO:0001619")
UPSTREAM_GENE_VARIANT =  VariantEffect('upstream_gene_variant', "SO:0001631")
DOWNSTREAM_GENE_VARIANT =  VariantEffect('downstream_gene_variant', "SO:0001632")
TFBS_ABLATION =  VariantEffect('TFBS_ablation', "SO:0001895")
TFBS_AMPLIFICATION =  VariantEffect('TFBS_amplification', "SO:0001892")
TF_BINDING_SITE_VARIANT =  VariantEffect('TF_binding_site_variant', "SO:0001782")
REGULATORY_REGION_ABLATION =  VariantEffect('regulatory_region_ablation', "SO:0001894")
REGULATORY_REGION_AMPLIFICATION =  VariantEffect('regulatory_region_amplification', "SO:0001891")
FEATURE_ELONGATION =  VariantEffect('feature_elongation', "SO:0001907")
REGULATORY_REGION_VARIANT =  VariantEffect('regulatory_region_variant', "SO:0001566")
FEATURE_TRUNCATION =  VariantEffect('feature_truncation', "SO:0001906")
INTERGENIC_VARIANT =  VariantEffect('intergenic_variant',  "SO:0001628")
SEQUENCE_VARIANT = VariantEffect('sequence_variant', "SO:0001060")