import abc
import os
import json

with open(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'variant_effect.json'), 'r') as json_file:
    VE_LIST = json.load(json_file)


class VariantEffect():

    def __init__(self, effect_info: dict) -> None:
        self._id = effect_info.get('SO_term')
        self._descrip = effect_info.get('description')
        self._so_id = effect_info.get('SO_accession')
        self._name = effect_info.get('label')

    @staticmethod
    def from_json(name:str):
        if name is None:
            raise ValueError('Valid Variant Effect must be given.')
        for effect in VE_LIST:
            if effect.get('SO_term') == name:
                return VariantEffect(effect)

    @property
    def effect_id(self):
        return self._id

    @property
    def effect_name(self):
        return self._name

    @property
    def effect_description(self):
        return self._descrip

    @property
    def SO_id(self):
        return self._so_id

TRANSCRIPT_ABLATION = VariantEffect.from_json('transcript_ablation')
SPLICE_ACCEPTOR_VARIANT = VariantEffect.from_json('splice_acceptor_variant')
SPLICE_DONOR_VARIANT =  VariantEffect.from_json('splice_donor_variant' )
STOP_GAINED =  VariantEffect.from_json('stop_gained' )
FRAMESHIFT_VARIANT =  VariantEffect.from_json('frameshift_variant' )
STOP_LOST =  VariantEffect.from_json('stop_lost' )
START_LOST =  VariantEffect.from_json('start_lost' )
TRANSCRIPT_AMPLIFICATION =  VariantEffect.from_json('transcript_amplification' )
INFRAME_INSERTION =  VariantEffect.from_json('inframe_insertion' )
INFRAME_DELETION =  VariantEffect.from_json('inframe_deletion' )
MISSENSE_VARIANT =  VariantEffect.from_json('missense_variant' )
PROTEIN_ALTERING_VARIANT =  VariantEffect.from_json('protein_altering_variant' )
SPLICE_REGION_VARIANT =  VariantEffect.from_json('splice_region_variant' )
SPLICE_DONOR_5TH_BASE_VARIANT =  VariantEffect.from_json('splice_donor_5th_base_variant' )
SPLICE_DONOR_REGION_VARIANT =  VariantEffect.from_json('splice_donor_region_variant' )
SPLICE_POLYPYRIMIDINE_TRACT_VARIANT =  VariantEffect.from_json('splice_polypyrimidine_tract_variant' )
INCOMPLETE_TERMINAL_CODON_VARIANT =  VariantEffect.from_json('incomplete_terminal_codon_variant' )
START_RETAINED_VARIANT =  VariantEffect.from_json('start_retained_variant' )
STOP_RETAINED_VARIANT =  VariantEffect.from_json('stop_retained_variant' )
SYNONYMOUS_VARIANT =  VariantEffect.from_json('synonymous_variant' )
CODING_SEQUENCE_VARIANT =  VariantEffect.from_json('coding_sequence_variant' )
MATURE_MIRNA_VARIANT =  VariantEffect.from_json('mature_miRNA_variant' )
FIVE_PRIME_UTR_VARIANT =  VariantEffect.from_json('5_prime_UTR_variant' )
THREE_PRIME_UTR_VARIANT =  VariantEffect.from_json('3_prime_UTR_variant' )
NON_CODING_TRANSCRIPT_EXON_VARIANT =  VariantEffect.from_json('non_coding_transcript_exon_variant' )
INTRON_VARIANT =  VariantEffect.from_json('intron_variant' )
NMD_TRANSCRIPT_VARIANT =  VariantEffect.from_json('NMD_transcript_variant' )
NON_CODING_TRANSCRIPT_VARIANT =  VariantEffect.from_json('non_coding_transcript_variant' )
UPSTREAM_GENE_VARIANT =  VariantEffect.from_json('upstream_gene_variant' )
DOWNSTREAM_GENE_VARIANT =  VariantEffect.from_json('downstream_gene_variant' )
TFBS_ABLATION =  VariantEffect.from_json('TFBS_ablation' )
TFBS_AMPLIFICATION =  VariantEffect.from_json('TFBS_amplification' )
TF_BINDING_SITE_VARIANT =  VariantEffect.from_json('TF_binding_site_variant' )
REGULATORY_REGION_ABLATION =  VariantEffect.from_json('regulatory_region_ablation' )
REGULATORY_REGION_AMPLIFICATION =  VariantEffect.from_json('regulatory_region_amplification' )
FEATURE_ELONGATION =  VariantEffect.from_json('feature_elongation' )
REGULATORY_REGION_VARIANT =  VariantEffect.from_json('regulatory_region_variant' )
FEATURE_TRUNCATION =  VariantEffect.from_json('feature_truncation' )
INTERGENIC_VARIANT =  VariantEffect.from_json('intergenic_variant' )