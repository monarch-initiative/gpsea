from curses.ascii import isdigit
import varcode as vc
import pyensembl
import re

class Variant:
    def __init__(self,ref, genoInterp, transcript):
        self._genoInterp = genoInterp
        self._variant = self.__find_variant(reference = ref)
        self._transcript = transcript
        
    def __find_variant(self, reference):
        varInterp = self._genoInterp.variant_interpretation.variation_descriptor
        self._allelic_state = varInterp.allelic_state.label
        contig = re.sub(r'[^0-9MXY]', '', varInterp.vcf_record.chrom)
        if len(contig) == 0 or (contig.isdigit() and (int(contig) == 0 or int(contig) >= 24)):
            ## Chromosome can only be values 1-23, X, Y, or M
            raise ValueError(f"Contig did not work: {varInterp.vcf_record.chrom}")
        start = varInterp.vcf_record.pos
        ref = varInterp.vcf_record.ref
        alt = varInterp.vcf_record.alt
        if reference.lower() == 'hg37' or reference.lower() == 'grch37' or reference.lower() == 'hg19':
            ens = pyensembl.ensembl_grch37
        elif reference.lower() == 'hg38' or reference.lower() == 'grch38':
            ens = pyensembl.ensembl_grch38
        else:
            raise ValueError(f'Unknown Reference {reference}. Please use hg19 or hg38.')
        myVar = vc.Variant(contig, start, ref, alt, ensembl = ens)
        return myVar


    @property
    def variant(self):
        return self._variant

    @property
    def variant_string(self):
        return self.variant.short_description

    @property
    def variant_types(self):
        all_types = []
        if not self.top_effect.short_description.endswith("*") and self.variant.is_snv:
            all_types.append('missense')
        if self.top_effect.short_description.endswith("*") and self.variant.is_snv:
            all_types.append('nonsense')
        if 'dup' in self.top_effect.short_description:
            all_types.append('duplication')
        if self.variant.is_deletion:
            all_types.append('deletion')
        if self.variant.is_insertion:
            all_types.append('insertion')
        if self.variant.is_transition:
            all_types.append('transition')
        if self.variant.is_transversion:
            all_types.append('transversion')
        if self.variant.is_indel:
            all_types.append('indel')
        return all_types

    @property
    def genotype(self):
        if self._allelic_state is not None:
            return self._allelic_state
        else:
            return None

    @property
    def genomic_location(self):
        if self.variant is not None:
            return self.variant.start
        else:
            return None

    @property
    def all_effects(self):
        if self.variant is not None:
            return self.variant.effects().short_string()
        else:
            return None

    @property
    def top_effect(self):
        if self.all_effects is not None and self._transcript is None:
            return self.variant.effects().top_priority_effect()
        elif self.all_effects is not None and self._transcript is not None:
            return self.variant.effects().top_priority_effect_per_transcript_id().get(self._transcript)
        else:
            return None

    @property
    def all_transcripts(self):
        if self._variant is not None:
            return self.variant.transcripts
        else:
            return None

    @property
    def all_transcript_ids(self):
        if self.variant.transcripts is not None:
            return self.variant.transcript_ids
        else:
            return None

    @property
    def top_effect_transcript(self):
        if self.top_effect is not None and self._transcript is None:
            return self.top_effect.transcript
        elif self.top_effect is not None and self._transcript is not None:
            return [trans for trans in self.variant.transcripts if trans.id == self._transcript][0]
        else:
            return None

    @property
    def top_effected_protein(self):
        if self.top_effect_transcript is not None:
            return self.top_effect_transcript.protein_id
        else:
            return None

    @property
    def protein_effect_location(self):
        # Currently only works with single amino acid substitutions
        loc = None
        if self.top_effected_protein is not None:
            pattern = re.compile(r'p\.[A-Z](\d+)[A-Z]')
            if pattern.match(self.top_effect.short_description):
                loc = int(re.sub(pattern, '\\g<1>', self.top_effect.short_description))
        return loc
            
