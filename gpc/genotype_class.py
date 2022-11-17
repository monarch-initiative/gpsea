import pyensembl
import varcode as vc

class Genotype:
    def __init__(self, genInterp):
        contig = genInterp.variant_interpretation.variation_descriptor.vcf_record.chrom
        start = genInterp.variant_interpretation.variation_descriptor.vcf_record.pos
        ref = genInterp.variant_interpretation.variation_descriptor.vcf_record.ref
        alt = genInterp.variant_interpretation.variation_descriptor.vcf_record.alt
        self._myVar = vc.Variant(str(contig), start, ref, alt, ensembl = pyensembl.ensembl_grch37)
        self._variant = self._myVar.short_description
        self._topVar = self._myVar.effects().top_priority_effect()
        
    @property
    def variant(self):
        return self._variant
    
    @property
    def top_var_effect(self):
        return self._topVar.short_description
    
    @property
    def is_missense(self):
        if self._topVar.short_description.endswith("*") or not self._myVar.is_snv:
            return False
        else:
            return True
    
    @property
    def is_nonsense(self):
        if self._topVar.short_description.endswith("*") and self._myVar.is_snv:
            return True
        else:
            return False
        
    @property 
    def is_deletion(self):
        return self._myVar.is_deletion
    
    @property
    def is_insertion(self):
        return self._myVar.is_insertion
    
    @property
    def is_transition(self):
        return self._myVar.is_transition
    
    @property
    def is_transversion(self):
        return self._myVar.is_transversion
    
    @property
    def is_indel(self):
        return self._myVar.is_indel
    
    @property
    def is_duplication(self):
        if 'dup' in self._myVar.short_description:
            return True
        else:
            return False
    
