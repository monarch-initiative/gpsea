import varcode as vc
import pyensembl

class Variant:
    def __init__(self, phenopack = None):
        if phenopack is not None:
            self._phenopack = phenopack
            self._variant = self.__find_variant()
        else:
            self._variant = None
        
    def __find_variant(self):
        if len(self._phenopack.interpretations) != 0:
            Interp = self._phenopack.interpretations[0]
            genInterp = Interp.variant_interpretation.variation_descriptor.vcf_record
            contig = genInterp.chrom
            start = genInterp.pos
            ref = genInterp.ref
            alt = genInterp.alt
            myVar = vc.Variant(str(contig), start, ref, alt, ensembl = pyensembl.ensembl_grch37)
            return myVar
        else:
            return None

    @property
    def variant(self):
        return self._variant

    @property
    def variant_string(self):
        if self._variant is not None:
            return self.variant.short_description
        else:
            return None

    @property
    def effect(self):
        if self._variant is not None:
            return self.variant.effects().top_priority_effect()
        else:
            return None
    