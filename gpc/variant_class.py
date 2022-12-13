import varcode as vc
import pyensembl
from .proteins_class import Protein

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
            try:
                genInterp = Interp.diagnosis.genomic_interpretations[0]
                varInterp = genInterp.variant_interpretation.variation_descriptor.vcf_record
                contig = varInterp.chrom
                start = varInterp.pos
                ref = varInterp.ref
                alt = varInterp.alt
                myVar = vc.Variant(str(contig), start, ref, alt, ensembl = pyensembl.ensembl_grch37)
                return myVar
            except AttributeError:
                print('Could not find Variants')
                return None
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
    def effects(self):
        if self._variant is not None:
            return self.variant.effects()
        else:
            return None

    @property
    def top_effect(self):
        if self.effects is not None:
            return self.effects.top_priority_effect()
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
        if self.transcripts is not None:
            return self.variant.transcript_ids
        else:
            return None

    @property
    def top_effect_transcript(self):
        if self.top_effect is not None:
            return self.top_effect.transcript
        else:
            return None
    
    @property
    def protein(self):
        if self.top_effect_transcript is not None:
            if self.top_effect_transcript.is_protein_coding:
                return Protein(self.top_effect_transcript)
            else:
                print('WARNING - Transcript associated with top effect does not code for a protein.')
                print('Will use top coding effect for protein analysis.')
                effect = self.variant.effects().drop_silent_and_noncoding().top_priority_effect()
                return Protein(effect.transcript)