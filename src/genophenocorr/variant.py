from curses.ascii import isdigit
import re
import requests
import pandas as pd

URL = 'https://rest.ensembl.org/vep/human/hgvs/%s:%s?protein=1&variant_class=1&numbers=1&LoF=1&canonical=1&domains=1&hgvs=1&mutfunc=1&refseq=1&transcript_version=1'

class Variant:
    def __init__(self, genoInterp, transcript):
        self._genoInterp = genoInterp
        self._variant_json = self.__find_variant()
        self.__set_transcript(transcript)

    def __find_variant(self):
        varInterp = self._genoInterp.variant_interpretation.variation_descriptor
        self._allelic_state = varInterp.allelic_state.label
        chrom = varInterp.vcf_record.chrom
        contig = re.sub(r'[^0-9MXY]', '', chrom)
        if len(contig) == 0 or (contig.isdigit() and (int(contig) == 0 or int(contig) >= 24)):
            ## Chromosome can only be values 1-23, X, Y, or M
            raise ValueError(f"Chromosome unknown: {chrom}")
        HGVS = varInterp.expressions[1]
        var_input = HGVS.value.split(':')[1]
        api_url = URL % (chrom, var_input)
        r = requests.get(api_url, headers={'Content-Type':'application/json'})
        if not r.ok:
            r.raise_for_status()
        varJson = r.json()
        return varJson[0]

## TODO: Change All Variable to work with new JSON version
    @property
    def variant_json(self):
        return self._variant_json

    @property
    def variant_string(self):
        return self.variant_json.get('id')

    @property
    def genotype(self):
        if self._allelic_state is not None:
            return self._allelic_state
        else:
            return None

    @property
    def genomic_location(self):
        return self._variant_json.get('start')

    @property
    def variant_types(self):
        all_types = []
        all_types.append(self._variant_json.get('variant_class'))
        for typ in self._canonical.get('consequence_terms'):
            all_types.append(typ)
        return all_types

    @property
    def transcript(self):
        return self._transcript

    @property
    def gene_name(self):
        return self._canonical.get('gene_symbol')

    @property
    def effected_protein(self):
        return self._canonical.get('protein_id')

    @property
    def protein_effect_location(self):
        return self._canonical.get('protein_start')

    @property
    def protein_effect(self):
        return self._canonical.get('hgvsp').split(':')[1]


    def list_all_variant_effects_df(self):
        allSeries = []
        varID = self._variant_json.get('id')
        varClass = self._variant_json.get('variant_class')
        for trans in self._variant_json.get('transcript_consequences'):
            transID = trans.get('transcript_id')
            if not transID.startswith('NM'):
                continue
            exonEffected = trans.get('exon').split('/')[0]
            proteinEffect = trans.get('hgvsp').split(':')[1]
            consequence = trans.get('consequence_terms')
            geneName = trans.get('gene_symbol')
            proteinID = trans.get('protein_id')
            siftScore = trans.get('sift_score')
            polyphenScore = trans.get('polyphen_score')
            if trans.get('canonical') == 1:
                canon = True
            else:
                canon = False
            allSeries.append(pd.Series([consequence, geneName, proteinID, proteinEffect, exonEffected, siftScore, polyphenScore, canon], name = transID, index=['Variant Type', 'Gene Name', 'Protein ID', 'Effect on Protein', 'Exon Effected', 'Sift Score', 'Polyphen Score', 'Current Transcript']))
        results = pd.concat(allSeries, axis=1)
        results = results.transpose()
        results.Name = varID + ' - Class: ' + varClass
        return results 
            
    def __set_transcript(self, transcript = None):
        if transcript is not None:
            self._transcript = transcript
            for trans in self._variant_json.get('transcript_consequences'):
                if trans.get('transcript_id') == transcript:
                    self._canonical = trans
        else:
            for trans in self._variant_json.get('transcript_consequences'):
                if trans.get('canonical') == 1:
                    self._transcript = trans.get('transcript_id')
                    self._canonical = trans
        return None