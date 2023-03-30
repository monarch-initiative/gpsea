from os.path import isfile
from phenopackets import OntologyClass
from phenopackets import Phenopacket 
from collections import defaultdict
from google.protobuf.json_format import Parse
import json
import pyensembl
import pandas as pd
from .disease import Disease
from .phenotype import Phenotype 
from .variant import Variant
from .proteins import Protein


class Patient:
    def __init__(self, phenopackJson, ref):
        if not isfile(phenopackJson):
            raise FileNotFoundError("Could not find phenopacket")
            
        with open(phenopackJson) as f:
            data = f.read()
        jsondata = json.loads(data)
        phenopack = Parse(json.dumps(jsondata), Phenopacket())
        
        self._id = phenopack.id
        self._phenopack = phenopack
        self._phenotype = self.__get_hpids()
        self._reference = ref
        self._protein = []
        if len(phenopack.diseases) != 0:
            dis = phenopack.diseases[0]
            self._diseases = Disease(dis.term.id, dis.term.label)
        elif len(phenopack.interpretations) != 0:
            if len(phenopack.interpretations[0].diagnosis.disease.id) > 0:
                dis = phenopack.interpretations[0].diagnosis.disease
                self._diseases = Disease(dis.id, dis.label)
            else:
                self._diseases = None
        else:
            self._diseases = None
        self._variants = self.__get_vars()
        for var in self._variants:
            self._protein.append(Protein(var.top_effected_protein))
        
    def __get_hpids(self):
        hp_ids = defaultdict(Phenotype)
        for x in self._phenopack.phenotypic_features:
            if not x.excluded:
                hp_ids[x.type.id] = Phenotype(x.type.id)
            elif x.excluded:
                hp_ids[x.type.id] = Phenotype(x.type.id, excluded=True)
        return hp_ids

    def __get_vars(self):
        allVars = []
        if len(self._phenopack.interpretations) > 0:
            Interp = self._phenopack.interpretations[0]
            if len(Interp.diagnosis.genomic_interpretations) > 0:
                for genoInterp in Interp.diagnosis.genomic_interpretations:
                    allVars.append(Variant(ref = self._reference, genoInterp = genoInterp))
            else: raise ValueError('No genomic interpretations found in phenopacket.')
        else: ValueError('No interpretations found in phenopacket.')
        return allVars

      
    @property
    def id(self):
        return self._phenopack.id

    @property
    def disease_id(self):
        if self._diseases is not None:
            return self._diseases.id
        else:
            return None

    @property
    def hg_reference(self):
        return self._reference
    
    @property
    def disease_label(self):
        if self._diseases is not None:
            return self._diseases.label
        else:
            return None

    @property
    def diseases(self):
        return self._diseases
    
    @property
    def phenopacket(self):
        return self._phenopack
    
    @property
    def phenotype_ids(self):
        if self._phenotype is not None:
            return [phenotype.id for phenotype in self._phenotype.values() if not phenotype.excluded]
        else:
            return None
    
    @property
    def phenotype_labels(self):
        if self._phenotype is not None:
            return [phenotype.label for phenotype in self._phenotype.values() if not phenotype.excluded]
        else:
            return None
    
    @property
    def phenotypes(self):
        return self._phenotype
    
    @property
    def variants(self):
        return self._variants

    @property
    def variant_strings(self):
        return [var.variant_string for var in self.variants]
    
    @property
    def variant_types(self):
        return [var.variant_types for var in self.variants]

    @property
    def proteins(self):
        return self._protein

    @property
    def protein_ids(self):
        return [p.id for p in self.proteins]

    @property
    def genes(self):
        return [var.variant.genes for var in self.variants]

    def get_patient_description_df(self): 
        stats = pd.Series({
            "ID": self.id,
            "Disease ID": self.disease_id,
            "Disease Label": self.disease_label,
            "HPO IDs": str(self.phenotype_ids),
            "HPO Terms": str(self.phenotype_labels),
            "Variant": self.variant_strings,
            "Variant Type": self.variant_types,
            "Gene Affected": [g.gene_name for gs in self.genes for g in gs],
            "Gene ID": [g.gene_id for gs in self.genes for g in gs],
            "Effect of Variant": [v.top_effect.short_description for v in self.variants],
            "Transcript ID": [v.top_effect_transcript.id for v in self.variants],
            "Protein Affected": [p.label for p in self.proteins],
            "Protein ID": [p.id for p in self.proteins]
                })
        return stats.T


    def has_hpo(self, hpo, all_hpo):
        if not isinstance(all_hpo, defaultdict):
            for h in self.phenotype_ids:
                if h == hpo:
                    return True
            return False
        else:
            for h in self.phenotype_ids:
                if h in all_hpo.get(hpo):
                    return True
            return False