from os.path import isfile
from phenopackets import OntologyClass
from phenopackets import Phenopacket 
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
        if len(phenopack.diseases) != 0:
            dis = phenopack.diseases[0]
            self._diseases = Disease(dis.term.id, dis.term.label)
        elif isinstance(phenopack.interpretations[0].diagnosis.disease, OntologyClass):
            dis = phenopack.interpretations[0].diagnosis.disease
            self._diseases = Disease(dis.id, dis.label)
        else:
            self._diseases = None
        self._variant = Variant(ref = ref, phenopack = phenopack)
        self._protein = Protein(self._variant.top_effected_protein)
        
    def __get_hpids(self):
        hp_ids = []
        for x in self._phenopack.phenotypic_features:
             if not x.excluded:
                hp_ids.append(Phenotype(x.type.id, x.type.label))
        return hp_ids
      
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
            return [phenotype.id for phenotype in self._phenotype]
        else:
            return None
    
    @property
    def phenotype_labels(self):
        if self._phenotype is not None:
            return [phenotype.label for phenotype in self._phenotype]
        else:
            return None
    
    @property
    def phenotypes(self):
        return self._phenotype
    
    @property
    def variant(self):
        return self._variant

    @property
    def protein(self):
        return self._protein

    @property
    def gene(self):
        return self._variant.variant.genes[0]

    def describe(self):
        stats = pd.Series({
            "ID": self.id,
            "Disease ID": self.disease_id,
            "Disease Label": self.disease_label,
            "HPO IDs": str(self.phenotype_ids),
            "HPO Terms": str(self.phenotype_labels),
            "Variant": self.variant.variant_string,
            "Variant Type": str(self.variant.variant_type),
            "Gene Affected": self.gene.gene_name,
            "Gene ID": self.gene.gene_id,
            "Effect of Variant": self.variant.top_effect.short_description,
            "Transcript ID": self.variant.top_effect_transcript.id,
            "Protein Affected": self.protein.label,
            "Protein ID": self.protein.id
                })
        return stats.T