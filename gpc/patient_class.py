from os.path import isfile
from phenopackets import Phenopacket
from google.protobuf.json_format import Parse
import json
import pyensembl
from .disease_class import Disease
from .phenotype_class import Phenotype 
from .variant_class import Variant

class Patient:
    def __init__(self, phenopackJson):
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
            self._diseases = Disease(phenopack.diseases[0])
        else:
            self._diseases = None
        self._variant = Variant(phenopack)
        
    def __get_hpids(self):
        hp_ids = []
        for x in self._phenopack.phenotypic_features:
             if not x.excluded:
                hp_ids.append(Phenotype(x))
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

    def describe(self):
        stats = {
            "ID": self.id,
            "Disease": self.disease_label,
            "Phenotypic Features": self.phenotype_labels,
            "Variant": self.variant.variant_string,
            "Effect of Variant": self.variant.effects
        }
        return stats