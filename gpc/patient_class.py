from os.path import isfile
from phenopackets import Phenopacket
from google.protobuf.json_format import Parse
import varcode as vc
import json
import pyensembl
from .genotype_class import Genotype
from .disease_class import Disease
from .phenotype_class import Phenotype 

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
        if len(phenopack.interpretations) != 0:
            self._genotype = self.__get_variants()
        else:
            #print('No interpretations found')
            self._genotype = None
                
        
    def __get_hpids(self):
        hp_ids = []
        for x in self._phenopack.phenotypic_features:
             if not x.excluded:
                hp_ids.append(Phenotype(x))
        return hp_ids
    
    def __get_variants(self):
        interp = self._phenopack.interpretations[0]
        genotypes = []
        for geno in interp.diagnosis.genomic_interpretations:
            genotypes.append(Genotype(geno))
        return genotypes
            
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
        if self._genotype is not None:
            return [g.variant for g in self._genotype]
        else:
            return None

    @property
    def genotypes(self):
        return self._genotype
    
    @property
    def var_effect(self):
        if self._genotype is not None:
            return [g.top_var_effect for g in self._genotype]
        else:
            return None
    
    @property
    def var_is_missense(self):
        if self._genotype is not None:
            return [g.is_missense for g in self._genotype]
        else:
            return None
    
    @property
    def var_is_nonsense(self):
        if self._genotype is not None:
            return [g.is_nonsense for g in self._genotype]
        else:
            return None
        
    @property
    def var_is_deletion(self):
        if self._genotype is not None:
            return [g.is_deletion for g in self._genotype]
        elif self._diseases is not None:
            return [self._diseases.is_deletion]
        else:
            return [None]
        
    @property
    def var_is_duplication(self):
        if self._genotype is not None:
            return [g.is_duplication for g in self._genotype]
        elif self._diseases is not None:
            return [self._diseases.is_duplication]
        else:
            return [None]
    
    def describe(self):
        stats = {
            "ID": self.id,
            "Disease": self.disease_label,
            "Phenotypic Features": self.phenotype_labels,
            "Variant": self.variant,
            "Primary Effect of Variant": self.var_effect,
            "Is Missense Mutation?": self.var_is_missense,
            "Is Nonsense Mutation?": self.var_is_nonsense,
            "Is Duplication?": self.var_is_duplication,
            "Is Deletion?": self.var_is_deletion
        }
        return stats