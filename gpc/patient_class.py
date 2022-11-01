from os.path import isfile
from phenopackets import Phenopacket
from google.protobuf.json_format import Parse
import varcode as vc
import json
import pyensembl


class Patient:
    def __init__(self, phenopackJson):
        if not isfile(phenopackJson):
            raise FileNotFoundError("Could not find phenopacket")
            
        with open(phenopackJson) as f:
            data = f.read()
        jsondata = json.loads(data)
        phenopack = Parse(json.dumps(jsondata), Phenopacket())
        
        self._phenopack = phenopack
        self._phenotype = self.__get_hpids()
        self._genotype = []
        if len(phenopack.interpretations) != 0:
            for i in range(len(phenopack.interpretations)):
                self._genotype.append(self.__get_variants(i))
        else:
            print('WARNING: No interpretations found for: '+ phenopackJson )
            self._genotype = None
                
        
    
    def __get_hpids(self):
        hp_ids = [x.type for x in self._phenopack.phenotypic_features if not x.excluded]
        return hp_ids
    
    def __get_variants(self, n):
        interp = self._phenopack.interpretations[n] 
        contig = []
        start = []
        ref = []
        alt = []
        for i in interp.diagnosis.genomic_interpretations:
            try:
                var_des = i.variant_interpretation.variation_descriptor.vcf_record
                contig.append(int(var_des.chrom.split('_')[1].split('.')[0]))
                ref.append(var_des.ref)
                alt.append(var_des.alt)
                start.append(var_des.pos)
            except:
                print("WARNING: Variant has no vcf_record. Cannot be processed.")
                continue
        myVars = []
        for i in range(len(contig)):
            myVar = vc.Variant(str(contig[i]), start[i], ref[i], alt[i], ensembl = pyensembl.ensembl_grch38)
            myVars.append(myVar)
        return myVars
    
    @property
    def get_phenopacket(self):
        return self._phenopack
    
    @property
    def get_phenotypes(self):
        return self._phenotype
    
    @property
    def get_genotypes(self):
        geno = []
        if self._genotype == None:
            return None
        for i in self._genotype:
            for e in i:
                geno.append(e.short_description)
        return geno
    
    @property
    def get_var_effects(self):
        if self._genotype == None:
            return None
        effected = []
        for i in self._genotype:
            for e in i:
                effected.extend(e.effects().effects)
        return effected
    
    def is_missense(self):
        if self._genotype == None:
            return None
        miss = []
        for i in self.get_var_effects:
            if i.short_description.endswith("*") or not i.variant.is_snv:
                miss.append(False)
            else:
                miss.append(True)
        return miss
    
    def is_nonsense(self):
        if self._genotype == None:
            return None
        non = []
        for i in self.get_var_effects():
            if i.short_description.endswith("*") and i.variant.is_snv:
                non.append(True)
            else:
                non.append(False)
        return non

    def describe(self):
        stats = {
            "ID": self._phenopack.id,
            "Phenotypic Features": self.get_phenotypes,
            "Variants": self.get_genotypes,
            "Effects of Variants": self.get_var_effects,
            "Number Missense": sum(self.is_missense())
        }
        return stats