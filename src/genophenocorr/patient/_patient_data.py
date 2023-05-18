



class Patient:
    def __init__(self, patient_id, phenotypes, variants, proteins, disease):
        self._id = patient_id
        self._phenotypes = phenotypes
        self._variants = variants
        self._proteins = proteins
        self._diseases = disease

    @property
    def patient_id(self):
        return self._id

    @property
    def diseases(self):
        return self._diseases
    
    @property
    def phenotypes(self):
        return self._phenotypes
    
    @property
    def variants(self):
        return self._variants

    @property
    def proteins(self):
        return self._proteins


    # def has_hpo(self, hpo, all_hpo):
    #     if not isinstance(all_hpo, defaultdict):
    #         for h in self.phenotype_ids:
    #             if h == hpo:
    #                 return True
    #         return False
    #     else:
    #         for h in self.phenotype_ids:
    #             if h in all_hpo.get(hpo):
    #                 return True
    #         return False


    def __str__(self) -> str:
        return f"Patient(patient_id:{self.patient_id}, " \
            f"diseases:{[dis.identifier for dis in self.diseases]}, " \
            f"variants:{[var.variant_string for var in self.variants]}, " \
            f"phenotypes:{[pheno.identifier for pheno in self.phenotypes]}, " \
            f"proteins:{[prot.protein_id for prot in self.proteins]})"

    def __repr__(self) -> str:
        return str(self)

    def __eq__(self, other) -> bool:
        return isinstance(other, Patient) \
            and self.patient_id == other.patient_id \
            and self.diseases == other.diseases \
            and self.variants == other.variants \
            and self.phenotypes == other.phenotypes \
            and self.proteins == other.proteins
    
    def __hash__(self) -> int:
        return hash((self.patient_id, self.diseases, self.variants, self.phenotypes, self.proteins))