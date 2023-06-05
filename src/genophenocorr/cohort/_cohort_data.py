import tabulate

class Cohort:
    """
    This class creates a collection of patients and makes it easier to determine overlapping diseases, 
    phenotypes, variants, and proteins among the patients. If a list of JSON files is given, it will
    add each file as a patient into the grouping, otherwise patients can be added individually with
    the self.add(Patient) function. 
    
    """
    def __init__(self, patient_set, phenotype_dict, variant_dict, protein_dict, recessive = False):
        self._patient_set = patient_set
        self._phenotype_dict = phenotype_dict
        self._protein_dict = protein_dict
        self._variant_dict = variant_dict
        self._recessive = recessive

    @property
    def all_patients(self):
        return self._patient_set

    @property
    def all_phenotypes(self):
        return self._phenotype_dict.keys()

    @property
    def all_variants(self):
        return self._variant_dict.keys()

    @property
    def total_patient_count(self):
        return len(self.all_patients)

    @property
    def all_proteins(self):
        return self._protein_dict.keys()

    def list_all_patients(self):
        return [pat.patient_id for pat in self.all_patients]

    def list_all_phenotypes(self):
        headers = ["HPO ID", "Total Patients"]
        print(tabulate.tabulate(sorted([(k.identifier.value, v) for k, v in self._phenotype_dict.items() if k.observed == True], key=lambda row: row[1], reverse=True), headers = headers))

    def list_all_variants(self):
        headers = ["Variant", "Total Patients"]
        print(tabulate.tabulate(sorted([(k.variant_string, v) for k, v in self._variant_dict.items()], key=lambda row: row[1], reverse=True), headers = headers))

    def list_all_proteins(self):
        headers = ["Protein ID", "Total Patients"]
        print(tabulate.tabulate(sorted([(k.protein_id, v) for k, v in self._protein_dict.items()], key=lambda row: row[1], reverse=True), headers = headers))
