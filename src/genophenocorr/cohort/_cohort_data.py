from collections import Counter
class Cohort:
    """This class creates a collection of patients and makes it easier to determine overlapping diseases, 
    phenotypes, variants, and proteins among the patients. If a list of JSON files is given, it will
    add each file as a patient into the grouping.
    
    Attributes:
        all_patients (Sequence[Patient]): A set of all Patient objects in the Cohort
        all_phenotypes (Sequence[Phenotype]): A set of all Phenotype objects in the Cohort
        all_variants (Sequence[Variant]): A set of all Variant objects in the Cohort
        all_proteins (Sequence[ProteinMetadata]): A set of all ProteinMetadata objects in the Cohort
        all_transcripts (Sequence[string]): A set of all transcript IDs referenced in all the Variant objects
        total_patient_count (integer): The total number of Patient objects
    Methods:
        list_all_patients(): A list of all patient IDs 
        list_all_phenotypes(top:Optional[integer]): A list of all the top phenotype IDs (or all IDs if top is None) and the count of how many patients have it. 
        list_all_variants(top:Optional[integer]): A list of all the top variants (or all variants if top is None) and the count of how many patients have it. 
        list_all_proteins(top:Optional[integer]): A list of all the top protein IDs (or all IDs if top is None) and the count of how many patients have it. 
        list_data_by_tx(transcript:Optional[string]): A list and count of all the variants effects found for all transcripts or a given transcript if transcript is not None.
    """
    def __init__(self, patient_set, phenotype_set, variant_set, protein_set, counts_dict, recessive = False): 
        """Constructs all necessary attributes for a Cohort object
        
        Args:
            patient_set (Sequence[Patient]): A set of all Patient objects in the Cohort
            phenotype_set (Sequence[Phenotype]): A set of all Phenotype objects in the Cohort
            variant_set (Sequence[Variant]): A set of all Variant objects in the Cohort
            protein_set (Sequence[ProteinMetadata]): A set of all ProteinMetadata objects in the Cohort
            counts_dict (Dictionary{String, Counter()}): A Dictionary with counts for Phenotypes, Variant, and Proteins objects represented in all Patients
            recessive (boolean, Optional): True if the Cohort is focused on a recessive allele. Defaults to False.
            """
        self._patient_set = patient_set
        self._phenotype_set = phenotype_set
        self._protein_set = protein_set
        self._variant_set = variant_set
        self._all_counts_dict = counts_dict
        self._recessive = recessive

    @property
    def all_patients(self):
        """
        Returns:
            set: A set of all the Patient objects in the Cohort
        """
        return self._patient_set

    @property
    def all_phenotypes(self):
        """
        Returns:
            set: A set of all the Phenotype objects in the Cohort
        """
        return self._phenotype_set

    @property
    def all_variants(self):
        """
        Returns:
            set: A set of all the Variant objects in the Cohort
        """
        return self._variant_set

    @property
    def all_proteins(self):
        """
        Returns:
            set: A set of all the ProteinMetadata objects in the Cohort
        """
        return self._protein_set

    @property
    def all_transcripts(self):
        """
        Returns:
            set: A set of all the transcript IDs in the Cohort
        """
        all_trans = set()
        for var in self.all_variants:
            all_trans.update([trans.transcript_id for trans in var.tx_annotations])
        return all_trans

    @property
    def total_patient_count(self):
        """
        Returns:
            integer: The count of all the Patient objects
        """
        return self._all_counts_dict.get('patients')

    def list_all_patients(self):
        """
        Returns:
            list: A list of all the patient IDs in the Cohort
        """
        return [pat.patient_id for pat in self.all_patients]

    def list_all_phenotypes(self, top = None):
        """
        Args:
            top (integer, Optional): If not given, lists all phenotypes. Otherwise, lists only the `top` highest counts
        Returns:
            list: A list of tuples, formatted (phenotype ID, number of patients with that phenotype) 
        """
        return self._all_counts_dict.get('phenotypes').most_common(top)

    def list_all_variants(self, top = None):
        """
        Args:
            top (integer, Optional): If not given, lists all variants. Otherwise, lists only the `top` highest counts
        Returns:
            list: A list of tuples, formatted (variant string, number of patients with that variant) 
        """
        return self._all_counts_dict.get('variants').most_common(top)

    def list_all_proteins(self, top = None):
        """
        Args:
            top (integer, Optional): If not given, lists all proteins. Otherwise, lists only the `top` highest counts
        Returns:
            list: A list of tuples, formatted (protein ID, number of patients with that protein) 
        """
        return self._all_counts_dict.get('proteins').most_common(top)

    def list_data_by_tx(self, transcript = None):
        """
        Args:
            transcript (string, Optional): If not given, lists all transcripts. Otherwise, will only list the given transcript
        Returns:
            dictionary: Each transcript ID references a Counter(), with the variant effect as the key and total variants with that effect as the count value
        """
        if transcript is not None:
            var_type_dict = {transcript:Counter()}
        else:
            var_type_dict = {tx_id: Counter() for tx_id in self.all_transcripts}
        for var in self.all_variants:
            for trans in var.tx_annotations:
                if trans.transcript_id in var_type_dict:
                    var_type_dict.get(trans.transcript_id).update(trans.variant_effects)
        too_small = []
        for tx_id, var_effect_counter in var_type_dict.items():
            if len(var_effect_counter) <= 2:
                too_small.append(tx_id)
        for tx_id in too_small:
            del var_type_dict[tx_id]
        return var_type_dict
