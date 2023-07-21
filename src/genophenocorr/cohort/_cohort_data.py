from collections import Counter


class Cohort:
    # TODO(lnrekerle) - update the class doc string.
    """
    This class creates a collection of patients and makes it easier to determine overlapping diseases, 
    phenotypes, variants, and proteins among the patients. If a list of JSON files is given, it will
    add each file as a patient into the grouping, otherwise patients can be added individually with
    the self.add(Patient) function. 
    
    """

    def __init__(self, patient_set, phenotype_set, variant_set, protein_set, counts_dict, recessive=False):
        self._patient_set = patient_set
        self._phenotype_set = phenotype_set
        self._protein_set = protein_set
        self._variant_set = variant_set
        self._all_counts_dict = counts_dict
        self._recessive = recessive

    @property
    def all_patients(self):
        return self._patient_set

    @property
    def all_phenotypes(self):
        return self._phenotype_set

    @property
    def all_variants(self):
        return self._variant_set

    @property
    def all_proteins(self):
        return self._protein_set

    @property
    def all_transcripts(self):
        all_trans = set()
        for var in self.all_variants:
            all_trans.update([trans.transcript_id for trans in var.tx_annotations])
        return all_trans

    @property
    def total_patient_count(self):
        return self._all_counts_dict.get('patients')

    def list_all_patients(self):
        return [pat.patient_id for pat in self.all_patients]

    def list_all_phenotypes(self, top=None):
        return self._all_counts_dict.get('phenotypes').most_common(top)

    def list_all_variants(self, top=None):
        return self._all_counts_dict.get('variants').most_common(top)

    def list_all_proteins(self, top=None):
        return self._all_counts_dict.get('proteins').most_common(top)

    def list_data_by_tx(self, transcript=None):
        if transcript is not None:
            var_type_dict = {transcript: Counter()}
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
