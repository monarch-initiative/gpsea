import scipy
import numpy as np

def has_hpo(pat, hpo):
    for h in pat.phenotype_ids:
        if h == hpo[0]:
            return True
    return False

def is_vartype(pat, vartype):
    if vartype == 'missense':
        return pat.var_is_missense[0]
    elif vartype == 'nonsense':
        return pat.var_is_nonsense[0]
    elif vartype == 'deletion':
        if pat.var_is_deletion[0] is not None:
            return pat.var_is_deletion[0]
        else:
            return False
    elif vartype == 'duplication':
        if pat.var_is_duplication[0] is not None:
            return pat.var_is_duplication[0]
        else:
            return False


def RunStats(PatientList, VarType1, VarType2):
    all_hpo = PatientList.list_all_phenotypes()
    print(f"Total hpo terms: {len(all_hpo)}")
    for hpo_id in all_hpo:
        var1_with_hpo = len([ pat for pat in PatientList.all_patients.values() if has_hpo(pat, hpo_id) and is_vartype(pat, VarType1)])
        var1_without_hpo = len([ pat for pat in  PatientList.all_patients.values() if not has_hpo(pat, hpo_id) and is_vartype(pat, VarType1)])
        var2_with_hpo = len([ pat  for pat in PatientList.all_patients.values() if has_hpo(pat, hpo_id) and is_vartype(pat, VarType2)])
        var2_without_hpo = len([ pat for pat in PatientList.all_patients.values() if not has_hpo(pat, hpo_id) and is_vartype(pat, VarType2)])
        print(f"HPO {hpo_id}:\n {VarType1} with HPO: {var1_with_hpo}\n {VarType1} without HPO: {var1_without_hpo}\n {VarType2} with HPO: {var2_with_hpo}\n {VarType2} without HPO: {var2_without_hpo}")
        table = np.array([[var1_with_hpo, var1_without_hpo], [var2_with_hpo, var2_without_hpo]])
        oddsr, p =  scipy.stats.fisher_exact(table, alternative='two-sided') ##Add option for chi2
        print(f"p for {hpo_id}: {p}")
        