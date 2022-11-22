import scipy
import numpy as np
import pandas as pd

def has_hpo(pat, hpo):
    for h in pat.phenotype_ids:
        if h == hpo[0]:
            return True
    return False

def RunStats(PatientList, Fun1, Fun2): 
    all_hpo = PatientList.list_all_phenotypes()
    print(f"Total hpo terms: {len(all_hpo)}")
    results = pd.DataFrame({'counts':['1 w/ hpo', '1 w/o hpo', '2 w/ hpo', '2 w/o hpo', 'pval']})
    for hpo_id in all_hpo:
        var1_with_hpo = len([ pat for pat in PatientList.all_patients.values() if has_hpo(pat, hpo_id) and Fun1(pat)])
        var1_without_hpo = len([ pat for pat in  PatientList.all_patients.values() if not has_hpo(pat, hpo_id) and Fun1(pat)])
        var2_with_hpo = len([ pat  for pat in PatientList.all_patients.values() if has_hpo(pat, hpo_id) and Fun2(pat)])
        var2_without_hpo = len([ pat for pat in PatientList.all_patients.values() if not has_hpo(pat, hpo_id) and Fun2(pat)])
        table = np.array([[var1_with_hpo, var1_without_hpo], [var2_with_hpo, var2_without_hpo]])
        oddsr, p =  scipy.stats.fisher_exact(table, alternative='two-sided') ##Add option for chi2
        results.insert(1, hpo_id[0] + '-'+hpo_id[1], [var1_with_hpo, var1_without_hpo, var2_with_hpo, var2_without_hpo, p])
    results = results.set_index('counts')
    return results

        