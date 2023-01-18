from scipy import stats 
from statsmodels.stats.multitest import multipletests
import numpy as np
import pandas as pd

def has_hpo(pat, hpo):
    for h in pat.phenotype_ids:
        if h == hpo[0]:
            return True
    return False

# TODO change to run_stats
def RunStats(PatientList, Fun1, Fun2, extraVar_1 = None, extraVar_2 = None, show = 5, chi2 = False): 
    all_hpo = PatientList.list_all_phenotypes()
    print(f"Total hpo terms: {len(all_hpo)}")
    print(f"BonCorrected p-value = {0.05/len(all_hpo)}")
    allSeries = []
    for hpo_id in all_hpo:
        var1_with_hpo = len([ pat for pat in PatientList.all_patients.values() if has_hpo(pat, hpo_id) and Fun1(pat, extraVar_1)])
        var1_without_hpo = len([ pat for pat in  PatientList.all_patients.values() if not has_hpo(pat, hpo_id) and Fun1(pat,extraVar_1)])
        var2_with_hpo = len([ pat  for pat in PatientList.all_patients.values() if has_hpo(pat, hpo_id) and Fun2(pat,extraVar_2)])
        var2_without_hpo = len([ pat for pat in PatientList.all_patients.values() if not has_hpo(pat, hpo_id) and Fun2(pat,extraVar_2)])
        table = np.array([[var1_with_hpo, var1_without_hpo], [var2_with_hpo, var2_without_hpo]])
        if chi2:
            chi, p, dof, exp = stats.chi2_contingency(table)
        else:
            oddsr, p =  stats.fisher_exact(table, alternative='two-sided') 
        allSeries.append(pd.Series([var1_with_hpo, var1_without_hpo, var2_with_hpo, var2_without_hpo, p], name= hpo_id[0] + '-'+hpo_id[1], index=['1 w/ hpo', '1 w/o hpo', '2 w/ hpo', '2 w/o hpo', 'pval']))
    
    results = pd.concat(allSeries, axis=1)
    results = results.transpose()
    results = results.sort_values(['pval'])
    
    print(results.head(show))
    #pvalues = results.loc[:,'pval'].to_list()
    return results

        