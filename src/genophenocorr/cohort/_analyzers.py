from ._cohort_data import Cohort


class CohortAnalysis():

    def __init__(self, cohort, recessive = False, include_unmeasured = True,  include_large_SV = True, min_perc_patients_w_hpo = 10) -> None:
        if not isinstance(cohort, Cohort):
            raise ValueError(f"cohort must be type Cohort but was type {type(cohort)}")
        self._cohort = cohort
        self._recessive = recessive
        self._in_unmeasured = include_unmeasured
        self._include_SV = include_large_SV
        if not isinstance(min_perc_patients_w_hpo, int) or min_perc_patients_w_hpo <= 0 or min_perc_patients_w_hpo > 100:
            raise ValueError(f"hpo_percent_of_patients must be an integer between 1 & 100 but was {min_perc_patients_w_hpo}")
        self._percent_pats = min_perc_patients_w_hpo
        self._testing_hpo_terms = self._remove_low_hpo_terms()
        
    @property
    def analysis_cohort(self):
        return self._cohort

    @property
    def is_recessive(self):
        return self._recessive

    @property
    def include_unmeasured(self):
        return self._in_unmeasured

    @property
    def include_large_SV(self):
        return self._include_SV

    @property
    def min_perc_patients_w_hpo(self):
        return self._percent_pats

    def _remove_low_hpo_terms(self):
        final_hpo_list = set()
        if self.include_large_SV:
            all_patients = self.analysis_cohort.all_patients
        else:
            all_patients = [pat for pat in self.analysis_cohort.all_patients if not any(var.variant_coordinates.is_structural for var in pat.variants)]
        for hpo in self.analysis_cohort.all_phenotypes:
            total_pats_with = 0
            total_pats_measured = 0
            for pat in all_patients:
                if hpo in pat.phenotypes:
                    if hpo.observed == True:
                        total_pats_with += 1
                        total_pats_measured += 1
                    elif hpo.observed == False:
                        total_pats_measured += 1
                elif self.include_unmeasured:
                    total_pats_measured += 1
            if (total_pats_with / total_pats_measured) >= (self.min_perc_patients_w_hpo / 100):
                final_hpo_list.add(hpo)
        if len(final_hpo_list) == 0:
            raise ValueError(f"No HPO terms found in over {self.min_perc_patients_w_hpo}% of patients.")
        return final_hpo_list

    def compare_by_variant_type(self, var_type1, var_type2 = None):
        final_dict = dict()
        if not self.is_recessive:
            for hpo in self._testing_hpo_terms:
                with_hpo = [pat for pat in self.analysis_cohort.all_patients if hpo in pat.phenotypes]
                not_hpo = [pat for pat in self.analysis_cohort.all_patients if hpo not in pat.phenotypes]
                with_hpo_with_var1_count = len([pat for pat in with_hpo if var_type1 in [var.selected_transcript.variant_effects for var in pat.variants]])
                not_hpo_with_var1_count = len([pat for pat in not_hpo if var_type1 in [var.selected_transcript.variant_effects for var in pat.variants]])
                if var_type2 is None:
                    with_hpo_with_var2_count = len([pat for pat in with_hpo if var_type1 not in [var.selected_transcript.variant_effects for var in pat.variants]])
                    not_hpo_with_var2_count = len([pat for pat in not_hpo if var_type1 not in [var.selected_transcript.variant_effects for var in pat.variants]])
                else:
                    with_hpo_with_var2_count = len([pat for pat in with_hpo if var_type2 in [var.selected_transcript.variant_effects for var in pat.variants]])
                    not_hpo_with_var2_count = len([pat for pat in not_hpo if var_type2 in [var.selected_transcript.variant_effects for var in pat.variants]])
                final_dict[hpo] = ((with_hpo_with_var1_count, not_hpo_with_var1_count), (with_hpo_with_var2_count, not_hpo_with_var2_count))
               
        return final_dict

    # def _count_patients_with(self, variable, without = False):


    #     # if not recessive:
    #     #     allSeries = self.__dominant_test(testing_patients, compare_hpos_d, hpo_counts, remove_not_measured, Function_1, Func1_Variable, Function_2, Func2_Variable)
    #     # elif recessive:
    #     #     allSeries = self.__recessive_test(testing_patients, compare_hpos_d, hpo_counts, remove_not_measured, Function_1, Func1_Variable, Function_2, Func2_Variable)
    #     # results = pd.concat(allSeries, axis=1)
    #     # results = results.transpose() 
    #     # results = results.sort_values(['pval'])
    #     # pval_adjusted = multipletests(results['pval'].to_list(), alpha=0.05, method=adjusted_pval_method) 
    #     # results['adjusted pval'] = np.array(pval_adjusted[1])
    #     # results = results.convert_dtypes()
        
    #     # return results

    # def __dominant_test(self, testing_hpos, ):
    #     allSeries = []
    #     for hpo_id in compare_hpos_d:
    #         if remove_not_measured:
    #             tested_patients_d = [patient for patient in testing_patients_d if hpo_id in patient[1].phenotypes.keys()]
    #         else:
    #             tested_patients_d = testing_patients_d
    #         var1_with_hpo = len([ pat for pat in tested_patients_d if pat[1].has_hpo(hpo_id, compare_hpos_d) and True in [Function_1(var, pat[1], Func1_Variable) for var in pat[1].variants]])
    #         var1_without_hpo = len([ pat for pat in  tested_patients_d if not pat[1].has_hpo(hpo_id, compare_hpos_d) and True in [Function_1(var, pat[1], Func1_Variable) for var in pat[1].variants]])
    #         var2_with_hpo = len([ pat  for pat in tested_patients_d if pat[1].has_hpo(hpo_id, compare_hpos_d) and True in [Function_2(var, pat[1], Func2_Variable) for var in pat[1].variants]])
    #         var2_without_hpo = len([ pat for pat in tested_patients_d if not pat[1].has_hpo(hpo_id, compare_hpos_d) and True in [Function_2(var, pat[1], Func2_Variable) for var in pat[1].variants]])
    #         table = np.array([[var1_with_hpo, var1_without_hpo], [var2_with_hpo, var2_without_hpo]])
    #         oddsr, p =  stats.fisher_exact(table, alternative='two-sided') 
    #         allSeries.append(pd.Series([var1_with_hpo, var1_without_hpo, var2_with_hpo, var2_without_hpo, p], name= hpo_id + ' - ' + hpo_counts.at[hpo_id, 'Class'].label, index=['1 w/ hpo', '1 w/o hpo', '2 w/ hpo', '2 w/o hpo', 'pval']))
    #     return allSeries

    # # def __recessive_test(self, testing_patients_d, compare_hpos_d, hpo_counts, remove_not_measured, Function_1, Func1_Variable, Function_2, Func2_Variable):
    # #     allSeries = []
    # #     for hpo_id in compare_hpos_d:
    # #         AA_with_hpo = 0
    # #         AA_without_hpo = 0
    # #         AB_with_hpo = 0
    # #         AB_without_hpo = 0
    # #         BB_with_hpo = 0
    # #         BB_without_hpo = 0
    # #         if remove_not_measured:
    # #             tested_patients_d = [patient for patient in testing_patients_d if hpo_id in patient[1].phenotypes.keys()]
    # #         else:
    # #             tested_patients_d = testing_patients_d
    # #         for pat in tested_patients_d:
    # #             if len(pat[1].variants) == 2:
    # #                 var1 = pat[1].variants[0]
    # #                 var2 = pat[1].variants[1]
    # #                 if pat[1].has_hpo(hpo_id, compare_hpos_d) and (Function_1(var1, pat[1], Func1_Variable) != Function_1(var2, pat[1], Func1_Variable)) and (Function_2(var1, pat[1], Func2_Variable) != Function_2(var2, pat[1], Func2_Variable)):
    # #                     AB_with_hpo += 1
    # #                 if not pat[1].has_hpo(hpo_id, compare_hpos_d) and (Function_1(var1, pat[1], Func1_Variable) != Function_1(var2, pat[1], Func1_Variable)) and (Function_2(var1, pat[1], Func2_Variable) != Function_2(var2, pat[1], Func2_Variable)):
    # #                     AB_without_hpo += 1
    # #             elif len(pat[1].variants) == 1:
    # #                 var1 = pat[1].variants[0]
    # #                 var2 = pat[1].variants[0]
    # #             if pat[1].has_hpo(hpo_id, compare_hpos_d) and Function_1(var1, pat[1], Func1_Variable) and Function_1(var2, pat[1], Func1_Variable):
    # #                 AA_with_hpo += 1
    # #             if not pat[1].has_hpo(hpo_id, compare_hpos_d) and Function_1(var1, pat[1], Func1_Variable) and Function_1(var2, pat[1], Func1_Variable):
    # #                 AA_without_hpo += 1
    # #             if pat[1].has_hpo(hpo_id, compare_hpos_d) and Function_2(var1, pat[1],Func2_Variable) and Function_2(var2, pat[1], Func2_Variable):
    # #                 BB_with_hpo += 1
    # #             if not pat[1].has_hpo(hpo_id, compare_hpos_d) and Function_2(var1, pat[1],Func2_Variable) and Function_2(var2, pat[1], Func2_Variable):
    # #                 BB_without_hpo += 1
    # #         table = np.array([[AA_with_hpo, AA_without_hpo],[AB_with_hpo, AB_without_hpo], [BB_with_hpo, BB_without_hpo]])
    # #         p = fisher_exact(table)
    # #         #oddsr, p =  stats.fisher_exact(table, alternative='two-sided') 
    # #         allSeries.append(pd.Series([AA_with_hpo, AA_without_hpo, AB_with_hpo, AB_without_hpo, BB_with_hpo, BB_without_hpo, p], name= hpo_id + ' - ' + hpo_counts.at[hpo_id, 'Class'].label, index=['AA w/ hpo', 'AA w/o hpo', 'AB w/ hpo', 'AB w/o hpo', 'BB w/ hpo', 'BB w/o hpo', 'pval']))
    # #     return allSeries

