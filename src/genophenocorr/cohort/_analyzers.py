import logging
from ._cohort_data import Cohort
import typing
from genophenocorr.constants import VariantEffect
from genophenocorr.protein import FeatureType, ProteinFeature
from genophenocorr.predicate import VariantEffectPredicate, HPOPresentPredicate, VariantPredicate, ExonPredicate, ProtFeatureTypePredicate, ProtFeaturePredicate
from genophenocorr.variant import Variant
from scipy import stats
from pandas import DataFrame
from collections import Counter, namedtuple


class CohortAnalysis():

    def __init__(self, cohort, transcript, recessive = False, include_unmeasured = True,  include_large_SV = True, min_perc_patients_w_hpo = 10) -> None:
        if not isinstance(cohort, Cohort):
            raise ValueError(f"cohort must be type Cohort but was type {type(cohort)}")
        ## IF TRANSCRIPT DOES NOT EXIST, GIVE ERROR
        self._cohort = cohort
        self._transcript = transcript
        self._recessive = recessive
        self._in_unmeasured = include_unmeasured
        self._include_SV = include_large_SV
        self._patient_list = cohort.all_patients if include_large_SV else [pat for pat in cohort.all_patients if not all([var.variant_coordinates.is_structural() for var in pat.variants])]
        if not isinstance(min_perc_patients_w_hpo, int) or min_perc_patients_w_hpo <= 0 or min_perc_patients_w_hpo > 100:
            raise ValueError(f"hpo_percent_of_patients must be an integer between 1 & 100 but was {min_perc_patients_w_hpo}")
        self._percent_pats = min_perc_patients_w_hpo
        self._hpo_present_test = HPOPresentPredicate()
        self._testing_hpo_terms = self._remove_low_hpo_terms()
        self._patients_by_hpo = self._sort_patients_by_hpo()
        self._logger = logging.getLogger(__name__)
        
    @property
    def analysis_cohort(self):
        return self._cohort
    
    @property
    def analysis_patients(self):
        return self._patient_list
    
    @property
    def analysis_transcript(self):
        return self._transcript

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
        all_hpo_counts = Counter()
        all_no_hpo_counts = Counter()
        for pat in self.analysis_patients:
            all_hpo_counts.update([hpo for hpo in pat.phenotypes if hpo.observed == True])
            if not self.include_unmeasured:
                all_no_hpo_counts.update([hpo for hpo in pat.phenotypes if hpo.observed == False])
        if self.include_unmeasured:
            for hpo, count in all_hpo_counts.items():
                count = count if not None else 0
                if (count / len(self.analysis_patients)) >= (self.min_perc_patients_w_hpo / 100):
                    final_hpo_list.add(hpo)
        else:
            for hpo, count in all_hpo_counts.items():
                count = count if count is not None else 0
                no_count = all_no_hpo_counts.get(hpo) if all_no_hpo_counts.get(hpo) is not None else 0
                if (count / (count + no_count)) >= (self.min_perc_patients_w_hpo / 100):
                    final_hpo_list.add(hpo)
        if len(final_hpo_list) == 0:
            raise ValueError(f"No HPO terms found in over {self.min_perc_patients_w_hpo}% of patients.")
        return final_hpo_list

    def _sort_patients_by_hpo(self):
        all_with_hpo = {}
        all_without_hpo = {}
        for hpo in self._testing_hpo_terms:
            if self.include_unmeasured:
                with_hpo = [pat for pat in self.analysis_patients if self._hpo_present_test.test(pat, hpo).name == 'Observed']
                not_hpo = [pat for pat in self.analysis_patients if self._hpo_present_test.test(pat, hpo).name in ('NotObserved', 'NotMeasured')]
            else:
                with_hpo = [pat for pat in self.analysis_patients if self._hpo_present_test.test(pat, hpo).name == 'Observed']
                not_hpo = [pat for pat in self.analysis_patients if self._hpo_present_test.test(pat, hpo).name == 'NotObserved']
            all_with_hpo[hpo] = with_hpo
            all_without_hpo[hpo] = not_hpo
        patientsByHPO = namedtuple('patientByHPO', field_names=['all_with_hpo', 'all_without_hpo'])
        return patientsByHPO(all_with_hpo, all_without_hpo)

    def _run_analysis(self, predicate, variable1, variable2):
        final_dict = dict()
        if not self.is_recessive:
            for hpo in self._testing_hpo_terms:
                with_hpo_var1_count = len([pat for pat in self._patients_by_hpo.all_with_hpo.get(hpo) if predicate.test(pat, variable1)])
                not_hpo_var1_count = len([pat for pat in self._patients_by_hpo.all_without_hpo.get(hpo) if predicate.test(pat, variable1)])
                if variable2 is None:
                    with_hpo_var2_count = len([pat for pat in self._patients_by_hpo.all_with_hpo.get(hpo) if not predicate.test(pat, variable1)])
                    not_hpo_var2_count = len([pat for pat in self._patients_by_hpo.all_without_hpo.get(hpo) if not predicate.test(pat, variable1)])
                else:
                    with_hpo_var2_count = len([pat for pat in self._patients_by_hpo.all_with_hpo.get(hpo) if predicate.test(pat, variable2)])
                    not_hpo_var2_count = len([pat for pat in self._patients_by_hpo.all_without_hpo.get(hpo) if predicate.test(pat, variable2)])
                if with_hpo_var1_count + not_hpo_var1_count == 0 or with_hpo_var2_count + not_hpo_var2_count == 0:
                    self._logger.warning(f"Error with HPO {hpo.identifier.value}, not included in this analysis.")
                    continue
                p_val = self._run_fisher_exact([[with_hpo_var1_count, not_hpo_var1_count], [with_hpo_var2_count, not_hpo_var2_count]])
                final_dict[f"{hpo.identifier.value} ({hpo.name})"] = [with_hpo_var1_count, (with_hpo_var1_count/(with_hpo_var1_count + not_hpo_var1_count)) * 100, with_hpo_var2_count, (with_hpo_var2_count/(with_hpo_var2_count + not_hpo_var2_count)) * 100, p_val]     
        return final_dict

    def compare_by_variant_type(self, var_type1:VariantEffect, var_type2:VariantEffect = None):
        var_effect_test = VariantEffectPredicate(self.analysis_transcript)
        final_dict = self._run_analysis(var_effect_test, var_type1, var_type2)
        if var_type2 is None:
            col_name = f'without {var_type1.effect_name}'
        else:
            col_name = f'with {var_type2.effect_name}'
        final_df = DataFrame.from_dict(final_dict, orient='index', columns=['with count', f'% with {var_type1.effect_name}', 'count', f'% {col_name}', 'p-value'])
        return final_df.sort_values('p-value')

    def compare_by_variant(self, variant1:str, variant2:str = None):
        variant_test = VariantPredicate(self.analysis_transcript)
        final_dict = self._run_analysis(variant_test, variant1, variant2)
        if variant2 is None:
            col_name = f'without {variant1}'
        else:
            col_name = f'with {variant2}'
        final_df = DataFrame.from_dict(final_dict, orient='index', columns=['with count', f'% with {variant1}', 'count' ,f'% {col_name}', 'p-value'])
        return final_df.sort_values('p-value')

    def compare_by_exon(self, exon1:int, exon2:int = None):
        exon_test = ExonPredicate(self.analysis_transcript)
        final_dict = self._run_analysis(exon_test, exon1, exon2)
        if exon2 is None:
            col_name = f'outside exon {exon1}'
        else:
            col_name = f'inside exon {exon2}'
        final_df = DataFrame.from_dict(final_dict, orient='index', columns=['with count', f'% inside exon {exon1}', 'count', f'% {col_name}', 'p-value'])
        return final_df.sort_values('p-value')

    def compare_by_protein_feature_type(self, feature1:FeatureType, feature2:FeatureType = None):
        feat_type_test = ProtFeatureTypePredicate(self.analysis_transcript)
        final_dict = self._run_analysis(feat_type_test, feature1, feature2)
        if feature2 is None:
            col_name = f'outside exon {feature1}'
        else:
            col_name = f'inside exon {feature2}'
        final_df = DataFrame.from_dict(final_dict, orient='index', columns=['with count', f'% inside {feature1}', 'count', f'% {col_name}', 'p-value'])
        return final_df.sort_values('p-value')


    def compare_by_protein_feature(self, feature1:str, feature2:str = None):
        domain_test = ProtFeaturePredicate(self.analysis_transcript)
        final_dict = self._run_analysis(domain_test, feature1, feature2)
        if feature2 is None:
            col_name = f'outside exon {feature1}'
        else:
            col_name = f'inside exon {feature2}'
        final_df = DataFrame.from_dict(final_dict, orient='index', columns=['with count', f'% inside {feature1}', 'count', f'% {col_name}', 'p-value'])
        return final_df.sort_values('p-value')
    

    def _run_fisher_exact(self, two_by_two_table: typing.Sequence[typing.Sequence[int]]):
        oddsr, p = stats.fisher_exact(two_by_two_table, alternative='two-sided')
        return p



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

