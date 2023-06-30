import logging
from ._cohort_data import Cohort
import typing
from genophenocorr.constants import VariantEffect
from genophenocorr.protein import FeatureType, ProteinFeature
from genophenocorr.predicate import VariantEffectPredicate, HPOPresentPredicate, \
    VariantPredicate, ExonPredicate, ProtFeatureTypePredicate, ProtFeaturePredicate, HasVariantResults
from genophenocorr.variant import Variant
from scipy import stats
from pandas import DataFrame, MultiIndex
from collections import Counter, namedtuple
from enum import Flag
from FisherExact import fisher_exact


class CohortAnalysis():

    def __init__(self, cohort, transcript, recessive = False, include_unmeasured = True,  
                include_large_SV = True, min_perc_patients_w_hpo = 10) -> None:
        if not isinstance(cohort, Cohort):
            raise ValueError(f"cohort must be type Cohort but was type {type(cohort)}")
        ## IF TRANSCRIPT DOES NOT EXIST, GIVE ERROR
        if transcript not in cohort.all_transcripts:
            raise ValueError(f"Transcript {transcript} not found in Cohort")
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

    def _run_dom_analysis(self, predicate, variable1, variable2):
        final_dict = dict()
        if not self.is_recessive:
            for hpo in self._testing_hpo_terms:
                with_hpo_var1_count = len([pat for pat in self._patients_by_hpo.all_with_hpo.get(hpo) if predicate.test(pat, variable1) in HasVariantResults.DOMINANTVARIANTS])
                not_hpo_var1_count = len([pat for pat in self._patients_by_hpo.all_without_hpo.get(hpo) if predicate.test(pat, variable1) in HasVariantResults.DOMINANTVARIANTS])
                if variable2 is None:
                    with_hpo_var2_count = len([pat for pat in self._patients_by_hpo.all_with_hpo.get(hpo) if predicate.test(pat, variable1) == HasVariantResults.NOVARIANT])
                    not_hpo_var2_count = len([pat for pat in self._patients_by_hpo.all_without_hpo.get(hpo) if predicate.test(pat, variable1) == HasVariantResults.NOVARIANT])
                else:
                    with_hpo_var2_count = len([pat for pat in self._patients_by_hpo.all_with_hpo.get(hpo) if predicate.test(pat, variable2) in HasVariantResults.DOMINANTVARIANTS])
                    not_hpo_var2_count = len([pat for pat in self._patients_by_hpo.all_without_hpo.get(hpo) if predicate.test(pat, variable2) in HasVariantResults.DOMINANTVARIANTS])
                if with_hpo_var1_count + not_hpo_var1_count == 0 or with_hpo_var2_count + not_hpo_var2_count == 0:
                    self._logger.warning(f"Divide by 0 error with HPO {hpo.identifier.value}, not included in this analysis.")
                    continue
                p_val = self._run_fisher_exact([[with_hpo_var1_count, not_hpo_var1_count], [with_hpo_var2_count, not_hpo_var2_count]])
                final_dict[f"{hpo.identifier.value} ({hpo.name})"] = [with_hpo_var1_count, "{:0.2f}%".format((with_hpo_var1_count/(with_hpo_var1_count + not_hpo_var1_count)) * 100), with_hpo_var2_count, "{:0.2f}%".format((with_hpo_var2_count/(with_hpo_var2_count + not_hpo_var2_count)) * 100), p_val]
        else:
            return ValueError(f"Run a recessive analysis")
        return final_dict

    def _run_rec_analysis(self, predicate, variable1, variable2):
        final_dict = dict()
        if self.is_recessive:
            for hpo in self._testing_hpo_terms:
                if variable2 is None:
                    with_hpo_var1_var1 = len([pat for pat in self._patients_by_hpo.all_with_hpo.get(hpo) if predicate.test(pat, variable1) == HasVariantResults.HOMOVARIANT])
                    no_hpo_var1_var1 = len([pat for pat in self._patients_by_hpo.all_without_hpo.get(hpo) if predicate.test(pat, variable1) == HasVariantResults.HOMOVARIANT])
                    with_hpo_var1_var2 = len([pat for pat in self._patients_by_hpo.all_with_hpo.get(hpo) if predicate.test(pat, variable1) == HasVariantResults.HETEROVARIANT])
                    no_hpo_var1_var2 = len([pat for pat in self._patients_by_hpo.all_without_hpo.get(hpo) if predicate.test(pat, variable1) == HasVariantResults.HETEROVARIANT])
                    with_hpo_var2_var2 = len([pat for pat in self._patients_by_hpo.all_with_hpo.get(hpo) if predicate.test(pat, variable1) == HasVariantResults.NOVARIANT])
                    no_hpo_var2_var2 = len([pat for pat in self._patients_by_hpo.all_without_hpo.get(hpo) if predicate.test(pat, variable1) == HasVariantResults.NOVARIANT])
                else:
                    with_hpo_var1_var1 = len([pat for pat in self._patients_by_hpo.all_with_hpo.get(hpo) if predicate.test(pat, variable1) == HasVariantResults.HOMOVARIANT])
                    no_hpo_var1_var1 = len([pat for pat in self._patients_by_hpo.all_without_hpo.get(hpo) if predicate.test(pat, variable1) == HasVariantResults.HOMOVARIANT])
                    with_hpo_var1_var2 = len([pat for pat in self._patients_by_hpo.all_with_hpo.get(hpo) if predicate.test(pat, variable1) == HasVariantResults.HETEROVARIANT and predicate.test(pat, variable2) == HasVariantResults.HETEROVARIANT])
                    no_hpo_var1_var2 = len([pat for pat in self._patients_by_hpo.all_without_hpo.get(hpo) if predicate.test(pat, variable1) == HasVariantResults.HETEROVARIANT and predicate.test(pat, variable2) == HasVariantResults.HETEROVARIANT])
                    with_hpo_var2_var2 = len([pat for pat in self._patients_by_hpo.all_with_hpo.get(hpo) if predicate.test(pat, variable2) == HasVariantResults.HOMOVARIANT])
                    no_hpo_var2_var2 = len([pat for pat in self._patients_by_hpo.all_without_hpo.get(hpo) if predicate.test(pat, variable2) == HasVariantResults.HOMOVARIANT])
                if with_hpo_var1_var1 + no_hpo_var1_var1 == 0 or with_hpo_var1_var2 + no_hpo_var1_var2 == 0 or with_hpo_var2_var2 + no_hpo_var2_var2 == 0:
                    self._logger.warning(f"Divide by 0 error with HPO {hpo.identifier.value}, not included in this analysis.")
                    continue
                p_val = self._run_recessive_fisher_exact([[with_hpo_var1_var1, no_hpo_var1_var1], [with_hpo_var1_var2, no_hpo_var1_var2], [with_hpo_var2_var2, no_hpo_var2_var2]] )
                final_dict[f"{hpo.identifier.value} ({hpo.name})"] = [with_hpo_var1_var1, "{:0.2f}%".format((with_hpo_var1_var1/(with_hpo_var1_var1 + no_hpo_var1_var1)) * 100), with_hpo_var1_var2, "{:0.2f}%".format((with_hpo_var1_var2/(no_hpo_var1_var2 + with_hpo_var1_var2)) * 100), with_hpo_var2_var2, "{:0.2f}%".format((with_hpo_var2_var2/(with_hpo_var2_var2 + no_hpo_var2_var2)) * 100), p_val] 
        else:
            return ValueError("Run a dominant analysis")
        return final_dict


    def compare_by_variant_type(self, var_type1:VariantEffect, var_type2:VariantEffect = None):
        var_effect_test = VariantEffectPredicate(self.analysis_transcript)
        if self.is_recessive:
            final_dict = self._run_rec_analysis(var_effect_test, var_type1, var_type2)
            if var_type2 is None:
                col_name1 = f'Heterozygous {var_type1.effect_name}'
                col_name2 = f'No {var_type1.effect_name}'
            else:
                col_name1 = f'Heterozygous {var_type1.effect_name} and {var_type2.effect_name}'
                col_name2 = f'Homozygous {var_type2.effect_name}'
            index = MultiIndex.from_product([[f'Homozygous {var_type1.effect_name}', col_name1, col_name2], ['Count', 'Percent']])
            final_df = DataFrame.from_dict(final_dict, orient='index', columns=index.insert(6, ('', "p-value")))
        else:
            final_dict = self._run_dom_analysis(var_effect_test, var_type1, var_type2)
            if var_type2 is None:
                col_name = f'Without {var_type1.effect_name}'
            else:
                col_name = f'With {var_type2.effect_name}'
            index = MultiIndex.from_product([[f'With {var_type1.effect_name}', col_name], ['Count', 'Percent']])
            final_df = DataFrame.from_dict(final_dict, orient='index', columns=index.insert(4, ('', 'p-value')))
        return final_df.sort_values(('', 'p-value'))

    def compare_by_variant(self, variant1:str, variant2:str = None):
        variant_test = VariantPredicate(self.analysis_transcript)
        if self.is_recessive:
            final_dict = self._run_rec_analysis(variant_test, variant1, variant2)
            if variant2 is None:
                col_name1 = f'Heterozygous {variant1}'
                col_name2 = f'No {variant1}'
            else:
                col_name1 = f'Heterozygous {variant1} and {variant2}'
                col_name2 = f'Homozygous {variant2}'
            index = MultiIndex.from_product([[f'Homozygous {variant1}', col_name1, col_name2], ['Count', 'Percent']])
            final_df = DataFrame.from_dict(final_dict, orient='index', columns=index.insert(6, ('', 'p-value')))
        else:
            final_dict = self._run_dom_analysis(variant_test, variant1, variant2)  
            if variant2 is None:
                col_name = f'Without {variant1}'
            else:
                col_name = f'With {variant2}'
            index = MultiIndex.from_product([[f'With {variant1}', col_name], ['Count', 'Percent']])
            final_df = DataFrame.from_dict(final_dict, orient='index', columns=index.insert(4, ('', 'p-value')))
        return final_df.sort_values(('', 'p-value'))

    def compare_by_exon(self, exon1:int, exon2:int = None):
        exon_test = ExonPredicate(self.analysis_transcript)
        if self.is_recessive:
            final_dict = self._run_rec_analysis(exon_test, exon1, exon2)
            if exon2 is None:
                col_name1 = f'Heterozygous {exon1}'
                col_name2 = f'No {exon1}'
            else:
                col_name1 = f'Heterozygous {exon1} and {exon2}'
                col_name2 = f'Homozygous {exon2}'
            index = MultiIndex.from_product([[f'Homozygous {exon1}', col_name1, col_name2], ['Count', 'Percent']])
            final_df = DataFrame.from_dict(final_dict, orient='index', columns=index.insert(6, ('', 'p-value')))
        else: 
            final_dict = self._run_dom_analysis(exon_test, exon1, exon2)
            if exon2 is None:
                col_name = f'Outside Exon {exon1}'
            else:
                col_name = f'Inside Exon {exon2}'
            index = MultiIndex.from_product([[f'Inside Exon {exon1}', col_name], ['Count', 'Percent']])
            final_df = DataFrame.from_dict(final_dict, orient='index', columns=index.insert(4, ('', 'p-value')))
        return final_df.sort_values(('', 'p-value'))

    def compare_by_protein_feature_type(self, feature1:FeatureType, feature2:FeatureType = None):
        feat_type_test = ProtFeatureTypePredicate(self.analysis_transcript)
        if self.is_recessive:
            final_dict = self._run_rec_analysis(feat_type_test, feature1, feature2)
            if feature2 is None:
                col_name1 = f'Heterozygous {feature1}'
                col_name2 = f'No {feature1}'
            else:
                col_name1 = f'Heterozygous {feature1} and {feature2}'
                col_name2 = f'Homozygous {feature2}'
            index = MultiIndex.from_product([[f'Homozygous {feature1}', col_name1, col_name2], ['Count', 'Percent']])
            final_df = DataFrame.from_dict(final_dict, orient='index', columns=index.insert(6, ('', 'p-value')))
        else:
            final_dict = self._run_dom_analysis(feat_type_test, feature1, feature2)
            if feature2 is None:
                col_name = f'Outside {feature1.name}'
            else:
                col_name = f'Inside {feature2.name}'
            index = MultiIndex.from_product([[f'Inside {feature1.name}', col_name], ['Count', 'Percent']])
            final_df = DataFrame.from_dict(final_dict, orient='index', columns=index.insert(4, ('', 'p-value')))
        return final_df.sort_values(('', 'p-value'))


    def compare_by_protein_feature(self, feature1:str, feature2:str = None):
        domain_test = ProtFeaturePredicate(self.analysis_transcript)
        if self.is_recessive:
            final_dict = self._run_rec_analysis(domain_test, feature1, feature2)
            if feature2 is None:
                col_name1 = f'Heterozygous {feature1}'
                col_name2 = f'No {feature1}'
            else:
                col_name1 = f'Heterozygous {feature1} and {feature2}'
                col_name2 = f'Homozygous {feature2}'
            index = MultiIndex.from_product([[f'Homozygous {feature1}', col_name1, col_name2], ['Count', 'Percent']])
            final_df = DataFrame.from_dict(final_dict, orient='index', columns=index.insert(6, ('', 'p-value')))
        else:
            final_dict = self._run_dom_analysis(domain_test, feature1, feature2)
            if feature2 is None:
                col_name = f'Outside {feature1}'
            else:
                col_name = f'Inside {feature2}'
            index = MultiIndex.from_product([[f'Inside {feature1}', col_name], ['Count', 'Percent']])
            final_df = DataFrame.from_dict(final_dict, orient='index', columns=index.insert(4, ('', 'p-value')))
        return final_df.sort_values(('', 'p-value'))
    

    def _run_fisher_exact(self, two_by_two_table: typing.Sequence[typing.Sequence[int]]):
        oddsr, p = stats.fisher_exact(two_by_two_table, alternative='two-sided')
        return p

    def _run_recessive_fisher_exact(self, two_by_three_table: typing.Sequence[typing.Sequence[int]]):
        return fisher_exact(two_by_three_table)
