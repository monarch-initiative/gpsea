import abc
import logging
import math

import numpy as np

from ._cohort_data import Cohort
import typing
from genophenocorr.constants import VariantEffect
from genophenocorr.protein import FeatureType, ProteinFeature
from genophenocorr.predicate import VariantEffectPredicate, HPOPresentPredicate, \
    VariantPredicate, ExonPredicate, ProtFeatureTypePredicate, ProtFeaturePredicate, HasVariantResults, PatientCategory
from genophenocorr.variant import Variant
from scipy import stats
from pandas import DataFrame, MultiIndex
from collections import Counter, namedtuple
from enum import Flag


class CohortAnalysis():
    """This class can be used to run different analyses of a given Cohort. 

    Attributes:
        analysis_cohort (Cohort): The Cohort object given for this analysis
        analysis_patients (list): A subset of patients with no structural variants if they are being removed from the analysis. Otherwise, all patients in the cohort.
        analysis_transcript (string): The given transcript ID that will be used in this analysis
        is_recessive (boolean): True if the variant is recessive. Default is False. 
        include_unmeasured (boolean): True if we want to assume a patient without a specific phenotype listed does not have that phenotype. Otherwise, will only include those that explicitely say the phenotype was not observed. Defaults to True. 
        include_large_SV (boolean): True if we want to include large structural variants in the analysis. Defaults to True.
        min_perc_patients_w_hpo (integer): Will only test phenotypes that are listed in at least this percent of patients. Defaults to 10 (10%).
    Methods:
        compare_by_variant_type(var_type1:VariantEffect, var_type2:Optional VariantEffect): Runs Fisher Exact analysis, finds any correlation between given variant effects and phenotypes
        compare_by_variant_type(variant1:string, variant2:Optional string): Runs Fisher Exact analysis, finds any correlation between given variants and phenotypes
        compare_by_variant_type(exon1:integer, exon2:Optional integer): Runs Fisher Exact analysis, finds any correlation between given affected exons and phenotypes
        compare_by_variant_type(feature1:FeatureType, feature2:Optional FeatureType): Runs Fisher Exact analysis, finds any correlation between given protein feature types and phenotypes
        compare_by_variant_type(feature1:string, feature2:Optional string): Runs Fisher Exact analysis, finds any correlation between given protein features and phenotypes
    """
    def __init__(self, cohort, transcript, recessive = False, include_unmeasured = True,  
                include_large_SV = True, min_perc_patients_w_hpo = 10) -> None:
        """Constructs all necessary attributes for a CohortAnalysis object 

        Args:
            cohort (Cohort): The Cohort object that will be analyzed
            transcript (string): The transcript ID that will be used in this analysis
            recessive (boolean): True if the variant is recessive. Default is False. 
            include_unmeasured (boolean): True if we want to assume a patient without a specific phenotype listed does not have that phenotype. Otherwise, will only include those that explicitely say the phenotype was not observed. Defaults to True. 
            include_large_SV (boolean): True if we want to include large structural variants in the analysis. Defaults to True.
            min_perc_patients_w_hpo (integer): Will only test phenotypes that are listed in at least this percent of patients. Defaults to 10 (10%).
        """
        if not isinstance(cohort, Cohort):
            raise ValueError(f"cohort must be type Cohort but was type {type(cohort)}")
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
        """
        Returns:
            Cohort: The Cohort object being analyzed
        """
        return self._cohort
    
    @property
    def analysis_patients(self):
        """
        Returns:
            list: A subset of patients with no structural variants if they are being removed from the analysis. Otherwise, all patients in the cohort.
        """
        return self._patient_list
    
    @property
    def analysis_transcript(self):
        """
        Returns:
            string: The given transcript ID that will be used in this analysis
        """
        return self._transcript

    @property
    def is_recessive(self):
        """
        Returns:
            boolean: True if the variant is recessive.
        """
        return self._recessive

    @property
    def include_unmeasured(self):
        """
        Returns:
            boolean: True if we want to assume a patient without a specific phenotype listed does not have that phenotype. 
            Otherwise, will only include those that explicitely say the phenotype was not observed. 
        """
        return self._in_unmeasured

    @property
    def include_large_SV(self):
        """
        Returns:
            boolean: True if we want to include large structural variants in the analysis.
        """
        return self._include_SV

    @property
    def min_perc_patients_w_hpo(self):
        """
        Returns:
            integer: Will only test phenotypes that are listed in at least this percent of patients. 
        """
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
                with_hpo = [pat for pat in self.analysis_patients if self._hpo_present_test.test(pat, hpo) in PatientCategory.OBSERVED]
                not_hpo = [pat for pat in self.analysis_patients if self._hpo_present_test.test(pat, hpo) in PatientCategory.NOTINCLUDED]
            else:
                with_hpo = [pat for pat in self.analysis_patients if self._hpo_present_test.test(pat, hpo) in PatientCategory.OBSERVED]
                not_hpo = [pat for pat in self.analysis_patients if self._hpo_present_test.test(pat, hpo) in PatientCategory.NOTOBSERVED]
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
        """Runs Fisher Exact analysis, finds any correlation between given variant effects across phenotypes.
        
        Args:
            var_type1 (VariantEffect): 
            var_type2 (VariantEffect, Optional): If None, we compare between those with var_type1 and those without var_type1
        Returns:
            DataFrame: A pandas DataFrame showing the results of the analysis
        """
        var_effect_test = VariantEffectPredicate(self.analysis_transcript)
        if self.is_recessive:
            final_dict = self._run_rec_analysis(var_effect_test, var_type1, var_type2)
            if var_type2 is None:
                col_name1 = f'Heterozygous {str(var_type1)}'
                col_name2 = f'No {str(var_type1)}'
            else:
                col_name1 = f'Heterozygous {str(var_type1)} and {str(var_type2)}'
                col_name2 = f'Homozygous {str(var_type2)}'
            index = MultiIndex.from_product([[f'Homozygous {str(var_type1)}', col_name1, col_name2], ['Count', 'Percent']])
            final_df = DataFrame.from_dict(final_dict, orient='index', columns=index.insert(6, ('', "p-value")))
        else:
            final_dict = self._run_dom_analysis(var_effect_test, var_type1, var_type2)
            if var_type2 is None:
                col_name = f'Without {str(var_type1)}'
            else:
                col_name = f'With {str(var_type2)}'
            index = MultiIndex.from_product([[f'With {str(var_type1)}', col_name], ['Count', 'Percent']])
            final_df = DataFrame.from_dict(final_dict, orient='index', columns=index.insert(4, ('', 'p-value')))
        return final_df.sort_values(('', 'p-value'))

    def compare_by_variant(self, variant1:str, variant2:str = None):
        """Runs Fisher Exact analysis, finds any correlation between given variants across phenotypes.
        
        Args:
            variant1 (string): 
            variant2 (string, Optional): If None, we compare between those with variant1 and those without variant1
        Returns:
            DataFrame: A pandas DataFrame showing the results of the analysis
        """
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
        """Runs Fisher Exact analysis, finds any correlation between given exons across phenotypes.
        
        Args:
            exon1 (integer): 
            exon2 (integer, Optional): If None, we compare between those within exon1 and those outside exon1
        Returns:
            DataFrame: A pandas DataFrame showing the results of the analysis
        """
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
        """Runs Fisher Exact analysis, finds any correlation between given feature type across phenotypes.
        
        Args:
            feature1 (FeatureType): 
            feature2 (FeatureType, Optional): If None, we compare between those with feature1 and those without feature1
        Returns:
            DataFrame: A pandas DataFrame showing the results of the analysis
        """
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
        """Runs Fisher Exact analysis, finds any correlation between given feature and phenotypes.
        
        Args:
            feature1 (string): 
            feature2 (string, Optional): If None, we compare between those within feature1 and those outside feature1
        Returns:
            DataFrame: A pandas DataFrame showing the results of the analysis
        """
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
        a = np.array(two_by_three_table, dtype=np.int64)
        test_class = PythonMultiFisherExact()
        val = test_class.calculate(a)
        return val


class MultiFisherExact(metaclass=abc.ABCMeta):

    @abc.abstractmethod
    def calculate(self, a: np.ndarray) -> float:
        """
        :param a: a 2x3 int array with counts
        :returns: a p value calculated with Fisher's exact test
        """
        pass

    @staticmethod
    def _check_input(a: np.ndarray):
        if not isinstance(a, np.ndarray):
            raise ValueError(f'Expected a numpy array but got {type(a)}')
        if not a.shape == (3, 2):
            raise ValueError(f'Shape of the array must be (3, 2) but got {a.shape}')
        if np.array_equal(a, np.zeros_like(a)):
            raise ValueError(f'Array is all zeros, cannot run analysis')
        


class PythonMultiFisherExact(MultiFisherExact):

    def calculate(self, a: np.ndarray) -> float:
        MultiFisherExact._check_input(a)
        return self._fisher_exact(a)

    def _fisher_exact(self, table):
        row_sum = []
        col_sum = []

        for i in range(len(table)):
            temp = 0
            for j in range(len(table[0])):
                temp += table[i][j]
            row_sum.append(temp)

        for j in range(len(table[0])):
            temp = 0
            for i in range(len(table)):
                temp += table[i][j]
            col_sum.append(temp)

        mat = [[0] * len(col_sum)] * len(row_sum)
        pos = (0, 0)

        p_0 = 1

        for x in row_sum:
            p_0 *= math.factorial(x)
        for y in col_sum:
            p_0 *= math.factorial(y)

        n = 0
        for x in row_sum:
            n += x
        p_0 /= math.factorial(n)

        for i in range(len(table)):
            for j in range(len(table[0])):
                p_0 /= math.factorial(table[i][j])

        p = [0]
        self._dfs(mat, pos, row_sum, col_sum, p_0, p)

        return p[0]

    def _dfs(self, mat, pos, r_sum, c_sum, p_0, p):

        (xx, yy) = pos
        (r, c) = (len(r_sum), len(c_sum))

        mat_new = []

        for i in range(len(mat)):
            temp = []
            for j in range(len(mat[0])):
                temp.append(mat[i][j])
            mat_new.append(temp)

        if xx == -1 and yy == -1:
            for i in range(r - 1):
                temp = r_sum[i]
                for j in range(c - 1):
                    temp -= mat_new[i][j]
                mat_new[i][c - 1] = temp
            for j in range(c - 1):
                temp = c_sum[j]
                for i in range(r - 1):
                    temp -= mat_new[i][j]
                mat_new[r - 1][j] = temp
            temp = r_sum[r - 1]
            for j in range(c - 1):
                temp -= mat_new[r - 1][j]
            if temp < 0:
                return
            mat_new[r - 1][c - 1] = temp

            p_1 = 1
            for x in r_sum:
                p_1 *= math.factorial(x)
            for y in c_sum:
                p_1 *= math.factorial(y)

            n = 0
            for x in r_sum:
                n += x
            p_1 /= math.factorial(n)

            for i in range(len(mat_new)):
                for j in range(len(mat_new[0])):
                    p_1 /= math.factorial(mat_new[i][j])
            if p_1 <= p_0 + 0.00000001:
                # print(mat_new)
                # print(p_1)
                p[0] += p_1
        else:
            max_1 = r_sum[xx]
            max_2 = c_sum[yy]
            for j in range(c):
                max_1 -= mat_new[xx][j]
            for i in range(r):
                max_2 -= mat_new[i][yy]
            for k in range(min(max_1, max_2) + 1):
                mat_new[xx][yy] = k
                if xx == r - 2 and yy == c - 2:
                    pos_new = (-1, -1)
                elif xx == r - 2:
                    pos_new = (0, yy + 1)
                else:
                    pos_new = (xx + 1, yy)
                self._dfs(mat_new, pos_new, r_sum, c_sum, p_0, p)
