from hpotk import MinimalOntology
from collections import defaultdict
from genophenocorr.model import Cohort



CSS_CODE = """
table {
  border-collapse: collapse; 
  margin: 25px 0;
    font-size: 0.9em;
    font-family: sans-serif;
    min-width: 400px;
    box-shadow: 0 0 20px rgba(0, 0, 0, 0.15);
}


th {
  background-color: #f2f2f2; 
  border: 1px solid #dddddd; 
  text-align: left;
  padding: 2px;
}

tr {
  border: 1px solid #dddddd; 
}

td {
  padding: 2px; 
  font-weight: bold;
}

tr:nth-child(even) {
  background-color: #f2f2f2;
}
"""



class CohortViewer:
    def __init__(self, hpo):
        if not isinstance(hpo, MinimalOntology):
            raise ValueError(f"hpo argument must be a MinimalOntology but we got a {type(hpo)} object")
        self._hpo = hpo

    @staticmethod
    def html_row(items):
        line = "</td><td>".join(items)
        return f"<tr><td>{line}</td></tr>"

    @staticmethod
    def anchor(url, payload):
        return f"<a href=\"{url}\">{payload}</a>"

    @staticmethod
    def hpo_anchor(hpo_id, label):
        url = f"https://hpo.jax.org/app/browse/term/{hpo_id}"
        hpo_label = f"{label} ({hpo_id})"
        return CohortViewer.anchor(url, hpo_label)

    def hpo_term_counts_table(self, cohort, min_count=1) -> str:
        """
        Generate HTML code designed to be displayed on a Jupyter notebook using ipython/display/HTML
        Show the number of annotations per HPO terms. Provide an explanation.
        
        :param cohort: A cohort of patients to be analyzed
        :type cohort: Cohort
        :param min_count: Minimum number of annotations to be displayed in the table
        :type min_count: int
        :returns: HTML code for display
        """
        if not isinstance(cohort, Cohort):
            raise ValueError(f"cohort argument must be a Cohort object but we got a {type(cohort)} object")
        rows = list()
        rows.append(f"<style>\n{CSS_CODE}</style>\n")
        rows.append("<table>")
        header_items = ["HPO Term", "Count"]
        rows.append(CohortViewer.html_row(header_items))
        hpo_list = cohort.list_all_phenotypes()
        n_individuals = cohort.total_patient_count
        capt = f"<caption>Counts of annotations to HPO terms for the {n_individuals} in the cohort</caption>"
        rows.append(capt)
        min_terms = set()
        for hpo_tuple in hpo_list:
            row_items = []
            hpo_id = hpo_tuple[0]
            count = hpo_tuple[1]
            term = self._hpo.get_term(hpo_id)
            if term is None:
                hpo_label = "n/a"
            else:
                hpo_label = term.name
            if count >= min_count:
                row_items.append(CohortViewer.hpo_anchor(hpo_id, hpo_label))
                row_items.append(str(count))
                rows.append(CohortViewer.html_row(row_items))
            else:
                min_terms.add(CohortViewer.hpo_anchor(hpo_id, hpo_label))
        rows.append("</table>")
        if len(min_terms) > 0:
            term_str = "; ".join(min_terms)
            rows.append(f"<p>Additionally, the following terms were observed {min_count-1} or fewer times:</p>")
            rows.append(f"<p>{term_str}.</p>")
        return "\n".join(rows)


    def variants_table(self, cohort, preferred_transcript, min_count=2) -> str:
        """
        Generate HTML code designed to be displayed on a Jupyter notebook using ipython/display/HTML
        Show genotype-phenotype correlation tests that could be run.
        
        :param cohort: A cohort of patients to be analyzed
        :type cohort: Cohort
        :param min_count: Minimum number of annotations to be displayed in the table
        :type min_count: int
        :returns: HTML code for display
        """
        if not isinstance(cohort, Cohort):
            raise ValueError(f"cohort argument must be a Cohort object but we got a {type(cohort)} object")
        rows = list()
        # Get a list of variants for the preferred transcript
        variant_count_d = defaultdict(int)
        variant_to_effect_d = defaultdict()
        variant_to_key = defaultdict()
        # key, variant string, value int
        all_variant_tuple_list = cohort.list_all_variants()
        if not isinstance(all_variant_tuple_list, list):
            raise ValueError(f"all_variant_tuple_list is not a list but is a {type(all_variant_tuple_list)}")
        all_variant_counter = {x[0]:x[1] for x in all_variant_tuple_list}
        for variant in cohort.all_variants:
            var_count = all_variant_counter[variant.variant_string]
            targets = [txa for txa in variant.tx_annotations if txa.transcript_id == preferred_transcript]
            if len(targets) == 1:
                target_txa = targets[0]
                if target_txa.hgvs_cdna is not None:
                    hgvs_cdna = target_txa.hgvs_cdna
                else:
                    hgvs_cdna = "NA"
                # split out the variant
                fields = hgvs_cdna.split(":")
                if len(fields) == 2:
                    hgvs = fields[1]
                else:
                    hgvs = hgvs_cdna
                effect_tuple = [var_eff.name for var_eff in target_txa.variant_effects]
                variant_count_d[hgvs] = var_count
                variant_to_effect_d[hgvs] = effect_tuple[0] # for simplicity, just display first effect
                variant_to_key[hgvs] = variant.variant_string
            else:
                print(f"[WARN] could not identify a single variant for target transcript (got {len(targets)}), variant {variant.variant_string}")
                # could not find an entry for our transcript, so just show the genomic coordinates
                variant_count_d[variant.variant_string] = var_count
                variant_to_key[variant.variant_string] = variant.variant_string
        # sort the variants by count and put variants below mininum count for display into a separate set
        sorted_vars = sorted(variant_count_d.items(), key=lambda x:x[1], reverse=True) # sort descending by values
        sorted_vars = [x[0] for x in sorted_vars] # take the first item from each sorted tuple
        below_threshold_vars = set()
        rows.append(f"<style>\n{CSS_CODE}</style>\n")
        rows.append("<table>")
        header_items = ["Variant", "Effect", "Count", "Key"]
        rows.append(CohortViewer.html_row(header_items))
        for var in sorted_vars:
            items = []
            var_count = variant_count_d.get(var)
            #print(f"{var} - {var_count}")
            if var_count >= min_count: 
                variant_key = variant_to_key.get(var)
                items.append(var)
                items.append(variant_to_effect_d.get(var, "n/a"))
                items.append(str(var_count))
                items.append(variant_key)
                rows.append(CohortViewer.html_row(items))
            else:
                below_threshold_vars.add(var)
        rows.append("</table>")
        if len(below_threshold_vars) > 0:
            var_str = "; ".join(below_threshold_vars)
            rows.append(f"<p>Additionally, the following variants were observed {min_count-1} or fewer times: ")
            rows.append(f"{var_str}.</p>")
        rows.append("<p>Use the entry in the \"Key\" column to investigate whether specific variants display genotype-phenotype correlations</p>")
        return "\n".join(rows)


    def cohort_summary_table(self, cohort, min_count=1) -> str:
        """
        Generate HTML code designed to be displayed on a Jupyter notebook using ipython/display/HTML
        Show the number of annotations per HPO terms. Provide an explanation.
        
        :param cohort: A cohort of patients to be analyzed
        :type cohort: Cohort
        :param min_count: Minimum number of annotations to be displayed in the table
        :type min_count: int
        :returns: HTML code for display
        """
        if not isinstance(cohort, Cohort):
            raise ValueError(f"cohort argument must be a Cohort object but we got a {type(cohort)} object")
        rows = list()
        rows.append(f"<style>\n{CSS_CODE}</style>\n")
        rows.append("<table>")
        header_items = ["Item", "Description"]
        rows.append(CohortViewer.html_row(header_items))
        n_individuals = cohort.total_patient_count
        n_unique_hpo = len(cohort.all_phenotypes)
        n_unique_variants = len(cohort.all_variants)
        n_excluded_patients = cohort.get_excluded_count()
        excluded_patients = cohort.get_excluded_ids()

        cap = "Description of the cohort. "
        if n_excluded_patients > 0:
            cap =  cap + f"{n_excluded_patients} individuals were removed from the cohort because they had no HPO terms."

        capt = f"<caption>{cap}</caption>"
        rows.append(capt)
        
        rows.append(CohortViewer.html_row(["Total Individuals", str(n_individuals)]))
        #TODO: Add Diseases
        if n_excluded_patients > 0:
            rows.append(CohortViewer.html_row(["Excluded Individuals", f"{str(n_excluded_patients)}: {';'.join(excluded_patients)}"]))
        rows.append(CohortViewer.html_row(["Total Unique HPO Terms", str(n_unique_hpo)]))
        rows.append(CohortViewer.html_row(["Total Unique Variants", str(n_unique_variants)]))

        rows.append("</table>")
    
        return "\n".join(rows)


    def protein_features_table(self, cohort, preferred_transcript, min_count=2) -> str:
        """
        Generate HTML code designed to be displayed on a Jupyter notebook using ipython/display/HTML
        Show genotype-phenotype correlation tests that could be run.
        
        :param cohort: A cohort of patients to be analyzed
        :type cohort: Cohort
        :param min_count: Minimum number of annotations to be displayed in the table
        :type min_count: int
        :returns: HTML code for display
        """
        if not isinstance(cohort, Cohort):
            raise ValueError(f"cohort argument must be a Cohort object but we got a {type(cohort)} object")
        rows = list()
        # Get a list of variants for the preferred transcript
        variant_count_d = defaultdict(int)
        variant_to_effect_d = defaultdict()
        variant_to_key = defaultdict()
        # key, variant string, value int
        all_variant_tuple_list = cohort.list_all_variants()
        if not isinstance(all_variant_tuple_list, list):
            raise ValueError(f"all_variant_tuple_list is not a list but is a {type(all_variant_tuple_list)}")
        all_variant_counter = {x[0]:x[1] for x in all_variant_tuple_list}
        for variant in cohort.all_variants:
            var_count = all_variant_counter[variant.variant_string]
            targets = [txa for txa in variant.tx_annotations if txa.transcript_id == preferred_transcript]
            if len(targets) == 1:
                target_txa = targets[0]
                if target_txa.hgvs_cdna is not None:
                    hgvs_cdna = target_txa.hgvs_cdna
                else:
                    hgvs_cdna = "NA"
                # split out the variant
                fields = hgvs_cdna.split(":")
                if len(fields) == 2:
                    hgvs = fields[1]
                else:
                    hgvs = hgvs_cdna
                effect_tuple = [var_eff.name for var_eff in target_txa.variant_effects]
                variant_count_d[hgvs] = var_count
                variant_to_effect_d[hgvs] = effect_tuple[0] # for simplicity, just display first effect
                variant_to_key[hgvs] = variant.variant_string
            else:
                print(f"[WARN] could not identify a single variant for target transcript (got {len(targets)}), variant {variant.variant_string}")
                # could not find an entry for our transcript, so just show the genomic coordinates
                variant_count_d[variant.variant_string] = var_count
                variant_to_key[variant.variant_string] = variant.variant_string
        # sort the variants by count and put variants below mininum count for display into a separate set
        sorted_vars = sorted(variant_count_d.items(), key=lambda x:x[1], reverse=True) # sort descending by values
        sorted_vars = [x[0] for x in sorted_vars] # take the first item from each sorted tuple
        below_threshold_vars = set()
        rows.append(f"<style>\n{CSS_CODE}</style>\n")
        rows.append("<table>")
        header_items = ["Variant", "Effect", "Count", "Key"]
        rows.append(CohortViewer.html_row(header_items))
        for var in sorted_vars:
            items = []
            var_count = variant_count_d.get(var)
            #print(f"{var} - {var_count}")
            if var_count >= min_count: 
                variant_key = variant_to_key.get(var)
                items.append(var)
                items.append(variant_to_effect_d.get(var, "n/a"))
                items.append(str(var_count))
                items.append(variant_key)
                rows.append(CohortViewer.html_row(items))
            else:
                below_threshold_vars.add(var)
        rows.append("</table>")
        if len(below_threshold_vars) > 0:
            var_str = "; ".join(below_threshold_vars)
            rows.append(f"<p>Additionally, the following variants were observed {min_count-1} or fewer times: ")
            rows.append(f"{var_str}.</p>")
        rows.append("<p>Use the entry in the \"Key\" column to investigate whether specific variants display genotype-phenotype correlations</p>")
        return "\n".join(rows)