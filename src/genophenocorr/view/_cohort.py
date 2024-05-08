from hpotk import MinimalOntology
from collections import defaultdict
from genophenocorr.model import Cohort
from genophenocorr.model import Patient, Cohort
import os
from jinja2 import Environment, FileSystemLoader



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

class CohortViewable:
    """
    Class to create a viewable object that is uses a Jinja2 template to create an HTML element for display in the Jupyter notebook.
    """
    def __init__(self, cohort:Cohort, 
                 hpo_ontology:MinimalOntology, 
                 top_phenotype_count:int=10,
                 top_variant_count:int=10,
                 transcript_id:str=None) -> None:
        """
        Args:
            cohort (Cohort): The cohort being analyzed in the current Notebook
            hpo_ontology(MinimalOntology): An HPO ontology object from hpo-toolkit
            top_phenotype_count(int): Maximum number of HPO terms to display in the HTML table (default: 10)
            top_variant_count(int): Maximum number of variants to display in the HTML table (default: 10)
            transcript_id(str): The identifier (Only NCBI RefSeq supported!) to use for the HGVS-able variants
        """
        self._cohort = Cohort
        n_patients = len(cohort.all_patients)
        total_HPO_count = len(cohort.all_phenotypes)
        top_hpos = cohort.list_all_phenotypes(top=top_phenotype_count)
        hpo_counts = list()
        for hpo in top_hpos:
            hpo_id = hpo[0]
            individual_count = hpo[1]
            hpo_label = "n/a"
            if hpo_id in hpo_ontology:
                hpo_label = hpo_ontology.get_term_name(hpo_id)
            hpo_counts.append({"HPO": hpo_label, "ID": hpo_id, "Count": individual_count})
        top_vars = cohort.list_all_variants(top=top_variant_count)
        var_counts = list()
        for var in top_vars:
            chrom_var = var[0]
            var_count = var[1]
            # TODO get HGVS or human readable variant
            var_name = "todo"
            var_counts.append({"variant": chrom_var, "variant_name":var_name, "Count": var_count})
        diseases = cohort.list_all_diseases()
        n_diseases = len(diseases)
        disease_counts = list()
        for d in diseases:
            disease_id = d[0]
            disease_count = d[1]
            disease_counts.append({"disease_id": disease_id, "count": disease_count})
        var_effects_list = list()
        if transcript_id is not None:
            has_transcript = 1
            data_by_tx = cohort.list_data_by_tx(transcript=transcript_id)
            counter_d = data_by_tx.get(transcript_id)
            # data structure -- {'effect}': 'FRAMESHIFT_VARIANT', 'count': 175}, {'effect}': 'STOP_GAINED', 'count': 67}, 
            for k,v in counter_d.items():
                var_effects_list.append({"effect": k, "count": v})
        else:
            has_transcript = 0
        # The following dictionary is used by the Jinja2 HTML template
        self._context = {
            "n_individuals": n_patients,
            "n_excluded": cohort.get_excluded_count(),
            "total_hpo_count": total_HPO_count,
            "top_hpo_count": top_phenotype_count,
            "hpo_counts": hpo_counts,
            "unique_variant_count": len(cohort._variant_set),
            "top_var_count": top_variant_count,
            "var_counts": var_counts,
            "n_diseases": n_diseases,
            "disease_counts":disease_counts, 
            "has_transcript": has_transcript,
            "var_effects_list": var_effects_list,
        }
        


    def get_html(self) -> str:
        """
        This method uses a Jinja2 HTML template to create HTML that should be shown with display(HTML(..)) of the ipython package
        Returns:
            str: An HTML string with paragraphs and tables containing a summary of the cohort.
        """
        parent_dir = os.path.dirname(os.path.abspath(__file__))
        template_dir = os.path.join(parent_dir, "templates/")
        environment = Environment(loader=FileSystemLoader(template_dir))
        template = environment.get_template("cohort.html")
        return template.render(self._context)

