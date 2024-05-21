import typing

from hpotk import MinimalOntology
from jinja2 import Environment, PackageLoader

from genophenocorr.model import Cohort, Variant


class CohortViewable:
    """
    Class to create a viewable object that is uses a Jinja2 template to create an HTML element
    for display in the Jupyter notebook.
    """

    def __init__(
            self,
            hpo: MinimalOntology,
            top_phenotype_count: int = 10,
            top_variant_count: int = 10,
            transcript_id: str = None,
    ):
        """
        Args:
            hpo(MinimalOntology): An HPO ontology object from hpo-toolkit
            top_phenotype_count(int): Maximum number of HPO terms to display in the HTML table (default: 10)
            top_variant_count(int): Maximum number of variants to display in the HTML table (default: 10)
            transcript_id(str): The identifier (Only NCBI RefSeq supported!) to use for the HGVS-able variants
        """
        self._hpo = hpo
        self._top_phenotype_count = top_phenotype_count
        self._top_variant_count = top_variant_count
        self._tx_id = transcript_id

        environment = Environment(loader=(PackageLoader('genophenocorr.view', 'templates')))
        self._cohort_template = environment.get_template("cohort.html")

    def process(
            self,
            cohort: Cohort,
            transcript_id: str
    ) -> str:
        """
        Create an HTML that should be shown with display(HTML(..)) of the ipython package.

        Args:
            cohort (Cohort): The cohort being analyzed in the current Notebook
            transcript_id (str): the transcript that we map variants onto

        Returns:
            str: an HTML string with parameterized template for rendering
        """
        self._transcript_id = transcript_id
        self._variant_to_display_d = self.get_variant_description(cohort, transcript_id)
        context = self._prepare_context(cohort)
        return self._cohort_template.render(context)

    def _prepare_context(
            self,
            cohort: Cohort,
    ) -> typing.Mapping[str, typing.Any]:
        n_patients = len(cohort.all_patients)
        total_HPO_count = len(cohort.all_phenotypes())
        top_hpos = cohort.list_present_phenotypes(top=self._top_phenotype_count)
        hpo_counts = list()
        for hpo in top_hpos:
            hpo_id = hpo[0]
            individual_count = hpo[1]
            hpo_label = "n/a"
            if hpo_id in self._hpo:
                hpo_label = self._hpo.get_term(hpo_id).name
            hpo_counts.append({"HPO": hpo_label, "ID": hpo_id, "Count": individual_count})
        top_vars = cohort.list_all_variants(top=self._top_variant_count)
        var_counts = list()
        all_variant_set = cohort.all_variants()
        for var in top_vars:
            chrom_var = var[0]
            var_count = var[1]
            # get HGVS or human readable variant
            var_name = self._variant_to_display_d.get(chrom_var, chrom_var)
            # TODO get HGVS or human readable variant
            if self._tx_id is not None:
                variant = cohort.get_variant_by_key(chrom_var)
                var_name = variant.get_hgvs_cdna_by_tx(self._tx_id)
            else:
                var_name = "ID with transcript"
            var_counts.append({"variant": chrom_var, "variant_name": var_name, "Count": var_count})
        diseases = cohort.list_all_diseases()
        n_diseases = len(diseases)
        disease_counts = list()
        for d in diseases:
            disease_id = d[0]
            disease_count = d[1]
            disease_name = "Unknown"
            for dis in cohort.all_diseases():
                if dis.identifier == d[0]:
                    disease_name = dis.name
            disease_counts.append({"disease_id": disease_id, "disease_name": disease_name, "count": disease_count})
        var_effects_list = list()
        if self._tx_id is not None:
            has_transcript = 1
            data_by_tx = cohort.variant_effect_count_by_tx(tx_id=self._tx_id)
            counter_d = data_by_tx.get(self._tx_id)
            # data structure -- {'effect}': 'FRAMESHIFT_VARIANT', 'count': 175}, {'effect}': 'STOP_GAINED', 'count': 67},
            for k, v in counter_d.items():
                var_effects_list.append({"effect": k, "count": v})
        else:
            has_transcript = 0
        # The following dictionary is used by the Jinja2 HTML template
        return {
            "n_individuals": n_patients,
            "n_excluded": cohort.get_excluded_count(),
            "total_hpo_count": total_HPO_count,
            "top_hpo_count": self._top_phenotype_count,
            "hpo_counts": hpo_counts,
            "unique_variant_count": len(cohort.all_variants()),
            "top_var_count": self._top_variant_count,
            "var_counts": var_counts,
            "n_diseases": n_diseases,
            "disease_counts": disease_counts,
            "has_transcript": has_transcript,
            "var_effects_list": var_effects_list,
            "transcript_id": self._transcript_id,
        }
    
    @staticmethod
    def get_display(variant:Variant, transcript_id:str) -> str:
        for annot in variant.tx_annotations:
            if annot.transcript_id == transcript_id:
                return annot.hgvs_cdna
        return variant.variant_string
    
    def get_variant_description(self, cohort: Cohort,
            transcript_id: str,
            onlyHgvs:bool=True) -> typing.Dict[str, str]:
        """
        Get user-friendly strings (e.g., HGV for our target trancript) to match to the chromosomal strings
        Args:
            cohort (Cohort): The cohort being analyzed in the current Notebook
            transcript_id (str): the transcript that we map variants onto

        Returns:
            typing.Dict[str, str]: key: chromosomal, value: display (e.g. HGVS) string of variant
        """
        chrom_to_display = dict()
        all_var_set = cohort.all_variants()
        for var in all_var_set:
            var_string = var.variant_string
            display = CohortViewable.get_display(variant=var, transcript_id=transcript_id)
            if onlyHgvs:
                # do not show the transcript id
                fields = display.split(":")
                if len(fields)>1:
                    display = fields[1]
                else:
                    display = [0]
            chrom_to_display[var_string] = display
        return chrom_to_display
