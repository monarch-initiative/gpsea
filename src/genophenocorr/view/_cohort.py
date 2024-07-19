import typing

from hpotk import MinimalOntology
from jinja2 import Environment, PackageLoader

from genophenocorr.model import Cohort, Variant, VariantFormatter


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
    ):
        """
        Args:
            hpo(MinimalOntology): An HPO ontology object from hpo-toolkit
            top_phenotype_count(int): Maximum number of HPO terms to display in the HTML table (default: 10)
            top_variant_count(int): Maximum number of variants to display in the HTML table (default: 10)
        """
        self._hpo = hpo
        self._top_phenotype_count = top_phenotype_count
        self._top_variant_count = top_variant_count
        environment = Environment(loader=PackageLoader('genophenocorr.view', 'templates'))
        self._cohort_template = environment.get_template("cohort.html")

    def process(
            self,
            cohort: Cohort,
            transcript_id: typing.Optional[str] = None
    ) -> str:
        """
        Create an HTML that should be shown with display(HTML(..)) of the ipython package.

        Args:
            cohort (Cohort): The cohort being analyzed in the current Notebook
            transcript_id (str): the transcript that we map variants onto

        Returns:
            str: an HTML string with parameterized template for rendering
        """
        context = self._prepare_context(cohort, transcript_id=transcript_id)
        return self._cohort_template.render(context)

    def _prepare_context(
            self,
            cohort: Cohort,
            transcript_id: str
    ) -> typing.Mapping[str, typing.Any]:

        hpo_counts = list()
        for hpo in cohort.list_present_phenotypes(top=self._top_phenotype_count):
            hpo_id = hpo[0]
            individual_count = hpo[1]
            hpo_label = "n/a"
            if hpo_id in self._hpo:
                hpo_label = self._hpo.get_term(hpo_id).name
            hpo_counts.append({"HPO": hpo_label, "ID": hpo_id, "Count": individual_count})

        var_counts = list()
        variant_to_display_d = self.get_variant_description(cohort, transcript_id)
        for var in cohort.list_all_variants(top=self._top_variant_count):
            chrom_var = var[0]
            var_count = var[1]
            # get HGVS or human readable variant
            var_name = variant_to_display_d.get(chrom_var, chrom_var)
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
        if transcript_id is not None:
            has_transcript = True
            data_by_tx = cohort.variant_effect_count_by_tx(tx_id=transcript_id)
            # e.g., data structure -- {'effect}': 'FRAMESHIFT_VARIANT', 'count': 175}, {'effect}': 'STOP_GAINED', 'count': 67},
            for tx_id, counter in data_by_tx.items():
                if tx_id == transcript_id:
                    for effect, count in counter.items():
                        var_effects_list.append({"effect": effect, "count": count})
        else:
            has_transcript = False
        if transcript_id is None:
            transcript_id = "MANE transcript ID"
        # The following dictionary is used by the Jinja2 HTML template
        return {
            "n_individuals": len(cohort.all_patients),
            "n_excluded": cohort.get_excluded_count(),
            "total_hpo_count": len(cohort.all_phenotypes()),
            "top_hpo_count": self._top_phenotype_count,
            "hpo_counts": hpo_counts,
            "unique_variant_count": len(cohort.all_variants()),
            "top_var_count": self._top_variant_count,
            "var_counts": var_counts,
            "n_diseases": n_diseases,
            "disease_counts": disease_counts,
            "has_transcript": has_transcript,
            "var_effects_list": var_effects_list,
            "transcript_id": transcript_id,
        }

    @staticmethod
    def get_variant_description(
            cohort: Cohort,
            transcript_id: typing.Optional[str],
            only_hgvs: bool = True,
    ) -> typing.Mapping[str, str]:
        """
        Get user-friendly strings (e.g., HGVS for our target transcript) to match to the chromosomal strings
        Args:
            cohort (Cohort): The cohort being analyzed in the current Notebook
            transcript_id (str): the transcript that we map variants onto
            only_hgvs (bool): do not show the transcript ID part of the HGVS annotation, just the annotation.

        Returns:
            typing.Mapping[str, str]: key: chromosomal, value: display (e.g. HGVS) string of variant
        """
        chrom_to_display = dict()
        all_var_set = cohort.all_variants()
        var_formatter = VariantFormatter()
        for var in all_var_set:
            var_string = var.variant_coordinates.variant_key
            display = var_formatter.format_as_string(var, transcript_id)
            if only_hgvs:
                # do not show the transcript id
                fields = display.split(":")
                if len(fields) > 1:
                    display = fields[1]
                else:
                    display = [0]
            chrom_to_display[var_string] = display
        return chrom_to_display
