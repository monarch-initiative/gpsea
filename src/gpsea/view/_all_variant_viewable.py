import typing

from hpotk import MinimalOntology
from jinja2 import Environment, PackageLoader
from collections import namedtuple, defaultdict

from gpsea.model import Cohort
from gpsea.model._variant import Variant
from ._formatter import VariantFormatter


ToDisplay = namedtuple('ToDisplay', ['hgvs_cdna', 'hgvsp', 'variant_effects'])

VariantData = namedtuple('VariantData', ['variant_key', 'hgvs_cdna', 'hgvsp', 'chromosomal_description', 'variant_effects'] )

class AllVariantViewable:
    """
    Class to create a viewable object that is uses a Jinja2 template to create an HTML element
    for display in the Jupyter notebook.
    """

    def __init__(
            self,
            transcript_id: str
    ):
        """
        Args:
            hpo(MinimalOntology): An HPO ontology object from hpo-toolkit
            top_phenotype_count(int): Maximum number of HPO terms to display in the HTML table (default: 10)
            top_variant_count(int): Maximum number of variants to display in the HTML table (default: 10)
        """
        environment = Environment(loader=PackageLoader('gpsea.view', 'templates'))
        self._cohort_template = environment.get_template("all_variants.html")
        self._var_formatter = VariantFormatter(transcript_id)

    def process(
        self,
        cohort: Cohort,
        only_hgvs: bool = True
    ) -> str:
        """
        Create an HTML that should be shown with display(HTML(..)) of the ipython package.

        Args:
            cohort (Cohort): The cohort being analyzed in the current Notebook
            transcript_id (str): the transcript that we map variants onto

        Returns:
            str: an HTML string with parameterized template for rendering
        """
        context = self._prepare_context(cohort, transcript_id=self._transcript_id, only_hgvs=only_hgvs)
        return self._cohort_template.render(context)

    def _prepare_context(
        self,
        cohort: Cohort,
        transcript_id: typing.Optional[str],
        only_hgvs
    ) -> typing.Mapping[str, typing.Any]:
        variant_count_dictionary = defaultdict()
        for var in cohort.all_variants():
            vdata = self._get_variant_data(var, only_hgvs)
            variant_count_dictionary[vdata] += 1
        variant_counts = list()
        for var_data, count in variant_count_dictionary.items():
            variant_counts.append(
                {
                    "variant_key": var_data.variant_key,
                    "variant": var_data.hgvs_cdna, 
                    "variant_name": var_data.hgvs_cdna, 
                    "protein_name": var_data.hgvsp, 
                    "variant_effects": var_data.variant_effects,
                    "count": count,
                }
            )
        # The following dictionary is used by the Jinja2 HTML template
        return {
            "variant_count_list": variant_counts,
            "total_unique_allele_count": len(variant_counts)
        }


    def _get_variant_description(
        self,
        variant: Variant,
        only_hgvs: bool
    ) -> VariantData:
        """
        Get user-friendly strings (e.g., HGVS for our target transcript) to match to the chromosomal strings
        Args:
            cohort (Cohort): The cohort being analyzed in the current Notebook
            transcript_id (str): the transcript that we map variants onto
            only_hgvs (bool): do not show the transcript ID part of the HGVS annotation, just the annotation.

        Returns:
            typing.Mapping[str, ToDisplay]: key: variant key, value: namedtuple(display (e.g. HGVS) string of variant, hgvsp protein string of variant)
        """
        variant_key = variant.variant_info.variant_key
        display = self._var_formatter.format_as_string(variant)
        tx_annotation = variant.get_tx_anno_by_tx_id(self._transcript_id)
        if tx_annotation is not None:
            hgvsp = tx_annotation.hgvsp
            var_effects = [var_eff.name for var_eff in tx_annotation.variant_effects]
        else:
            hgvsp = None
            var_effects = None
        if only_hgvs:
            # do not show the transcript id
            fields_dna = display.split(":")
            fields_ps = hgvsp.split(":") if hgvsp is not None else [None]
            if len(fields_dna) > 1:
                display = fields_dna[1]
            else:
                display = fields_dna[0]
            if len(fields_ps) > 1:
                hgvsp = fields_ps[1]
            else:
                hgvsp = fields_ps[0]
            return VariantData(variant_key=variant_key, hgvs_cdna = fields_dna, hgvsp = hgvsp, chromosomal_description=None, variant_effects=var_effects)
