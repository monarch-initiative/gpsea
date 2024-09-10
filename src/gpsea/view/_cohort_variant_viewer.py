import typing

from jinja2 import Environment, PackageLoader
from collections import namedtuple, defaultdict

from gpsea.model import Cohort, Variant, VariantEffect
from ._formatter import VariantFormatter


ToDisplay = namedtuple('ToDisplay', ['hgvs_cdna', 'hgvsp', 'variant_effects'])

VariantData = namedtuple('VariantData', ['variant_key', 'hgvs_cdna', 'hgvsp', 'variant_effects'])


class CohortVariantViewer:
    """
    `CohortVariantViewer` creates an HTML report with the cohort variants.

    The report can be either written into an HTML file or displayed in a Jupyter notebook.

    See :ref:show-cohort-variants: for an example usage.
    """

    def __init__(
        self,
        tx_id: str
    ):
        """
        Args:
            tx_id (str): The transcript identifier (Usually, the MANE RefSeq transcript, that should start with "NM_")
        """
        environment = Environment(loader=PackageLoader('gpsea.view', 'templates'))
        self._cohort_template = environment.get_template("all_variants.html")
        self._var_formatter = VariantFormatter(tx_id)
        if not tx_id.startswith("NM"):
            print(f"[WARNING] Non-RefSeq transcript id: {tx_id}")
        self._transcript_id = tx_id

    def process(
        self,
        cohort: Cohort,
        only_hgvs: bool = True
    ) -> str:
        """
        Create an HTML that should be shown with ``display(HTML(..))`` of the ipython package.

        Args:
            cohort (Cohort): The cohort being analyzed in the current notebook.
            only_hgvs (bool): Do not show the transcript ID part of the HGVS annotation, just the annotation.

        Returns:
            str: an HTML string with parameterized template for rendering
        """
        context = self._prepare_context(cohort, only_hgvs=only_hgvs)
        return self._cohort_template.render(context)

    def _prepare_context(
        self,
        cohort: Cohort,
        only_hgvs: bool,
    ) -> typing.Mapping[str, typing.Any]:
        variant_count_dictionary = defaultdict(int)
        variant_effect_count_dictionary = defaultdict(int)
        variant_key_to_variant = dict()
        for var in cohort.all_variants():
            var_key = var.variant_info.variant_key
            vdata = self._get_variant_data(var, only_hgvs)
            variant_key_to_variant[var_key] = vdata
            variant_count_dictionary[var_key] += 1
            for v_eff in vdata.variant_effects:
                variant_effect_count_dictionary[v_eff] += 1
        variant_counts = list()
        variant_effect_counts = list()
        for var_key, count in variant_count_dictionary.items():
            var_data = variant_key_to_variant[var_key]
            variant_counts.append(
                {
                    "variant_key": var_data.variant_key,
                    "variant": var_data.hgvs_cdna,
                    "variant_name": var_data.hgvs_cdna,
                    "protein_name": var_data.hgvsp,
                    "variant_effects": ", ".join(var_data.variant_effects),
                    "count": count,
                }
            )
        for v_effect, count in variant_effect_count_dictionary.items():
            variant_effect_counts.append(
                {
                    "effect": v_effect,
                    "count": count
                }
            )
        variant_counts = sorted(variant_counts, key=lambda row: row.get("count"), reverse=True)
        variant_effect_counts = sorted(variant_effect_counts, key=lambda row: row.get("count"), reverse=True)

        # The following dictionary is used by the Jinja2 HTML template
        return {
            "has_transcript": False,
            "variant_count_list": variant_counts,
            "variant_effect_count_list": variant_effect_counts,
            "total_unique_allele_count": len(variant_counts)
        }

    def _get_variant_data(
        self,
        variant: Variant,
        only_hgvs: bool
    ) -> VariantData:
        """
        Get user-friendly strings (e.g., HGVS for our target transcript) to match to the chromosomal strings
        Args:
            variant (Variant): The variant to be formatted.
            only_hgvs (bool): do not show the transcript ID part of the HGVS annotation, just the annotation.

        Returns:
            VariantData: a named tuple with variant data formatted for human consumption.
        """
        variant_key = variant.variant_info.variant_key
        if variant.variant_info.has_sv_info():
            sv_info = variant.variant_info.sv_info
            gene_symbol = sv_info.gene_symbol
            display = f"SV involving {gene_symbol}"
            effect = VariantEffect.structural_so_id_to_display(so_term=sv_info.structural_type)
            return VariantData(
                variant_key=variant_key,
                hgvs_cdna=display,
                hgvsp="p.?",
                variant_effects=[effect],
            )
        else:
            variant_key = variant.variant_info.variant_key
            display = self._var_formatter.format_as_string(variant)
            tx_annotation = variant.get_tx_anno_by_tx_id(self._transcript_id)
            if tx_annotation is not None:
                hgvsp = tx_annotation.hgvsp
                var_effects = [var_eff.to_display() for var_eff in tx_annotation.variant_effects]
            else:
                hgvsp = None
                var_effects = []
            if only_hgvs:
                # do not show the transcript id
                fields_dna = display.split(":")
                if len(fields_dna) > 1:
                    display_hgvs_cDNA = fields_dna[1]
                else:
                    display_hgvs_cDNA = fields_dna[0]
                
                fields_ps = hgvsp.split(":") if hgvsp is not None else (None,)
                if len(fields_ps) > 1:
                    hgvsp = fields_ps[1]
                else:
                    hgvsp = fields_ps[0]
            return VariantData(
                variant_key=variant_key,
                hgvs_cdna=display_hgvs_cDNA,
                hgvsp=hgvsp,
                variant_effects=var_effects,
            )
