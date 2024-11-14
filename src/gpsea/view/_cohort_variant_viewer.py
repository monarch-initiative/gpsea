import typing

from collections import namedtuple, Counter

from gpsea.model import Cohort, Variant, VariantEffect
from ._formatter import VariantFormatter
from ._base import BaseViewer
from ._report import GpseaReport, HtmlGpseaReport


ToDisplay = namedtuple('ToDisplay', ['hgvs_cdna', 'hgvsp', 'variant_effects'])

VariantData = namedtuple('VariantData', ['variant_key', 'hgvsc', 'hgvsp', 'variant_effects', 'exons'])


class CohortVariantViewer(BaseViewer):
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
        super().__init__()
        self._cohort_template = self._environment.get_template("all_variants.html")
        self._var_formatter = VariantFormatter(tx_id)
        if not tx_id.startswith("NM"):
            print(f"[WARNING] Non-RefSeq transcript id: {tx_id}")
        self._transcript_id = tx_id

    def process(
        self,
        cohort: Cohort,
        only_hgvs: bool = True
    ) -> GpseaReport:
        """
        Generate the variant report.

        Args:
            cohort (Cohort): The cohort being analyzed in the current notebook.
            only_hgvs (bool): Do not show the transcript ID part of the HGVS annotation, just the annotation.

        Returns:
            GpseaReport: a report that can be stored to a path or displayed in
                interactive environment such as Jupyter notebook.
        """
        context = self._prepare_context(cohort, only_hgvs=only_hgvs)
        html = self._cohort_template.render(context)
        return HtmlGpseaReport(html=html)

    def _prepare_context(
        self,
        cohort: Cohort,
        only_hgvs: bool,
    ) -> typing.Mapping[str, typing.Any]:
        variant_count_dictionary = Counter()
        variant_key_to_variant = dict()
        
        for var in cohort.all_variants():
            var_key = var.variant_info.variant_key
            vdata = self._get_variant_data(var, only_hgvs)
            variant_key_to_variant[var_key] = vdata
            variant_count_dictionary[var_key] += 1
        
        variant_counts = list()
        for var_key, count in variant_count_dictionary.items():
            var_data = variant_key_to_variant[var_key]
            variant_counts.append(
                {
                    "count": count,
                    "variant_key": var_data.variant_key,
                    "hgvs": f'{var_data.hgvsc} ({var_data.hgvsp})',
                    "variant_effects": ", ".join(var_data.variant_effects),
                    "exons": '-' if var_data.exons is None else ", ".join(map(str, var_data.exons)),
                }
            )
        
        variant_counts.sort(key=lambda row: row["count"], reverse=True)

        return {
            "variant_counts": variant_counts,
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
        
        if variant.variant_info.has_sv_info():
            sv_info = variant.variant_info.sv_info
            assert sv_info is not None
            gene_symbol = sv_info.gene_symbol
            display = f"SV involving {gene_symbol}"
            effect = VariantEffect.structural_so_id_to_display(so_term=sv_info.structural_type)
            return VariantData(
                variant_key=variant.variant_info.variant_key,
                hgvsc=display,
                hgvsp="p.?",
                variant_effects=(effect,),
                exons=(),
            )
        elif variant.variant_info.has_variant_coordinates():
            tx_annotation = variant.get_tx_anno_by_tx_id(self._transcript_id)
            if tx_annotation is None:
                hgvsc = '-'
                hgvsp = '-'
                var_effects = ()
                exons = ()
            else:                
                if only_hgvs:
                    if tx_annotation.hgvs_cdna is None:
                        hgvsc = '-'
                        hgvsp = '-'
                    else:
                        fields_dna = tx_annotation.hgvs_cdna.split(":")
                        if len(fields_dna) > 1:
                            hgvsc = fields_dna[1]
                        else:
                            hgvsc = fields_dna[0]
                        
                        fields_ps = tx_annotation.hgvsp.split(":") if tx_annotation.hgvsp is not None else ('-',)
                        if len(fields_ps) > 1:
                            hgvsp = fields_ps[1]
                        else:
                            hgvsp = fields_ps[0]
                else:
                    hgvsc = tx_annotation.hgvs_cdna
                
                var_effects = tuple(var_eff.to_display() for var_eff in tx_annotation.variant_effects)
                exons=tx_annotation.overlapping_exons
            
            return VariantData(
                variant_key=variant.variant_info.variant_key,
                hgvsc=hgvsc,
                hgvsp=hgvsp,
                variant_effects=var_effects,
                exons=exons
            )
        else:
            raise ValueError('Neither variant coordinates nor SV info are available')
        
        
