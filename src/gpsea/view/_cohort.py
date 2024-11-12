import typing

from collections import namedtuple, defaultdict

import hpotk

from gpsea.model import Cohort, Sex
from ._base import BaseViewer
from ._report import GpseaReport, HtmlGpseaReport
from ._formatter import VariantFormatter

ToDisplay = namedtuple('ToDisplay', ['hgvs_cdna', 'hgvsp', 'variant_effects'])


class CohortViewer(BaseViewer):
    """
    `CohortViewer` summarizes the most salient :class:`~gpsea.model.Cohort` aspects into an HTML report.
    """

    def __init__(
        self,
        hpo: hpotk.MinimalOntology,
        top_phenotype_count: int = 10,
        top_variant_count: int = 10,
    ):
        """
        Args:
            hpo(MinimalOntology): An HPO ontology object from hpo-toolkit
            top_phenotype_count(int): Maximum number of HPO terms to display in the HTML table (default: 10)
            top_variant_count(int): Maximum number of variants to display in the HTML table (default: 10)
        """
        super().__init__()
        self._hpo = hpo
        self._top_phenotype_count = top_phenotype_count
        self._top_variant_count = top_variant_count
        self._cohort_template = self._environment.get_template("cohort.html")

    def process(
        self,
        cohort: Cohort,
        transcript_id: typing.Optional[str] = None,
    ) -> GpseaReport:
        """
        Generate the report for a given `cohort`.

        Args:
            cohort (Cohort): The cohort being analyzed in the current Notebook
            transcript_id (str): the transcript that we map variants onto

        Returns:
            GpseaReport: a report that can be stored to a path or displayed in
                interactive environment such as Jupyter notebook.
        """
        context = self._prepare_context(cohort, transcript_id=transcript_id)
        report = self._cohort_template.render(context)
        return HtmlGpseaReport(html=report)

    def _prepare_context(
        self,
        cohort: Cohort,
        transcript_id: typing.Optional[str],
    ) -> typing.Mapping[str, typing.Any]:

        hpo_counts = list()
        for term_id, count in cohort.list_present_phenotypes(top=self._top_phenotype_count):
            label = self._hpo.get_term_name(term_id)
            if label is None:
                label = 'N/A'

            hpo_counts.append(
                {
                    "label": label,
                    "term_id": term_id,
                    "count": count,
                }
            )

        variant_counts = list()
        variant_to_display_d = CohortViewer._get_variant_description(cohort, transcript_id)
        for variant_key, count in cohort.list_all_variants(top=self._top_variant_count):
            # get HGVS or human readable variant
            if variant_key in variant_to_display_d:
                display = variant_to_display_d[variant_key]
                hgvsc = display.hgvs_cdna
                hgvsp = display.hgvsp
                effects = '' if display.variant_effects is None else ", ".join(display.variant_effects)
            else:
                hgvsc = ''
                hgvsp = ''
                effects = ''

            variant_counts.append(
                {
                    "count": count,
                    "key": variant_key,
                    "hgvsc": hgvsc,
                    "hgvsp": hgvsp,
                    "effects": effects,
                }
            )

        disease_counts = list()
        for disease_id, disease_count in cohort.list_all_diseases():
            label = "Unknown"
            for disease in cohort.all_diseases():
                if disease.identifier.value == disease_id:
                    label = disease.name
                    break
                
            disease_counts.append(
                {
                    "disease_id": disease_id,
                    "label": label,
                    "count": disease_count,
                }
            )

        variant_effects = list()
        has_transcript = transcript_id is not None
        if has_transcript:
            var_effects_d = dict()
            for tx_id, counter in cohort.variant_effect_count_by_tx(tx_id=transcript_id).items():
                # e.g., data structure
                #   -- {'effect}': 'FRAMESHIFT_VARIANT', 'count': 175},
                #   -- {'effect}': 'STOP_GAINED', 'count': 67},
                if tx_id == transcript_id:
                    for effect, count in counter.items():
                        var_effects_d[effect] = count
                    break
            total = sum(var_effects_d.values())
            # Sort in descending order based on counts
            for effect, count in sorted(var_effects_d.items(), key=lambda item: item[1], reverse=True):
                variant_effects.append(
                    {
                        "effect": effect,
                        "count": count,
                        "percent": round(count / total * 100),
                    }
                )
        else:
            transcript_id = "MANE transcript ID"

        n_has_disease_onset = 0
        hpo_id_to_has_cnset_count_d = defaultdict(int)
        for pat in cohort.all_patients:
            diseases = pat.diseases
            for d in diseases:
                if d.onset is not None:
                    n_has_disease_onset += 1
                    break  # for now, do this for any diseases. All of our current phenopackets habve but one disease
            terms_and_ancestors_with_onset_information = set()
            for pf in pat.present_phenotypes():
                if pf.onset is not None:
                    ancs = self._hpo.graph.get_ancestors(pf.identifier, include_source=True)
                    terms_and_ancestors_with_onset_information.update(ancs)
            for hpo_id in terms_and_ancestors_with_onset_information:
                hpo_id_to_has_cnset_count_d[hpo_id] += 1
        
        # When we get here, we want to present the counts of HPO terms that have onset information
        # onset_threshold = int(0.2 * len(cohort))  # only show terms with a decent amount of information
        # has_onset_information = list()
        # for hpo_id, count in hpo_id_to_has_cnset_count_d.items():
        #     if count < onset_threshold:
        #         continue
        #     label = self._hpo.get_term_name(hpo_id)
        #     if label is not None:
        #         display_label = f"{label} ({hpo_id})"
        #     else:
        #         display_label = hpo_id
        #     has_onset_information.append(
        #         {"HPO": display_label, "count": count})
        # n_has_onset_info = len(has_onset_information)

        # The following dictionary is used by the Jinja2 HTML template
        return {
            "cohort": cohort,
            "transcript_id": transcript_id,
            "top_hpo_count": self._top_phenotype_count,
            "hpo_counts": hpo_counts,
            "top_var_count": self._top_variant_count,
            "var_counts": variant_counts,
            # "n_diseases": n_diseases,
            "disease_counts": disease_counts,
            "has_transcript": has_transcript,
            "variant_effects": variant_effects,
        }

    @staticmethod
    def _get_variant_description(
        cohort: Cohort,
        transcript_id: typing.Optional[str],
        only_hgvs: bool = True,
    ) -> typing.Mapping[str, ToDisplay]:
        """
        Get user-friendly strings (e.g., HGVS for our target transcript) to match to the chromosomal strings
        Args:
            cohort (Cohort): The cohort being analyzed in the current Notebook
            transcript_id (str): the transcript that we map variants onto
            only_hgvs (bool): do not show the transcript ID part of the HGVS annotation, just the annotation.

        Returns:
            typing.Mapping[str, ToDisplay]:
              key: variant key,
              value: namedtuple(display (e.g. HGVS) string of variant, hgvsp protein string of variant)
        """
        chrom_to_display = dict()
        var_formatter = VariantFormatter(transcript_id)

        for var in cohort.all_variants():
            variant_key = var.variant_info.variant_key
            display = var_formatter.format_as_string(var)
            if transcript_id is None:
                tx_annotation = None
            else:
                tx_annotation = var.get_tx_anno_by_tx_id(transcript_id)

            if tx_annotation is None:
                hgvsp = None
                var_effects = None
            else:
                hgvsp = tx_annotation.hgvsp
                var_effects = [var_eff.name for var_eff in tx_annotation.variant_effects]

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
            chrom_to_display[variant_key] = ToDisplay(display, hgvsp, var_effects)

        return chrom_to_display
