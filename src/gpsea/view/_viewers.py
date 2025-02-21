import typing

from collections import namedtuple, defaultdict, Counter

import hpotk

from gpsea.model import Cohort, Variant, VariantEffect, ProteinMetadata
from gpsea.analysis.pcats import HpoTermAnalysisResult
from ._base import BaseViewer, GpseaReport, HtmlGpseaReport
from ._formatter import VariantFormatter


ToDisplay = namedtuple("ToDisplay", ["hgvs_cdna", "hgvsp", "variant_effects"])
IdentifiedCount = namedtuple("IdentifiedCount", ("term_id", "label", "count"))
VariantData = namedtuple(
    "VariantData", ["variant_key", "hgvsc", "hgvsp", "variant_effects", "exons"]
)


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
            tx_id(str): the transcript accession `str` (e.g. `NM_123456.7`)
            top_phenotype_count(int): Maximum number of phenotype items (HPO terms, measurements, ...) to display in the HTML table (default: 10)
            top_variant_count(int): Maximum number of variants to display in the HTML table (default: 10)
        """
        super().__init__()
        assert isinstance(hpo, hpotk.MinimalOntology)
        self._hpo = hpo
        assert isinstance(top_phenotype_count, int)
        self._top_phenotype_count = top_phenotype_count
        assert isinstance(top_variant_count, int)
        self._top_variant_count = top_variant_count
        self._cohort_template = self._environment.get_template("cohort.html")

    def process(
        self,
        cohort: Cohort,
        transcript_id: str,
    ) -> GpseaReport:
        """
        Generate the report for a given `cohort`.

        The variant effects will be reported with respect to the provided transcript.

        Args:
            cohort (Cohort): the cohort to visualize
            transcript_id (str): the accession of the target transcript (e.g. `NM_123456.7`)

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
        hpo_counts = self._summarize_hpo_counts(cohort)

        measurement_counts = self._summarize_measurement_counts(cohort)

        disease_counts = self._summarize_disease_counts(cohort)

        variant_counts = list()
        variant_to_display_d = CohortViewer._get_variant_description(
            cohort, transcript_id
        )
        for variant_key, count in cohort.list_all_variants(top=self._top_variant_count):
            # get HGVS or human readable variant
            if variant_key in variant_to_display_d:
                display = variant_to_display_d[variant_key]
                hgvsc = display.hgvs_cdna
                hgvsp = display.hgvsp
                effects = (
                    ""
                    if display.variant_effects is None
                    else ", ".join(display.variant_effects)
                )
            else:
                hgvsc = ""
                hgvsp = ""
                effects = ""

            variant_counts.append(
                {
                    "count": count,
                    "key": variant_key,
                    "hgvsc": hgvsc,
                    "hgvsp": hgvsp,
                    "effects": effects,
                }
            )

        variant_effects = list()
        has_transcript = transcript_id is not None
        if has_transcript:
            var_effects_d = dict()
            for tx_id, counter in cohort.variant_effect_count_by_tx(
                tx_id=transcript_id
            ).items():
                # e.g., data structure
                #   -- {'effect}': 'FRAMESHIFT_VARIANT', 'count': 175},
                #   -- {'effect}': 'STOP_GAINED', 'count': 67},
                if tx_id == transcript_id:
                    for effect, count in counter.items():
                        var_effects_d[effect] = count
                    break
            total = sum(var_effects_d.values())
            # Sort in descending order based on counts
            for effect, count in sorted(
                var_effects_d.items(), key=lambda item: item[1], reverse=True
            ):
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
                    ancs = self._hpo.graph.get_ancestors(
                        pf.identifier, include_source=True
                    )
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
            "has_transcript": has_transcript,
            "top_phenotype_count": self._top_phenotype_count,
            "hpo_counts": hpo_counts,
            "measurement_counts": measurement_counts,
            "disease_counts": disease_counts,
            "top_var_count": self._top_variant_count,
            "var_counts": variant_counts,
            "variant_effects": variant_effects,
        }

    def _summarize_hpo_counts(
        self,
        cohort: Cohort,
    ) -> typing.Sequence[IdentifiedCount]:
        counts = list()
        for term_id, count in cohort.list_present_phenotypes(
            top=self._top_phenotype_count,
        ):
            label = self._hpo.get_term_name(term_id)
            if label is None:
                label = "N/A"

            counts.append(
                IdentifiedCount(
                    term_id=term_id,
                    label=label,
                    count=count,
                )
            )

        return counts

    def _summarize_measurement_counts(
        self,
        cohort: Cohort,
    ) -> typing.Sequence[IdentifiedCount]:
        measurement2label = {
            m.identifier.value: m.name for m in cohort.all_measurements()
        }

        return [
            IdentifiedCount(
                term_id=measurement_id,
                label=measurement2label.get(measurement_id, "N/A"),
                count=count,
            )
            for measurement_id, count in cohort.list_measurements(
                top=self._top_phenotype_count,
            )
        ]

    def _summarize_disease_counts(
        self,
        cohort: Cohort,
    ) -> typing.Sequence[IdentifiedCount]:
        disease2label = {m.identifier.value: m.name for m in cohort.all_diseases()}

        return [
            IdentifiedCount(
                term_id=disease_id,
                label=disease2label.get(disease_id, "N/A"),
                count=count,
            )
            for disease_id, count in cohort.list_all_diseases(
                top=self._top_phenotype_count,
            )
        ]

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
                var_effects = [
                    var_eff.name for var_eff in tx_annotation.variant_effects
                ]

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


class DiseaseViewer(BaseViewer):
    """
    TODO
    """

    def __init__(
        self, hpo: hpotk.MinimalOntology, transcript_id: typing.Optional[str] = None
    ):
        super().__init__()
        self._hpo = hpo
        self._tx_id = transcript_id
        self._cohort_template = self._environment.get_template("disease.html")

    def process(
        self,
        cohort: Cohort,
    ) -> GpseaReport:
        context = self._prepare_context(cohort)
        html = self._cohort_template.render(context)
        return HtmlGpseaReport(html=html)

    def _prepare_context(
        self,
        cohort: Cohort,
    ) -> typing.Mapping[str, typing.Any]:
        diseases = cohort.list_all_diseases()
        n_diseases = len(diseases)
        disease_counts = list()
        disease_variants_counts = defaultdict(Counter)
        disease_effects_counts = {}
        for curie, count in diseases:
            disease_name = "Unknown"
            disease_variants_counts[curie] = Counter()
            disease_effects_counts[curie] = Counter()

            for dis in cohort.all_diseases():
                if dis.identifier.value == curie:
                    disease_name = dis.name
                    break

            disease_counts.append(
                {
                    "disease_id": curie,
                    "disease_name": disease_name,
                    "count": count,
                }
            )

            for pat in cohort.all_patients:
                for disease in pat.diseases:
                    variant_key_count = disease_variants_counts[
                        disease.identifier.value
                    ]
                    variant_key_count.update(
                        variant.variant_info.variant_key for variant in pat.variants
                    )
                    for variant in pat.variants:
                        for txa in variant.tx_annotations:
                            if self._tx_id == txa.transcript_id:
                                disease_effects_counts[curie].update(
                                    eff.name for eff in txa.variant_effects
                                )

        return {
            "n_diseases": n_diseases,
            "disease_counts": disease_counts,
            "disease_variant_counts": disease_variants_counts,
            "disease_effects_counts": disease_effects_counts,
        }


class CohortVariantViewer(BaseViewer):
    """
    `CohortVariantViewer` creates an HTML report with the cohort variants.

    The report can be either written into an HTML file or displayed in a Jupyter notebook.

    See :ref:show-cohort-variants: for an example usage.
    """

    def __init__(self, tx_id: str):
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

    def process(self, cohort: Cohort, only_hgvs: bool = True) -> GpseaReport:
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
                    "hgvs": f"{var_data.hgvsc} ({var_data.hgvsp})",
                    "variant_effects": ", ".join(var_data.variant_effects),
                    "exons": (
                        "-"
                        if var_data.exons is None
                        else ", ".join(map(str, var_data.exons))
                    ),
                }
            )

        variant_counts.sort(key=lambda row: row["count"], reverse=True)

        return {
            "variant_counts": variant_counts,
        }

    def _get_variant_data(self, variant: Variant, only_hgvs: bool) -> VariantData:
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
            effect = VariantEffect.structural_so_id_to_display(
                so_term=sv_info.structural_type
            )
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
                hgvsc = "-"
                hgvsp = "-"
                var_effects = ()
                exons = ()
            else:
                if only_hgvs:
                    if tx_annotation.hgvs_cdna is None:
                        hgvsc = "-"
                        hgvsp = "-"
                    else:
                        fields_dna = tx_annotation.hgvs_cdna.split(":")
                        if len(fields_dna) > 1:
                            hgvsc = fields_dna[1]
                        else:
                            hgvsc = fields_dna[0]

                        fields_ps = (
                            tx_annotation.hgvsp.split(":")
                            if tx_annotation.hgvsp is not None
                            else ("-",)
                        )
                        if len(fields_ps) > 1:
                            hgvsp = fields_ps[1]
                        else:
                            hgvsp = fields_ps[0]
                else:
                    hgvsc = tx_annotation.hgvs_cdna

                var_effects = tuple(
                    var_eff.to_display() for var_eff in tx_annotation.variant_effects
                )
                exons = tx_annotation.overlapping_exons

            return VariantData(
                variant_key=variant.variant_info.variant_key,
                hgvsc=hgvsc,
                hgvsp=hgvsp,
                variant_effects=var_effects,
                exons=exons,
            )
        else:
            raise ValueError("Neither variant coordinates nor SV info are available")


class ProteinVariantViewer(BaseViewer):
    """
    Class to create a pretty HTML table to display the protein information in the Jupyter notebook.
    """

    def __init__(
        self,
        protein_metadata: ProteinMetadata,
        tx_id: str,
    ):
        # TODO[v1.0.0] - remove `tx_id`.
        super().__init__()
        self._cohort_template = self._environment.get_template("protein.html")
        self._protein_meta = protein_metadata

    def process(self, cohort: Cohort) -> GpseaReport:
        """
        Summarize the data regarding the protein into a HTML table.

        Args:
            cohort (Cohort): the cohort of patients being analyzed

        Returns:
            GpseaReport: a report that can be stored to a path or displayed in
                interactive environment such as Jupyter notebook.
        """
        context = self._prepare_context(cohort)
        html = self._cohort_template.render(context)
        return HtmlGpseaReport(html=html)

    def _prepare_context(self, cohort: Cohort) -> typing.Mapping[str, typing.Any]:
        protein_id = self._protein_meta.protein_id
        protein_label = self._protein_meta.label
        protein_features = sorted(
            self._protein_meta.protein_features, key=lambda f: f.info.start
        )

        # collect variants that are located in the protein features as well as other variants that are located
        # "in-between" the features
        feature_to_variants = defaultdict(list)
        non_feature_count = 0
        n_variants_in_features = 0

        # This iterates over variants as recorded in individuals,
        # not over *unique* `VariantInfo`s
        for var in cohort.all_variants():
            target_annot = next(
                (x for x in var.tx_annotations if x.protein_id == self._protein_meta.protein_id), None
            )
            if target_annot is None:
                # structural variants do not have a transcript id, and we skip them
                # It should never happen that HGVS variants do not have the proper transcript id but we do not check
                # that here in this visualization function
                continue
            hgvs_p = target_annot.hgvsp
            if hgvs_p is not None:
                # hgvs_p is now a string such as NP_001305781.1:p.Tyr15Ter. We remove the accession id, which
                # is shown in the title and caption and is redundant here
                fields = hgvs_p.split(":")
                if len(fields) == 2:
                    hgvs_p = fields[1]
            else:
                continue
            target_region = target_annot.protein_effect_location
            if target_region is None:
                # can happen for certain variant classes such as splice variant. Not an error
                continue
            found_overlap = False
            for feature in protein_features:
                if feature.info.region.overlaps_with(target_region):
                    found_overlap = True
                    feature_to_variants[feature].append(hgvs_p)

            if found_overlap:
                n_variants_in_features += 1
            else:
                non_feature_count += 1

        feature_list = list()
        for feature in protein_features:
            variant_list = "; ".join(set(feature_to_variants[feature]))
            feature_list.append(
                {
                    "name": feature.info.name,
                    "type": feature.feature_type,
                    "start": feature.info.start,
                    "end": feature.info.end,
                    "count": len(feature_to_variants[feature]),
                    "variants": variant_list,
                }
            )

        return {
            "protein_label": protein_label,
            "protein_id": protein_id,
            "protein_length": self._protein_meta.protein_length,
            "feature_list": feature_list,
            "n_variants_in_features": n_variants_in_features,
            "non_feature_count": non_feature_count,
        }


class MtcStatsViewer(BaseViewer):
    """
    `MtcStatsViewer` uses a Jinja2 template to create an HTML element for showing in the Jupyter notebook
    or for writing into a standalone HTML file.
    """

    def __init__(self):
        super().__init__()
        self._cohort_template = self._environment.get_template("stats.html")

    def process(
        self,
        result: HpoTermAnalysisResult,
    ) -> GpseaReport:
        """
        Create an HTML to present MTC part of the :class:`~gpsea.analysis.pcats.HpoTermAnalysisResult`.

        Args:
            result (HpoTermAnalysisResult): the result to show

        Returns:
            GpseaReport: a report that can be stored to a path or displayed in
                interactive environment such as Jupyter notebook.
        """
        assert isinstance(result, HpoTermAnalysisResult)
        context = self._prepare_context(result)
        html = self._cohort_template.render(context)
        return HtmlGpseaReport(html=html)

    @staticmethod
    def _prepare_context(
        report: HpoTermAnalysisResult,
    ) -> typing.Mapping[str, typing.Any]:
        counts = Counter()
        for result in report.mtc_filter_results:
            if result.is_filtered_out():
                counts[result.mtc_issue] += 1

        n_skipped = 0
        issue_to_count = list()
        for mtc_issue, count in sorted(
            counts.items(),
            key=lambda issue2count: (issue2count[0].code, issue2count[0].reason),
        ):
            issue_to_count.append(
                {
                    "code": mtc_issue.code,
                    "reason": mtc_issue.reason,
                    "doclink": mtc_issue.doclink,
                    "count": count,
                }
            )
            n_skipped += count

        n_all = len(report.phenotypes)
        n_tested = n_all - n_skipped

        # The following dictionary is used by the Jinja2 HTML template
        return {
            "mtc_name": report.mtc_correction,
            "hpo_mtc_filter_name": report.mtc_filter_name,
            "skipped_hpo_count": n_skipped,
            "tested_hpo_count": n_tested,
            "total_hpo_count": n_all,
            "issue_to_count": issue_to_count,
        }
