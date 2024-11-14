import typing

from dataclasses import dataclass


from collections import defaultdict

from gpsea.model import Cohort, ProteinMetadata
from gpsea.model.genome import Region

from ._base import BaseViewer
from ._report import GpseaReport, HtmlGpseaReport


@dataclass(frozen=False)
class Feature:
    """
    A private dataclass for representing a table row.

    Any edits to the dataclass must also be followed by an update of the Jinja template.
    """
    name: str
    type: str
    region: Region
    variant_count: int


class ProteinVariantViewer(BaseViewer):
    """
    Class to create a pretty HTML table to display the protein information in the Jupyter notebook.
    """
    def __init__(
        self,
        protein_metadata: ProteinMetadata,
        tx_id: str,
    ):
        super().__init__()
        self._cohort_template = self._environment.get_template("protein.html")
        self._protein_meta = protein_metadata
        self._tx_id = tx_id

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
        protein_features = sorted(self._protein_meta.protein_features, key=lambda f: f.info.start)

        # collect variants that are located in the protein features as well as other variants that are located
        # "in-between" the features
        feature_to_variants = defaultdict(list)
        non_feature_count = 0
        n_variants_in_features = 0

        # This iterates over variants as recorded in individuals,
        # not over *unique* `VariantInfo`s
        for var in cohort.all_variants():
            target_annot = next((x for x in var.tx_annotations if x.transcript_id == self._tx_id), None)
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
                    'name': feature.info.name,
                    'type': feature.feature_type.name,
                    'start': feature.info.start,
                    'end': feature.info.end,
                    'count': len(feature_to_variants[feature]),
                    'variants': variant_list,
                }
            )

        return {
            'protein_label': protein_label,
            'protein_id': protein_id,
            'protein_length': self._protein_meta.protein_length,
            'feature_list': feature_list,
            'n_variants_in_features': n_variants_in_features,
            'non_feature_count': non_feature_count,
        }

