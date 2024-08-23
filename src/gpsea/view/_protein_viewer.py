import typing

from dataclasses import dataclass

from jinja2 import Environment, PackageLoader
from collections import namedtuple

from gpsea.model import Cohort
from gpsea.model.genome import Region
from ._protein_visualizable import ProteinVisualizable

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


class ProteinViewable:
    """
    Class to create a pretty HTML table to display the protein information in the Jupyter notebook.
    """
    def __init__(self) -> None:
        environment = Environment(loader=(PackageLoader('gpsea.view', 'templates')))
        self._cohort_template = environment.get_template("protein.html")
        
    def process(self, cohort: Cohort, pvis: ProteinVisualizable) -> str:
        """
        Summarize the data regarding the protein into a HTML table.

        Args:
            cohort (Cohort): the cohort of patients being analyzed
            pvis (ProteinVisualizable): The class that collects data from the UniProt API for a given protein ID

        Returns:
            str: an HTML document for showing in Jupyter notebook
        """
        context = self._prepare_context(cohort, pvis)
        return self._cohort_template.render(context)
    
    def _prepare_context(self, cohort: Cohort, pvis: ProteinVisualizable) -> typing.Mapping[str, typing.Any]:
        protein_id = pvis.protein_id
        
        protein_features = []
        
        for i in range(len(pvis.protein_feature_names)):
            feature = Feature(
                name=pvis.protein_feature_names[i], 
                type=pvis.protein_feature_types[i], 
                region=Region(pvis.protein_feature_starts[i], pvis.protein_feature_ends[i]), 
                variant_count=0,
            )
            protein_features.append(feature)
            
        for feature in protein_features:
            count = 0
            for var in cohort.all_variants():
                tx_anno = var.get_tx_anno_by_tx_id(pvis.transcript_id)
                if tx_anno is not None:
                    location = tx_anno.protein_effect_location
                    if location is not None and location.overlaps_with(feature.region):
                        count += 1
            
            feature.variant_count = count
            
        final_protein_features = sorted(protein_features, key=lambda f: f.region.start)

        return {
            'protein_id': protein_id,
            'protein_label': pvis.protein_metadata.label,
            'protein_features': final_protein_features
        }
    
