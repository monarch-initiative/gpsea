import typing

from jinja2 import Environment, PackageLoader
from collections import namedtuple

from genophenocorr.model import Cohort
from genophenocorr.model.genome import Region
from ._protein_visualizable import ProteinVisualizable


class ProteinViewable:
    """
    Class to create a viewable table for the protein information that is uses a Jinja2 template to 
    create an HTML element for display in the Jupyter notebook.
    """
    def __init__(self) -> None:
        environment = Environment(loader=(PackageLoader('genophenocorr.view', 'templates')))
        self._cohort_template = environment.get_template("protein.html")
        
    def process(self, cohort: Cohort, pvis: ProteinVisualizable) -> str:
        """ This organizes the necessary data found through the UniProt API into a 
        easy to read table about the given protein.

        Args:
            cohort (Cohort): the cohort of patients being analyzed
            pvis (ProteinVisualizable): The class that collects data from the UniProt API for a given protein ID

        Returns:
            str: an HTML string with parameterized template for rendering
        """
        context = self._prepare_context(cohort, pvis)
        return self._cohort_template.render(context)
    
    @staticmethod
    def _get_start(feat_tuple: namedtuple) -> int:
        return feat_tuple.region.start
    
    def _prepare_context(self, cohort: Cohort, pvis: ProteinVisualizable) -> typing.Mapping[str, typing.Any]:
        protein_id = pvis.protein_id
        Feature = namedtuple('Feature', ['feature_name', 'feature_type', 'region', 'variant_count'])
        protein_features = []
        
        for i in range(len(pvis.protein_feature_names)):
            protein_features.append(Feature(pvis.protein_feature_names[i], pvis._protein_feature_types[i], Region(pvis.protein_feature_starts[i],pvis.protein_feature_ends[i]), 0))
        
        final_protein_features = []
            
        for feat_list in protein_features:
            count = 0
            for var in cohort.all_variants():
                tx_anno = var.get_tx_anno_by_tx_id(pvis.transcript_id)
                if tx_anno is not None:
                    if tx_anno.protein_effect_location is not None and tx_anno.protein_effect_location.overlaps_with(feat_list.region):
                        count += 1
            final_protein_features.append(feat_list._replace(variant_count=count))
            
            
        final_protein_features = sorted(final_protein_features, key=self._get_start)

        return {
            'protein_id': protein_id,
            'protein_label': pvis.protein_metadata.label,
            'protein_features': final_protein_features
        }
    
