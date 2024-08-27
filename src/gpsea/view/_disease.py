import typing

from hpotk import MinimalOntology
from jinja2 import Environment, PackageLoader
from collections import Counter, defaultdict


class DiseaseViewable:
    """
    TODO
    """
    def __init__(self, hpo: MinimalOntology, transcript_id: str = None, ) -> None:
        self._hpo = hpo
        self._tx_id = transcript_id
        environment = Environment(loader=(PackageLoader('gpsea.view', 'templates')))
        self._cohort_template = environment.get_template("disease.html")
        
    def process(self, cohort) -> str:
        context = self._prepare_context(cohort)
        return self._cohort_template.render(context)
    
    def _prepare_context(self, cohort) -> typing.Mapping[str, typing.Any]:
        diseases = cohort.list_all_diseases()
        n_diseases = len(diseases)
        disease_counts = list()
        disease_variants_counts = defaultdict(Counter)
        disease_effects_counts = {}
        for d in diseases:
            disease_id = d[0]
            disease_count = d[1]
            disease_name = "Unknown"
            disease_variants_counts[disease_id.value] = Counter()
            disease_effects_counts[disease_id.value] = Counter()
            for dis in cohort.all_diseases():
                if dis.identifier == d[0]:
                    disease_name = dis.name
            for pat in cohort.all_patients:
                if disease_id in [disease.identifier for disease in pat.diseases]:
                    disease_variants_counts[disease_id.value].update(variant.variant_coordinates.variant_key for variant in pat.variants)
                    for var_effect in [txa.variant_effects for variant in pat.variants for txa in variant.tx_annotations if txa.transcript_id == self._tx_id]:
                        disease_effects_counts[disease_id.value].update(var_eff.name for var_eff in var_effect)
            disease_counts.append({"disease_id": disease_id, "disease_name": disease_name, "count": disease_count})
            
        return {
            'n_diseases': n_diseases,
            'disease_counts': disease_counts,
            'disease_variant_counts': disease_variants_counts,
            'disease_effects_counts': disease_effects_counts
        }
    
