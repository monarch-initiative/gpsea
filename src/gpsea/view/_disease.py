import typing

from hpotk import MinimalOntology
from jinja2 import Environment, PackageLoader
from collections import Counter, defaultdict

from gpsea.model import Cohort
from ._base import GpseaReport, HtmlGpseaReport


class DiseaseViewer:
    """
    TODO
    """
    def __init__(
        self,
        hpo: MinimalOntology,
        transcript_id: typing.Optional[str] = None
    ):
        self._hpo = hpo
        self._tx_id = transcript_id
        environment = Environment(loader=(PackageLoader('gpsea.view', 'templates')))
        self._cohort_template = environment.get_template("disease.html")
        
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
                    variant_key_count = disease_variants_counts[disease.identifier.value]
                    variant_key_count.update(variant.variant_info.variant_key for variant in pat.variants)
                    for variant in pat.variants:
                        for txa in variant.tx_annotations:
                            if self._tx_id == txa.transcript_id:
                                disease_effects_counts[curie].update(eff.name for eff in txa.variant_effects)
            
        return {
            'n_diseases': n_diseases,
            'disease_counts': disease_counts,
            'disease_variant_counts': disease_variants_counts,
            'disease_effects_counts': disease_effects_counts
        }
    
