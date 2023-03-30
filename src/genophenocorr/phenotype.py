import hpotk
from hpotk.ontology import Ontology
from hpotk.ontology.load.obographs import load_ontology
from hpotk.validate import ObsoleteTermIdsValidator
ont = load_ontology('https://raw.githubusercontent.com/obophenotype/human-phenotype-ontology/master/hp.json')

class Phenotype:
    def __init__(self, id, excluded=False, measured=True) -> None:
        if not id.startswith("HP:") and len(id) == 10:
            raise ValueError(f"Malformed HPO id: {id}")
        obso = ObsoleteTermIdsValidator(ont)
        Term = ont.get_term(id)
        obso_results = obso.validate([Term])
        if obso_results.is_ok: ## TODO: Add an elif for using an obsolete term and print a WARNING
            self._hpo_term = Term
            self._ancestors = hpotk.algorithm.get_ancestors(ont, self.id)
            self._descendants = hpotk.algorithm.get_descendants(ont, self.id)
            self._excluded = excluded
            self._measured = measured
        
    @property
    def id(self):
        return self._hpo_term.identifier.value
    
    @property
    def label(self):
        return self._hpo_term.name

    @property
    def excluded(self):
        return self._excluded

    @property
    def full_definition(self):
        return self._hpo_term.definition
    
    @property
    def alt_ids(self):
        return self._hpo_term.alt_term_ids

    @property
    def ancestors(self):
        return self._ancestors

    @property
    def descendants(self):
        return self._descendants

    def list_ancestors(self):
        anc = []
        for a in self.ancestors:
            anc.append(a.value)
        return anc
    
    def list_descendants(self):
        desc = []
        for d in self.descendants:
            desc.append(d.value)
        return desc
