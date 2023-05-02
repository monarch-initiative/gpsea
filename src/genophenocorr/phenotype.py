import hpotk
import typing
import logging

class Phenotype(hpotk.model.Identified, hpotk.model.Named):

    def __init__(self, term: hpotk.model.Term) -> None:
        
        if not isinstance(term, hpotk.model.Term):
            raise ValueError(f"term must be type 'Term' not type {type(term)}")
        self._term = term
        ## May re-add if needed
        # if not isinstance(excluded, bool):
        #     raise ValueError(f"Excluded variable must be type boolean, but is type {type(excluded)}")
        # self._excluded = excluded


    @property
    def identifier(self) -> hpotk.model.TermId:
        return self._term.identifier

    @property
    def name(self) -> str:
        return self._term.name


    def __eq__(self, other):
        return isinstance(other, Phenotype) \
            and self.identifier == other.identifier \
            and self.name == other.name 

    def __hash__(self):
        return hash((self.identifier, self.name))

    def __str__(self):
        return f"Phenotype(identifier={self.identifier}, " \
               f"name={self.name}"

    def __repr__(self):
        return str(self)


class PhenotypeCreator:

    def __init__(self, hpo: hpotk.ontology.MinimalOntology):
        self._logger = logging.getLogger(__name__)
        if not isinstance(hpo, hpotk.ontology.MinimalOntology):
            raise ValueError(f'hpo must be hpotk.ontology.MinimalOntology but was {type(hpo)}')
        self._hpo = hpo
        self._validator = hpotk.validate.ObsoleteTermIdsValidator(hpo)

    def create_phenotype(self, id: str) -> Phenotype:
        term = self._hpo.get_term(id)
        if not isinstance(term, hpotk.model.Term):
            raise ValueError(f"Term must be type Term but is type {type(term)}")
        validation_results = self._validator.validate([term])
        if validation_results.is_ok: ## TODO: Add an elif for using an obsolete term and print a WARNING
            self._hpo_term = term
        elif term.is_obsolete:
            self._logger.warning('Phenotype ID %s is obsolete. Please update your IDs', term.identifier)
            self._hpo_term = term
        return Phenotype(term)



