import hpotk


class Phenotype(hpotk.model.Identified, hpotk.model.Named):

    def __init__(self, term_id: hpotk.model.TermId, name: str, excluded=False) -> None:
        # TODO(lnrekerle) - instance/type check
        self._term_id = term_id
        self._name = name
        self._excluded = excluded
        # TODO - the checks below can be done but this is not the best place. The phenotype features should
        #  be checked by the code that loads phenopacket/cohort or whatever else and complain if the input
        #  does not look good.
        # BUG - the line below is equivalent to (not id.startswith('HP:')) and len(id)==10:
        #  which is likely NOT the expected logic
        # if not id.startswith("HP:") and len(id) == 10:
        #     raise ValueError(f"Malformed HPO id: {id}")

        # if obso_results.is_ok: ## TODO: Add an elif for using an obsolete term and print a WARNING
        #     self._hpo_term = Term
        #     self._ancestors = hpotk.algorithm.get_ancestors(ont, self.id)
        #     self._descendants = hpotk.algorithm.get_descendants(ont, self.id)
        #     self._excluded = excluded
        #     self._measured = measured

    @property
    def identifier(self) -> hpotk.model.TermId:
        return self._term_id

    @property
    def name(self) -> str:
        return self._name

    # @property
    # def id(self):
    #     return self._hpo_term.identifier.value
    #
    # @property
    # def label(self):
    #     return self._hpo_term.name

    @property
    def excluded(self) -> bool:
        return self._excluded

    def __eq__(self, other):
        return isinstance(other, Phenotype) \
            and self.identifier == other.identifier \
            and self.name == other.name \
            and self.excluded == other.excluded

    def __hash__(self):
        return hash((self.identifier, self.name, self.excluded))

    def __str__(self):
        return f"Phenotype(identifier={self.identifier}, " \
               f"name={self.name}, " \
               f"excluded={self.excluded})"

    def __repr__(self):
        return str(self)

    # @property
    # def full_definition(self):
    #     return self._hpo_term.definition
    
    # @property
    # def alt_ids(self):
    #     return self._hpo_term.alt_term_ids

    # @property
    # def ancestors(self):
    #     return self._ancestors
    #
    # @property
    # def descendants(self):
    #     return self._descendants
    #
    # def list_ancestors(self):
    #     anc = []
    #     for a in self.ancestors:
    #         anc.append(a.value)
    #     return anc
    #
    # def list_descendants(self):
    #     desc = []
    #     for d in self.descendants:
    #         desc.append(d.value)
    #     return desc


# class PhenotypeCreator:
#
#     def __init__(self, hpo: hpotk.ontology.MinimalOntology):
#         if not isinstance(hpo, hpotk.ontology.MinimalOntology):
#             raise ValueError(f'hpo must be hpotk.ontology.MinimalOntology but was {type(hpo)}')
#         self._hpo = hpo
#         self._validator = hpotk.validate.ObsoleteTermIdsValidator(hpo)
#
#     def create_phenotype(self, id, excluded=False, measured=True) -> Phenotype:
#         term = ont.get_term(id)
#         validation_results = self._validator.validate([term])
#         if validation_results.is_ok: ## TODO: Add an elif for using an obsolete term and print a WARNING
#             self._hpo_term = term
#             self._ancestors = hpotk.algorithm.get_ancestors(ont, self.id)
#             self._descendants = hpotk.algorithm.get_descendants(ont, self.id)
#             self._excluded = excluded
#             self._measured = measured


