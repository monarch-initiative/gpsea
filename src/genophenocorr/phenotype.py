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
               f"name={self.name})"

    def __repr__(self):
        return str(self)


class PhenotypeValidationException(BaseException):

    def __init__(self, issues):
        self._issues = issues

    @property
    def issues(self):
        return self._issues


class PhenotypeCreator:

    def __init__(self, hpo: hpotk.ontology.MinimalOntology,
                 validator: hpotk.validate.ValidationRunner):
        self._logger = logging.getLogger(__name__)
        if not isinstance(hpo, hpotk.ontology.MinimalOntology):
            raise ValueError(f'hpo must be hpotk.ontology.MinimalOntology but was {type(hpo)}')
        self._hpo = hpo
        if not isinstance(validator, hpotk.validate.ValidationRunner):
            raise ValueError(f'validator must be hpotk.validate.ValidationRunner but was {type(validator)}')
        self._validator = validator

    def create_phenotype(self, term_ids: typing.Iterable[str]) -> list[Phenotype]:
        """
        Create a list of Phenotype instances from term IDs and check if the term IDs satisfy the validation requirements.

        The method returns a list if the term IDs are valid or raises :class:`PhenotypeValidationException`
        with the validation issues.
        """
        terms = []
        for term_id in term_ids:
            term = self._hpo.get_term(term_id)
            if term is None:
                raise ValueError(f'Term ID {term_id} is not present in HPO v{self._hpo.version}')
            if not isinstance(term, hpotk.model.Term):
                raise ValueError(f"Term must be type Term but is type {type(term)}")
            terms.append(term)
        validation_results = self._validator.validate_all(terms)
        if validation_results.is_ok:
            return [Phenotype(term) for term in terms]
        else:
            # We return the messages for now. We may provide more details in future, if necessary.
            issues = [r.message for r in validation_results.results]
            raise PhenotypeValidationException(issues)
