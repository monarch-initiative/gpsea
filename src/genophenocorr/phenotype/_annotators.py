import hpotk
import typing
import logging
from ._phenotype_data import Phenotype


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
        self._hpo = hpotk.util.validate_instance(hpo, hpotk.ontology.MinimalOntology, 'hpo')
        self._validator = hpotk.util.validate_instance(validator, hpotk.validate.ValidationRunner, 'validator')

    def create_phenotype(self, term_ids: typing.Iterable[typing.Tuple[str, bool]]) -> list[Phenotype]:
        """
        Create a list of Phenotype instances from term IDs and check if the term IDs satisfy the validation requirements.

        The method returns a list if the term IDs are valid or raises :class:`PhenotypeValidationException`
        with the validation issues.
        """
        terms = []
        for term_id, observed in term_ids:
            term = self._hpo.get_term(term_id)
            if term is None:
                raise ValueError(f'Term ID {term_id} is not present in HPO v{self._hpo.version}')
            terms.append((term, observed))
        validation_results = self._validator.validate_all([term[0] for term in terms])
        if validation_results.is_ok:
            return [Phenotype(term, observed) for term, observed in terms]
        else:
            # We return the messages for now. We may provide more details in future, if necessary.
            issues = [r.message for r in validation_results.results]
            raise PhenotypeValidationException(issues)
