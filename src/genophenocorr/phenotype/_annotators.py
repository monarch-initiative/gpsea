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
    """A class that creates a Phenotype object 

    Methods:
        create_phenotype(term_ids:Iterable[Tuple[str, bool]]): Creates a list of Phenotype objects from a list of tuples.
                                                               Each tuple has the HPO ID and a boolean on if the phenotype is observed.
    """
    def __init__(self, hpo: hpotk.ontology.Ontology,
                 validator: hpotk.validate.ValidationRunner):
        """Constructs all necessary attributes for a PhenotypeCreator object

        Args:
            hpo (hpotk.ontology.Ontology): An Ontology object 
            validator (hpotk.validate.ValidationRunner): A ValidationRunner object 
        """
        self._logger = logging.getLogger(__name__)
        self._hpo = hpotk.util.validate_instance(hpo, hpotk.ontology.Ontology, 'hpo')
        self._validator = hpotk.util.validate_instance(validator, hpotk.validate.ValidationRunner, 'validator')

    def create_phenotype(self, term_ids: typing.Iterable[typing.Tuple[str, bool]]) -> list[Phenotype]:
        """Creates a list of Phenotype objects from term IDs and checks if the term IDs satisfy the validation requirements.

        Args:
            term_ids (Iterable[Tuple[str, bool]]): A list of Tuples, structured (HPO IDs, boolean- True if observed)
        Returns:
            list[Phenotype]: A list of Phenotype objects
        Error:
            PhenotypeValidationException: An instance of an issue with the ValidationRunner
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
