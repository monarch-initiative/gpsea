import hpotk
import typing
import logging

from hpotk.validate import ValidationLevel

from genophenocorr.model import Phenotype

from ._audit import Auditor, AuditReport, DataSanityIssue, Level


class PhenotypeCreator(Auditor[typing.Iterable[typing.Tuple[str, bool]], typing.Sequence[Phenotype]]):
    """
    `PhenotypeCreator` validates the input phenotype features and prepares them for the downstream analysis.

    The creator expects an iterable with tuples that contain a CURIE and status. The CURIE must correspond
    to a HPO term identifier and status must be a `bool`.

    The creator prunes CURIES with simple errors such as malformed CURIE or non-HPO terms and validates the rest
    with HPO toolkit's validator.

    The results are wrapped into :class:`AuditReport`.
    """

    def __init__(self, hpo: hpotk.MinimalOntology,
                 validator: hpotk.validate.ValidationRunner):
        self._logger = logging.getLogger(__name__)
        self._hpo = hpotk.util.validate_instance(hpo, hpotk.MinimalOntology, 'hpo')
        self._validator = hpotk.util.validate_instance(validator, hpotk.validate.ValidationRunner, 'validator')

    def process(self, inputs: typing.Iterable[typing.Tuple[str, bool]]) -> AuditReport[typing.Sequence[Phenotype]]:
        """Creates a list of Phenotype objects from term IDs and checks if the term IDs satisfy the validation requirements.

        Args:
            inputs (Iterable[Tuple[str, bool]]): A list of Tuples, structured (HPO IDs, boolean- True if observed)
        Returns:
            A sequence of Phenotype objects
        Error:
            PhenotypeValidationException: An instance of an issue with the ValidationRunner
        """
        issues = []
        phenotypes = []

        for curie, is_observed in inputs:
            # Check the CURIE is well-formed
            try:
                term_id = hpotk.TermId.from_curie(curie)
            except ValueError as ve:
                issues.append(DataSanityIssue(Level.WARN, ve.args[0],
                                              'Ensure the term ID consists of a prefix (e.g. `HP`) '
                                              'and id (e.g. `0001250`) joined by colon `:` or underscore `_`')
                              )
                continue

            # Check the term is an HPO concept
            if term_id.prefix != 'HP':
                issues.append(
                    DataSanityIssue(Level.WARN,
                                    f'{term_id} is not an HPO term',
                                    'Remove non-HPO concepts from the analysis input')
                )
                continue

            # Term must be present in HPO
            term = self._hpo.get_term(term_id)
            if term is None:
                issues.append(
                    DataSanityIssue(Level.WARN,
                                    f'{term_id} is not in HPO version `{self._hpo.version}`',
                                    'Correct the HPO term or use the latest HPO for the analysis')
                )
                continue

            if term.identifier != term_id:
                # Input includes an obsolete term ID. We emit a warning and update the term ID behind the scenes,
                # since `term.identifier` always returns the primary term ID.
                issues.append(
                    DataSanityIssue(Level.WARN,
                                    f'{term_id} is an obsolete identifier for {term.name}',
                                    f'Replace {term_id} with the primary term ID {term.identifier}')
                )

            phenotypes.append(Phenotype.from_term(term, is_observed))

        # Check we have some phenotype terms to work with.
        if len(phenotypes) == 0:
            issues.append(
                DataSanityIssue(
                    Level.WARN,
                    f'No phenotype terms were left after the validation',
                    'Revise the phenotype terms and try again')
            )
        else:
            vr = self._validator.validate_all(phenotypes)
            for result in vr.results:
                level = self._translate_level(result.level)
                if level is None:
                    # Should not happen. Please let the developers know about this issue!
                    raise ValueError(f'Unknown result validation level {result.level}')

                issues.append(
                    DataSanityIssue(level,
                                    result.message,
                                    solution='Correct the input data')
                )

        return AuditReport(phenotypes, issues)

    @staticmethod
    def _translate_level(lvl: ValidationLevel) -> typing.Optional[Level]:
        if lvl == ValidationLevel.WARNING:
            return Level.WARN
        elif lvl == ValidationLevel.ERROR:
            return Level.ERROR
        else:
            return None
