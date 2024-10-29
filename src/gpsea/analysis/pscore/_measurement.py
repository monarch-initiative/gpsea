import typing
import hpotk

from gpsea.model import Patient

from ._api import PhenotypeScorer


class MeasurementPhenotypeScorer(PhenotypeScorer):
    """
    `MeasurementPhenotypeScorer` uses a value of a measurement as a phenotype score.

    For instance, the amount of `Testosterone [Mass/volume] in Serum or Plasma <https://loinc.org/2986-8/>`_.

    
    Example
    ^^^^^^^

    Create a scorer that uses the level of testosterone represented by the
    `Testosterone [Mass/volume] in Serum or Plasma <https://loinc.org/2986-8/>`_
    LOINC code as a phenotype score.

    >>> from gpsea.analysis.pscore import MeasurementPhenotypeScorer
    >>> pheno_scorer = MeasurementPhenotypeScorer.from_measurement_id(
    ...     term_id="LOINC:2986-8",
    ...     label="Testosterone [Mass/volume] in Serum or Plasma",
    ... )
    >>> # use the scorer in the analysis ...
    """

    @staticmethod
    def from_measurement_id(
        term_id: typing.Union[str, hpotk.TermId],
        label: str,
    ) -> "MeasurementPhenotypeScorer":
        """
        Create `MeasurementPhenotypeScorer` from a measurement identifier.

        :param term_id: a `str` with CURIE or a :class:`~hpotk.TermId`
            representing the term ID of a measurement (e.g. `LOINC:2986-8`).
        :param label: a `str` with the measurement label (e.g. `Testosterone [Mass/volume] in Serum or Plasma`)
        """
        return MeasurementPhenotypeScorer(
            term_id=term_id,
            label=label,
        )

    def __init__(
        self,
        term_id: typing.Union[str, hpotk.TermId],
        label: str,
    ):
        if isinstance(term_id, str):
            self._identifier = hpotk.TermId.from_curie(term_id)
        elif isinstance(term_id, hpotk.TermId):
            self._identifier = term_id
        else:
            raise ValueError(
                f"`term_id` must be a `str` or a `hpotk.TermId` but was {type(term_id)}"
            )
        
        assert isinstance(label, str)
        self._label = label

    @property
    def name(self) -> str:
        return "Measurement"

    @property
    def description(self) -> str:
        return f"Value of {self._label} [{self._identifier.value}]"

    @property
    def variable_name(self) -> str:
        return self._identifier.value

    @property
    def term_id(self) -> hpotk.TermId:
        return self._identifier
    
    @property
    def label(self) -> str:
        return self._label

    def score(
        self,
        patient: Patient,
    ) -> float:
        """
        Compute the phenotype score.
        """
        m = patient.measurement_by_id(self._identifier)
        if m is not None:
            return m.test_result

        return float("nan")

    def __eq__(self, value: object) -> bool:
        return (
            isinstance(value, MeasurementPhenotypeScorer)
            and self._identifier == value._identifier
            and self._label == value._label
        )

    def __hash__(self) -> int:
        return hash((self._identifier, self._label))

    def __str__(self) -> str:
        return f"MeasurementPhenotypeScorer(term_id={self._identifier}, label={self._label})"

    def __repr__(self) -> str:
        return str(self)
