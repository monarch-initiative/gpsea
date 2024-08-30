import typing
import warnings

import hpotk


class Phenotype(hpotk.model.Identified, hpotk.model.ObservableFeature):
    """
    `Phenotype` represents a clinical sign or symptom represented as an HPO term.
    
    The phenotype can be either present in the patient or excluded.
    """

    @staticmethod
    def from_term(term: hpotk.model.MinimalTerm, is_observed: bool):
        return Phenotype.from_raw_parts(term.identifier, is_observed)

    @staticmethod
    def from_raw_parts(
        term_id: typing.Union[str, hpotk.TermId],
        is_observed: bool,
    ) -> "Phenotype":
        """
        Create `Phenotype` from a term ID and observation state.

        :param term_id: a `str` with CURIE (e.g. `HP:0001250`) or a :class:`~hpotk.TermId`.
        :param is_observed: `True` if the term ID was observed in patient or `False` if it was explicitly excluded.
        """
        if isinstance(term_id, str):
            term_id = hpotk.TermId.from_curie(term_id)
        elif isinstance(term_id, hpotk.TermId):
            pass
        else:
            raise ValueError('`term_id` must be either a `str` or a `hpotk.TermId`')
        
        return Phenotype(
            term_id,
            is_observed,
        )

    def __init__(
        self,
        term_id: hpotk.TermId,
        is_observed: bool
    ):
        self._term_id = hpotk.util.validate_instance(term_id, hpotk.TermId, 'term_id')
        self._observed = hpotk.util.validate_instance(is_observed, bool, 'is_observed')

    @property
    def identifier(self) -> hpotk.TermId:
        """
        Get the phenotype ID; an HPO term ID most of the time.
        """
        return self._term_id

    @property
    def is_present(self) -> bool:
        """
        Return `True` if the phenotype feature was observed in the subject or `False` if the feature's presence was explicitly excluded.
        """
        return self._observed

    @property
    def observed(self) -> typing.Optional[bool]:
        """Returns a boolean for whether the phenotype is observed.

        Returns:
            boolean: True if this phenotype was observed in the respective patient.
        """
        warnings.warn('`observed` property was deprecated and will be removed in `v0.3.0`. '
                      'Use `is_present` instead', DeprecationWarning, stacklevel=2)
        return self.is_present

    @property
    def is_observed(self) -> bool:
        """
        Returns `True` if the phenotype was *present* in the respective patient.
        """
        warnings.warn('`is_observed` property was deprecated and will be removed in `v0.3.0`. '
                      'Use `is_present` instead', DeprecationWarning, stacklevel=2)
        return self.is_present

    def __eq__(self, other):
        return isinstance(other, Phenotype) \
            and self._term_id == other._term_id \
            and self._observed == other._observed

    def __hash__(self):
        return hash((self._term_id, self._observed))

    def __str__(self):
        return f"Phenotype(" \
               f"identifier={self._term_id}, " \
               f"is_present={self._observed})"

    def __repr__(self):
        return str(self)


class Disease(hpotk.model.Identified, hpotk.model.ObservableFeature, hpotk.model.Named):
    """
    Representation of a disease diagnosed (or excluded) in an investigated individual.
    """

    def __init__(
        self,
        term_id: hpotk.TermId,
        name: str,
        is_observed: bool,
    ):
        self._term_id = hpotk.util.validate_instance(term_id, hpotk.TermId, 'term_id')
        self._name = hpotk.util.validate_instance(name, str, 'name')
        self._observed = hpotk.util.validate_instance(is_observed, bool, 'is_observed')

    @property
    def identifier(self) -> hpotk.TermId:
        """
        Get the disease ID.
        """
        return self._term_id

    @property
    def name(self):
        """
        Get the disease label (e.g. `LEIGH SYNDROME, NUCLEAR; NULS`).
        """
        return self._name

    @property
    def is_present(self) -> bool:
        """
        Return `True` if the disease was diagnosed in the individual or `False` if it was excluded.
        """
        return self._observed

    def __eq__(self, other):
        return isinstance(other, Disease) \
            and self._term_id == other._term_id \
            and self._name == other._name \
            and self._observed == other._observed

    def __hash__(self):
        return hash((self._term_id, self._name, self._observed))

    def __str__(self):
        return f"Disease(" \
               f"identifier={self._term_id}, " \
               f"name={self._name}, " \
               f"is_present={self._observed})"

    def __repr__(self):
        return str(self)
