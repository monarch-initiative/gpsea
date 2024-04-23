import typing
import warnings

import hpotk


class Phenotype(hpotk.model.Identified, hpotk.model.ObservableFeature, hpotk.model.Named):
    """A class that represents an HPO verified phenotype

    Attributes:
        term_id (hpotk.model.Named): The HPO ID associated with this phenotype
        name (str): The official HPO name for this phenotype
        is_observed (bool): Is True if this phenotype was observed in the respective patient
    """

    @staticmethod
    def from_term(term: hpotk.model.MinimalTerm, is_observed: bool):
        return Phenotype(term.identifier, term.name, is_observed)

    def __init__(self, term_id: hpotk.TermId,
                 name: str,
                 is_observed: bool) -> None:
        self._term_id = hpotk.util.validate_instance(term_id, hpotk.TermId, 'term_id')
        self._name = hpotk.util.validate_instance(name, str, 'name')
        self._observed = hpotk.util.validate_instance(is_observed, bool, 'is_observed')

    @property
    def identifier(self) -> hpotk.TermId:
        """Returns an HPO ID unique to this Phenotype object.
        You can find more information about the phenotype by searching this ID at https://hpo.jax.org/app/

        Returns:
            hpotk.model.Named: HPO ID
        """
        return self._term_id

    @property
    def name(self):
        """Returns a string that describes this Phenotype object.

        Returns:
            string: phenotype name
        """
        return self._name

    @property
    def is_present(self) -> bool:
        """
        Returns:
            boolean: `True` if the phenotype feature was observed in the subject or `False` if the feature's presence
            was explicitly excluded.
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
            and self.identifier == other.identifier \
            and self.name == other.name \
            and self.is_present == other.is_present

    def __hash__(self):
        return hash((self.identifier, self.name, self.is_present))

    def __str__(self):
        return f"Phenotype(" \
               f"identifier={self.identifier}, " \
               f"name={self.name}, " \
               f"is_present={self._observed})"

    def __repr__(self):
        return str(self)


class Disease(hpotk.model.Identified, hpotk.model.ObservableFeature, hpotk.model.Named):
    """A class that represents a disease

    Attributes:
        term_id (hpotk.model.Named): The ID given by the user to reference this disease
        name (str): The name given by the user for this disease
        is_observed (bool): Is True if this disease was observed in the respective patient
    """

    def __init__(self, term_id: hpotk.TermId,
                 name: str,
                 is_observed: bool) -> None:
        self._term_id = hpotk.util.validate_instance(term_id, hpotk.TermId, 'term_id')
        self._name = hpotk.util.validate_instance(name, str, 'name')
        self._observed = hpotk.util.validate_instance(is_observed, bool, 'is_observed')

    @property
    def identifier(self) -> hpotk.TermId:
        """Returns an ID unique to this Disease object.
        
        Returns:
            hpotk.model.Named: Disease ID
        """
        return self._term_id

    @property
    def name(self):
        """Returns a string that describes this Disease object.

        Returns:
            string: disease name
        """
        return self._name

    @property
    def is_present(self) -> bool:
        """
        Returns:
            boolean: `True` if the disease was observed in the subject or `False` if the disease's presence
            was explicitly excluded.
        """
        return self._observed

    def __eq__(self, other):
        return isinstance(other, Disease) \
            and self.identifier == other.identifier \
            and self.name == other.name \
            and self.is_present == other.is_present

    def __hash__(self):
        return hash((self.identifier, self.name, self.is_present))

    def __str__(self):
        return f"Disease(" \
               f"identifier={self.identifier}, " \
               f"name={self.name}, " \
               f"is_present={self._observed})"

    def __repr__(self):
        return str(self)