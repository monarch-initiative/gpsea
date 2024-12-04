import abc
import typing
import warnings

import hpotk

from ._temporal import Age


class OnsetAware(metaclass=abc.ABCMeta):
    """
    A mixin for classes which may know about their :class:`~gpsea.model.Age` or onset.

    For instance, the onset of a phenotype or a disease in an individual.
    """

    @property
    @abc.abstractmethod
    def onset(self) -> typing.Optional[Age]:
        """
        Get the age of onset or `None` if unknown or otherwise not available.
        """
        pass


class Phenotype(hpotk.model.Identified, hpotk.model.ObservableFeature, OnsetAware):
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
        onset: typing.Optional[Age] = None,
    ) -> "Phenotype":
        """
        Create `Phenotype` from a term ID and observation state.

        :param term_id: a `str` with CURIE (e.g. `HP:0001250`) or a :class:`~hpotk.TermId`.
        :param is_observed: `True` if the term ID was observed in patient or `False` if it was explicitly excluded.
        :param onset: the :class:`~gpsea.model.Age` when the phenotype was first observed in patient
            or `None` if not available.
        """
        if isinstance(term_id, str):
            term_id = hpotk.TermId.from_curie(term_id)
        elif isinstance(term_id, hpotk.TermId):
            pass
        else:
            raise ValueError('`term_id` must be either a `str` or a `hpotk.TermId`')
        
        return Phenotype(
            term_id,
            is_observed=is_observed,
            onset=onset,
        )

    def __init__(
        self,
        term_id: hpotk.TermId,
        is_observed: bool,
        onset: typing.Optional[Age],
    ):
        assert isinstance(term_id, hpotk.TermId)
        self._term_id = term_id
        assert isinstance(is_observed, bool)
        self._observed = is_observed
        if onset is not None:
            assert isinstance(onset, Age)
        self._onset = onset

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
            bool: `True` if this phenotype was observed in the respective patient.
        """
        warnings.warn('`observed` property was deprecated and will be removed in `v0.3.0`. '
                      'Use `is_present` instead', DeprecationWarning, stacklevel=2)
        return self.is_present

    @property
    def onset(self) -> typing.Optional[Age]:
        """
        Get the age of onset of the phenotype or `None` if unknown or otherwise not available.
        """
        return self._onset

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
            and self._observed == other._observed \
            and self._onset == other._onset

    def __hash__(self):
        return hash((self._term_id, self._observed, self._onset))

    def __str__(self):
        return f"Phenotype(" \
               f"identifier={self._term_id}, " \
               f"is_present={self._observed}, " \
               f"onset={self._onset})"

    def __repr__(self):
        return str(self)


class Disease(hpotk.model.Identified, hpotk.model.ObservableFeature, hpotk.model.Named, OnsetAware):
    """
    Representation of a disease diagnosed (or excluded) in an investigated individual.
    """

    @staticmethod
    def from_raw_parts(
        term_id: typing.Union[str, hpotk.TermId],
        name: str,
        is_observed: bool,
        onset: typing.Optional[Age] = None,
    ) -> "Disease":
        if isinstance(term_id, str):
            term_id = hpotk.TermId.from_curie(term_id)
        
        return Disease(
            term_id=term_id,
            name=name,
            is_observed=is_observed,
            onset=onset,
        )

    def __init__(
        self,
        term_id: hpotk.TermId,
        name: str,
        is_observed: bool,
        onset: typing.Optional[Age],
    ):
        assert isinstance(term_id, hpotk.TermId)
        self._term_id = term_id
        assert isinstance(name, str)
        self._name = name
        assert isinstance(is_observed, bool)
        self._observed = is_observed
        if onset is not None:
            assert isinstance(onset, Age)
        self._onset = onset

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
    
    @property
    def onset(self) -> typing.Optional[Age]:
        """
        Get the age of onset of the disease or `None` if unknown or otherwise not available.
        """
        return self._onset

    def __eq__(self, other):
        return isinstance(other, Disease) \
            and self._term_id == other._term_id \
            and self._name == other._name \
            and self._observed == other._observed \
            and self._onset == other._onset

    def __hash__(self):
        return hash((self._term_id, self._name, self._observed, self._onset))

    def __str__(self):
        return f"Disease(" \
               f"identifier={self._term_id}, " \
               f"name={self._name}, " \
               f"is_present={self._observed}, " \
               f"onset={self._onset})"

    def __repr__(self):
        return str(self)
    

class Measurement(hpotk.model.Identified, hpotk.model.Named):
    """
    Representation of a GA4GH Phenopacket Measurement (numerical test result).
    An intended use case would be to perform a Student's t test on numerical measurements in individuals
    with two difference genotype classes.
    """

    def __init__(
        self,
        test_term_id: hpotk.TermId,
        test_name: str,
        test_result: float,
        unit: hpotk.TermId,
    ):
        assert isinstance(test_term_id, hpotk.TermId)
        self._term_id = test_term_id
        assert isinstance(test_name, str)
        self._name = test_name
        assert isinstance(test_result, float)
        self._test_result = test_result
        assert isinstance(unit, hpotk.TermId)
        self._unit = unit

    @property
    def identifier(self) -> hpotk.TermId:
        """
        Get the test ID, e.g. ``LOINC:2986-8``.
        """
        return self._term_id

    @property
    def name(self):
        """
        Get the test label (e.g. *Testosterone [Mass/volume] in Serum or Plasma* for ``LOINC 2986-8``).
        """
        return self._name

    @property
    def test_result(self) -> float:
        """
        Return `True` if the disease was diagnosed in the individual or `False` if it was excluded.
        """
        return self._test_result
    
    @property
    def unit(self) -> hpotk.TermId:
        """
        Return the unit for the test result encoded as a TermId, e.g., ``UCUM:ng/dL`` for *nanogram per deciliter*.
        """
        return self._unit

    def __eq__(self, other):
        return isinstance(other, Measurement) \
            and self._term_id == other._term_id \
            and self._name == other._name \
            and self._test_result == other._test_result \
            and self._unit == other._unit

    def __hash__(self):
        return hash((self._term_id, self._name, self._test_result, self._unit))

    def __str__(self):
        return "Measurement(" \
               f"identifier={self._term_id}, " \
               f"name={self._name}, " \
               f"test_result={self._test_result}), " \
               f"unit={self._unit})"

    def __repr__(self):
        return str(self)
