import hpotk
import typing


class Phenotype(hpotk.model.Identified, hpotk.model.Named):
    """A class that represents an HPO verified phenotype

    Attributes:
        identifier (hpotk.model.Named): The HPO ID associated with this phenotype
        name (str): The official HPO name for this phenotype
        observed (bool): Is True if this phenotype was observed in the respective patient
    """

    @staticmethod
    def from_term(term: hpotk.model.MinimalTerm, observed: typing.Optional[bool]):
        return Phenotype(term.identifier, term.name, observed)

    def __init__(self, term_id: hpotk.TermId,
                 name: str,
                 observed: typing.Optional[bool]) -> None:
        """Constructs all necessary attributes for a Phenotype object

        Args:
            term (hpotk.model.MinimalTerm): An hpotk object that has an HPO ID and Name
            observed (bool): Is True if this phenotype was observed in the respective patient
        """
        self._term_id = hpotk.util.validate_instance(term_id, hpotk.TermId, 'term_id')
        self._name = hpotk.util.validate_instance(name, str, 'name')

        if observed is not None:
            if not isinstance(observed, bool):
                raise ValueError(f"observed variable must be type boolean, but is type {type(observed)}")
        self._observed = observed

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
    def observed(self) -> typing.Optional[bool]:
        """Returns a boolean for whether the phenotype is observed.

        Returns:
            boolean: True if this phenotype was observed in the respective patient.
        """
        return self._observed

    def __eq__(self, other):
        return isinstance(other, Phenotype) \
            and self.identifier == other.identifier \
            and self.name == other.name \
            and self.observed == other.observed

    def __hash__(self):
        return hash((self.identifier, self.name, self.observed))

    def __str__(self):
        return f"Phenotype(" \
               "identifier={self.identifier}, " \
               f"name={self.name}, " \
               f"observed={self._observed})"

    def __repr__(self):
        return str(self)


