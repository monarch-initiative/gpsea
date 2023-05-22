import hpotk
import typing


class Phenotype(hpotk.model.Identified):

    def __init__(self, term: hpotk.model.Term, observed: typing.Optional[bool]) -> None:
        
        if not isinstance(term, hpotk.model.Term) and not isinstance(term, hpotk.model._term.DefaultMinimalTerm):
            raise ValueError(f"term must be type 'Term' not type {type(term)}")
        self._term = term
        if observed is not None:
            if not isinstance(observed, bool):
                raise ValueError(f"observed variable must be type boolean, but is type {type(observed)}")
        self._observed = observed

    @property
    def identifier(self) -> hpotk.model.TermId:
        return self._term.identifier

    @property
    def observed(self):
        return self._observed

    def __eq__(self, other):
        return isinstance(other, Phenotype) \
            and self.identifier == other.identifier \
            and self.observed == other.observed 

    def __hash__(self):
        return hash((self.identifier, self.observed))

    def __str__(self):
        return f"Phenotype(identifier={self.identifier}, " \
               f"observed={self.observed})"

    def __repr__(self):
        return str(self)


