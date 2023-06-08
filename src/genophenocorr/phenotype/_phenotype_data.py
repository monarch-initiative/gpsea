import hpotk
import typing


class Phenotype(hpotk.model.Identified):

    def __init__(self, term_id: typing.Optional[hpotk.model.Identified, hpotk.model.TermId],
                 observed: typing.Optional[bool]) -> None:
        if isinstance(term_id, hpotk.model.Identified):
            # Covers Term, MinimalTerm, well, anything that has a `TermId`.
            self._term_id = term_id.identifier
        elif isinstance(term_id, hpotk.model.TermId):
            self._term_id = term_id
        else:
            raise ValueError(f"`term_id` must be an instance of 'TermId' or `Identified` but it was {type(term_id)}")
        if observed is not None:
            if not isinstance(observed, bool):
                raise ValueError(f"observed variable must be type boolean, but is type {type(observed)}")
        self._observed = observed

    @property
    def identifier(self) -> hpotk.model.TermId:
        return self._term_id

    @property
    def observed(self) -> typing.Optional[bool]:
        return self._observed

    def __eq__(self, other):
        return isinstance(other, Phenotype) \
            and self.identifier == other.identifier \
            and self.observed == other.observed 

    def __hash__(self):
        return hash((self._term_id, self._observed))

    def __str__(self):
        return f"Phenotype(identifier={self._term_id}, " \
               f"observed={self._observed})"

    def __repr__(self):
        return str(self)


