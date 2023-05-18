import sys
import typing

import hpotk


class Disease(hpotk.model.Identified, hpotk.model.Named):

    def __init__(self, term_id: hpotk.model.TermId, label: str):
        if not isinstance(term_id, hpotk.model.TermId):
            raise ValueError(f'term_id must be type TermID but was type {type(term_id)}')
        self._term_id = term_id
        if not isinstance(label, str):
            raise ValueError(f'label must be type string but was type {type(label)}')
        self._label = label

    @property
    def identifier(self) -> hpotk.model.TermId:
        return self._term_id

    @property
    def name(self) -> str:
        return self._label

    def __eq__(self, other):
        return isinstance(other, Disease) \
            and self.identifier == other.identifier \
            and self.name == other.name

    def __hash__(self):
        return hash((self._term_id, self._label))

    def __str__(self):
        return f"Disease(identifier={self.identifier}, " \
               f"name={self.name})"

    def __repr__(self):
        return str(self)


# TODO(lnrekerle) - use the disease prefixes in the client code
# ['OMIM', 'MONDO', 'ORPHA', 'DECIPHER']
def create_disease(term_id_str: str,
                   label: str,
                   valid_prefixes: list[str] = ['OMIM', 'MONDO', 'ORPHA', 'DECIPHER']) -> typing.Optional[Disease]:
    """
    Create a `Disease` or raise an Exception if the `term_id` is invalid or if the `term_id` prefix
    is not present in `valid_prefixes`.
    """
    try:
        term_id = hpotk.model.TermId.from_curie(term_id_str)
    except ValueError as e:
        print(f'Cannot create a Disease: {e}', file=sys.stderr)
        return None

    if term_id.prefix not in valid_prefixes:
        print(f'Disease ID must be in {", ".join(valid_prefixes)} ID but was {term_id.prefix}', file=sys.stderr)
        return None
    return Disease(term_id, label)
