import logging
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


class DiseaseCreator:

    def __init__(self):
        self._logger = logging.getLogger(__name__)

    def create_disease(self, id, label):
        term_id = hpotk.model.TermId.from_curie(id)
        if term_id.prefix not in ['OMIM','MONDO', 'ORPHA', 'DECIPHER']:
            raise ValueError(f'Disease ID must be an \'OMIM\',\'MONDO\', \'ORPHA\', or \'DECIPHER\' ID but was {id}')
        return Disease(term_id, label)