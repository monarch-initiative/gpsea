import hpotk


class Disease(hpotk.model.Identified, hpotk.model.Named):

    def __init__(self, term_id: hpotk.model.TermId, label: str):
        # TODO(lnrekerle) - instance/type check
        self._term_id = term_id
        self._label = label

        # TODO(lnrekerle) - move the checks to a class that loads and/or prepares the diseases & samples
        # if id is not None or len(id) == 0:
        #     #if any(x in id for x in ['OMIM','MONDO', 'ORPHA', 'DECIPHER']):
        #     self._diseaseID = id
        #     self._diseaseLabel = label
        #     #else:
        #     #    raise ValueError("Disease ID must include 'OMIM','MONDO', 'ORPHA', or 'DECIPHER'.")
        # else:
        #     self._diseaseID = None
        #     self._diseaseLabel = None

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

    def __str__(self):
        return f"Disease(identifier={self.identifier}, " \
               f"name={self.name})"

    def __repr__(self):
        return str(self)
