


class Phenotype:
    def __init__(self, id, label, excluded=False, measured=True) -> None:
        if not id.startswith("HP:") and len(id) == 10:
            raise ValueError(f"Malformed HPO id: {id}") 
        self._hpid = id
        self._label = label
        self._excluded = excluded
        self._measured = measured
        
    @property
    def id(self):
        return self._hpid
    
    @property
    def label(self):
        return self._label