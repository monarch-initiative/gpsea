

class Disease:
    def __init__(self, id, label):
        if id is not None or len(id) == 0:
            if any(x in id for x in ['OMIM','MONDO', 'ORPHA', 'DECIPHER']):
                self._diseaseID = id
                self._diseaseLabel = label
            else:
                raise ValueError("Disease ID must include 'OMIM','MONDO', 'ORPHA', or 'DECIPHER'.")
        else:
            self._diseaseID = None
            self._diseaseLabel = None
        
    @property
    def id(self):
        return self._diseaseID
    
    @property
    def label(self):
        return self._diseaseLabel