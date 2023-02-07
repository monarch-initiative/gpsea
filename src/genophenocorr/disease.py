

class Disease:
    def __init__(self, id, label):
        #if any(x in id for x in ['OMIM','MONDO', 'ORPHA', 'DECIPHER']):
        self._diseaseID = id
        self._diseaseLabel = label
        #else:
        #    self._diseaseID = None
        #    self._diseaseLabel = None
        
    @property
    def id(self):
        return self._diseaseID
    
    @property
    def label(self):
        return self._diseaseLabel