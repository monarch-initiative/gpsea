

class Disease:
    def __init__(self, phenopacket):
        self._diseaseID = phenopacket.term.id
        self._diseaseLabel = phenopacket.term.label
        
    @property
    def id(self):
        return self._diseaseID
    
    @property
    def label(self):
        return self._diseaseLabel
    
    @property
    def is_duplication(self):
        if 'duplication' in self._diseaseLabel:
            return True
        else:
            return False
    
    @property
    def is_deletion(self):
        if 'deletion' in self._diseaseLabel:
            return True
        else:
            return False