
# TODO ADJUST LIKE Phenotype. Check that disease id is OMIM, MONDO, or ORPHA, or DECIPHER
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