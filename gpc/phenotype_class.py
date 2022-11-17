class Phenotype:
    def __init__(self, phenopacket):
        
        self._hpid = phenopacket.type.id
        self._label = phenopacket.type.label
        
    @property
    def id(self):
        return self._hpid
    
    @property
    def label(self):
        return self._label