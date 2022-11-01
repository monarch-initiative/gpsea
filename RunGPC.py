from collections import defaultdict
import glob
from gpc import Patient
import os

allPatients = defaultdict(Patient)

for file in glob.glob('phenopackets/*.json'):
    fileName = os.path.basename(file)
    print(fileName)
    current = Patient(file)
    
    if current.get_genotypes is not None and len(current.get_genotypes) != 0:
        allPatients[fileName] = current
    
for p in allPatients: print(allPatients[p].describe())