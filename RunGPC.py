from collections import defaultdict
import glob
from gpc import Patient, Disease, AllPatients
import os
from gpc import compare
from gpc import genotype_class as GC

allPatients = AllPatients()

for file in glob.glob('phenopackets/AutismTest/*.json'):
    fileName = os.path.basename(file)
    current = Patient(file)
    allPatients.add(current)



#print(allPatients.list_all_diseases())
#print(allPatients.list_all_phenotypes())
    
#for p in allPatients: print(allPatients[p].describe())

print(compare.RunStats(allPatients, GC.is_deletion, GC.is_duplication))