from collections import defaultdict
import glob
from gpc import Patient, Disease, AllPatients
import os
from gpc import compare
from gpc import compare_func as GC

allPatients = AllPatients()

for file in glob.glob('phenopackets/cohort/*.json'):
    fileName = os.path.basename(file)
    current = Patient(file)
    allPatients.add(current)


byDisease = allPatients.split_by_disease()
smallPatientList = byDisease['OMIM:616975']
#print(allPatients.list_all_diseases())
#print(allPatients.list_all_phenotypes())
    
#for p in allPatients: print(allPatients[p].describe())

print(compare.RunStats(smallPatientList, GC.in_domain, GC.in_domain, [103, 283], [284,387]))