from collections import defaultdict
import glob
from gpc import Patient, Disease, AllPatients
import os
from gpc import compare
from gpc import compare_func as GC

allPatients = AllPatients()

for file in glob.glob('phenopackets/PlatzerData/*.json'):
    fileName = os.path.basename(file)
    current = Patient(file)
    allPatients.add(current)

allPatients.count_vars_per_feature()

for prot in allPatients.all_proteins.values():
    feat = prot.features

for typ in allPatients.all_var_types:
    print(f'Comparing all in {typ} vs out of {typ}')
    compare.RunStats(allPatients, GC.is_var_type, GC.is_not_var_type, typ, typ)
    for typ2 in allPatients.all_var_types:
        print(f'Comparing all in {typ} vs all in {typ2}')
        compare.RunStats(allPatients, GC.is_var_type, GC.is_var_type, typ, typ2)

for var in allPatients.all_variants:
    print(f'Comparing affecting {var.variant_string} vs all not')
    compare.RunStats(allPatients, GC.is_var_match, GC.is_not_var_match, var, var)

compare.RunStats(allPatients,GC.in_feature, GC.not_in_feature, [feat, 'Domain'], [feat, 'Domain'] )