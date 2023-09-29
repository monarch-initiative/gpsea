[![Build status](https://github.com/monarch-initiative/genophenocorr/workflows/CI/badge.svg)](https://github.com/monarch-initiative/genophenocorr/actions/workflows/python_ci.yml)
[![GitHub release](https://img.shields.io/github/release/monarch-initiative/genophenocorr.svg)](https://github.com/monarch-initiative/genophenocorr/releases)
![PyPi downloads](https://img.shields.io/pypi/dm/genophenocorr.svg?label=Pypi%20downloads)
![PyPI - Python Version](https://img.shields.io/pypi/pyversions/genophenocorr)

Genophenocorr is a Python library for genotype-phenotype association analysis. 

An example of simple genotype-phenotype association analysis
```python
# Load HPO
import hpotk
hpo = hpotk.load_minimal_ontology('http://purl.obolibrary.org/obo/hp.json')

# Load a cohort of phenopackets 
from genophenocorr.data import get_toy_cohort
cohort = get_toy_cohort()

# Analyze genotype-phenotype associations 
from genophenocorr.analysis import CohortAnalysis
from genophenocorr.constants import VariantEffect

cohort_analysis = CohortAnalysis(cohort, 'NM_1234.5', hpo)
frameshift = cohort_analysis.compare_by_variant_type(VariantEffect.FRAMESHIFT_VARIANT)
print(frameshift)
```

prints a table with genotype-phenotype correlations:

```text
                            With frameshift_variant         Without frameshift_variant
                                              Count Percent                      Count Percent  p-value
HP:0001166 (Arachnodactyly)                       4  30.77%                         10  76.92%  0.04718
HP:0001250 (Seizure)                             11  84.62%                          9  69.23%  0.64472
HP:0001257 (Spasticity)                           8  61.54%                          9  69.23%  1.00000
```

## Documentation

Check out the User guide and the API reference for more info:

- [Stable documentation](https://monarch-initiative.github.io/genophenocorr/stable/) (last release on `main` branch)
- [Latest documentation](https://monarch-initiative.github.io/genophenocorr/latest) (bleeding edge, latest commit on `development` branch)
