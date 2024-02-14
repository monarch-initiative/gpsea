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
from genophenocorr.analysis import configure_cohort_analysis
from genophenocorr.analysis.predicate import BooleanPredicate
from genophenocorr.model import VariantEffect

cohort_analysis = configure_cohort_analysis(cohort, hpo)
frameshift = cohort_analysis.compare_by_variant_effect(VariantEffect.FRAMESHIFT_VARIANT, tx_id='NM_1234.5')

frameshift.summarize(hpo, phenotype_category=BooleanPredicate.YES)
```

provides a pandas data frame with genotype-phenotype correlations:

```text
FRAMESHIFT_VARIANT on NM_1234.5        No             Yes
                                    Count   Percent Count Percent   p value Corrected p value
    Arachnodactyly [HP:0001166]      1/26      3.84 13/26    50.0   0.00078          0.020299
    Spasticity [HP:0001257]          6/26     23.08 11/26    42.3   0.69245          1.000000
    Hypertonia [HP:0001276]          6/17     35.29 11/17    64.7   1.00000          1.000000
    ...                               ...       ...    ...    ...       ...               ...
```

## Documentation

Check out the User guide and the API reference for more info:

- [Stable documentation](https://monarch-initiative.github.io/genophenocorr/stable/) (last release on `main` branch)
- [Latest documentation](https://monarch-initiative.github.io/genophenocorr/latest) (bleeding edge, latest commit on `development` branch)
