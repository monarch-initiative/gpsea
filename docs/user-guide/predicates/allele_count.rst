.. _allele-count:

=====================
Group by allele count
=====================

Sometimes, we may want to compare individuals with different allele counts of a single variant category.
For instance, we may want to compare the survival of individuals harboring a mutation in *EGFR* (:math:`AC = 1`)
with those with no such mutation (:math:`AC = 0`). 
Alternatively, in some genes, heterozygous mutations (:math:`AC = 1`) and biallelic mutations (:math:`AC = 2`)
can lead to different diseases.

The allele count analysis differs from the :ref:`variant-category` analysis.
The allele count analysis partitions the individuals based on different allele counts of one genotype category,
while the variant category analysis partitions the individuals based on a fixed allele count of different genotype categories.

The allele count analysis can be done with :func:`~gpsea.analysis.predicate.genotype.allele_count` predicate.

********
Examples
********

Compare the individuals with *EGFR* mutation
============================================

First, let's create a :class:`~gpsea.analysis.predicate.genotype.VariantPredicate` to include
any *EGFR* variant:

>>> from gpsea.analysis.predicate.genotype import VariantPredicates
>>> affects_egfr = VariantPredicates.gene(symbol="EGFR")
>>> affects_egfr.get_question()
'affects EGFR'

Next, we create allele count predicate to partition the individuals
based on presence of zero or one *EGFR* mutation allele:

>>> from gpsea.analysis.predicate.genotype import allele_count
>>> gt_predicate = allele_count(
...     counts=({0,}, {1,}),
...     target=affects_egfr,
... )
>>> gt_predicate.display_question()
'Allele count: 0, 1'

We create the predicate with two arguments.
The `counts` argument takes a tuple of two sets, to partition the individuals
based on zero (``{0,}``) or one (``{1,}``) target variant allele.
The `target` takes a :class:`~gpsea.analysis.predicate.genotype.VariantPredicate`
for defining the target variants.

We can use the `gt_predicate` to partition a cohort along the genotype axis,
e.g. to compare the patient survivals in a `survival analysis <survival>`.


Compare the individuals with monoallelic and biallelic mutations
================================================================

Let's prepare a predicate for grouping individuals based on one or two alleles of a target mutation.

For this example, the target mutation is any mutation that affects *LMNA*

>>> from gpsea.analysis.predicate.genotype import VariantPredicates
>>> affects_lmna = VariantPredicates.gene(symbol="LMNA")
>>> affects_lmna.get_question()
'affects LMNA'

and we will compare the individuals with one allele with those with two alleles:

>>> gt_predicate = allele_count(
...     counts=({1,}, {2,}),
...     target=affects_lmna,
... )
>>> gt_predicate.display_question()
'Allele count: 1, 2'

The predicate will partition the individuals into two groups:
those with one *LMNA* variant allele and those with two *LMNA* variant alleles.
The individual with other allele counts (e.g. `0` or `3`) will be excluded
from the analysis.

