.. _predicates:

==========
Predicates
==========

Searching for genotype-phenotype associations usually requires that
the individuals are partitioned into two or more discrete groups to allow testing for the inter-group differences.
GPSEA reflects these requirements with its predicate API.


A **predicate** must be capable of partitioning the individuals into two or more groups.
The groups must be *exclusive* - each individual must be assigned at most into one group.
In general, it is desirable that the groups cover all or at least the majority of the cohort being analyzed to maximize statistical power.
However, the predicate is allowed to return `None` if the individual cannot be assigned.
As a result, the individual will be omitted from the downstream analysis.

Predicates serve both the *genotype* and *phenotype* prongs of the analysis.
Genotype predicates (:class:`~gpsea.analysis.predicate.genotype.GenotypePolyPredicate`)
assign the :class:`~gpsea.model.Patient`
into a group (mostly) based on the variant information, while the
phenotype predicates (:class:`~gpsea.analysis.predicate.phenotype.PhenotypePolyPredicate`)
use the HPO terms to assign a group.

All GPSEA analyses need at least one predicate (typically a *genotype* predicate) and many require both *genotype* and *phenotype* predicates.
The following pages provide more information.



.. toctree::
  :maxdepth: 1
  :caption: Contents:

  predicates/phenotype_predicates
  predicates/genotype_predicates


*******
Gallery
*******

Here we show examples of predicates used in some of our analyses.

TODO


**********
Need more?
**********

Please see :class:`~gpsea.analysis.predicate.genotype.VariantPredicates` 
and :class:`~gpsea.analysis.predicate.genotype.ProteinPredicates` 
for a list of the predicates available off the shelf.

However, feel free to open an issue on our `GitHub tracker <https://github.com/monarch-initiative/gpsea/issues>`_
if a predicate seems to be missing.
