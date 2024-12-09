.. _partitioning:

############
Partitioning
############

Searching for G/P associations usually requires to partition the individuals
into two or more discrete groups, to allow testing for the inter-group differences.
GPSEA reflects these requirements with its predicate API.


A *predicate* partitions an individual (:class:`~gpsea.model.Patient`) into one of several groups.
The groups must be *exclusive* - each individual must be assigned at most into one group.
In general, it is desirable that the groups cover all or at least the majority of the cohort being analyzed to maximize statistical power.
However, the predicate is allowed to return `None` if the individual cannot be assigned.
As a result, the individual will be omitted from the downstream analysis.

Predicates can be applied on both *genotype* and *phenotype*.
Genotype predicates (:class:`~gpsea.analysis.predicate.genotype.GenotypePolyPredicate`)
assign the individual into a group (mostly) based on the variant information,
while the phenotype predicates (:class:`~gpsea.analysis.predicate.phenotype.PhenotypePolyPredicate`)
decide on the group based on the HPO terms or a diagnosis.

Besides assigning an individual into a discrete group,
a continuous score (typically for phenotype) can also be computed.
A :class:`~gpsea.analysis.pscore.PhenotypeScorer` computes a phenotype score
for an individual, and the differences in score distributions of genotype groups
can be tested e.g. with Mann-Whitney U test.

Any GPSEA analysis needs a genotype predicate and a phenotype predicate/scorer.
The following sections show how to prepare the predicates and scorers
for your analyses. We also show a gallery with examples.


.. toctree::
  :maxdepth: 1
  :caption: Contents:

  genotype/index
  phenotype/index
  gallery
