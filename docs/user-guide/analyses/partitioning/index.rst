.. _partitioning:

############
Partitioning
############

Searching for G/P associations usually requires to assign the individuals
into two or more discrete classes, to allow testing for the inter-class differences.
GPSEA reflects these requirements with its classifier API.


A *classifier* assigns an individual (:class:`~gpsea.model.Patient`) into one of several classes.
The classes must be *exclusive* - each individual must be assigned at most into one class.
In general, it is desirable that the classes cover all or at least the majority
of the cohort being analyzed to maximize statistical power.
However, the classifier is allowed to return `None` if the individual cannot be assigned.
As a result, the individual will be omitted from the downstream analysis.

Classifiers can be applied on both *genotype* and *phenotype*.
Genotype classifiers (:class:`~gpsea.analysis.clf.GenotypeClassifier`)
classify the individual (mostly) based on the variant information,
while the phenotype classifiers (:class:`~gpsea.analysis.clf.PhenotypeClassifier`)
decide on the class based on the HPO terms or a diagnosis.

Besides assigning an individual into a discrete phenotype class,
a continuous score or a survival can also be computed.
A :class:`~gpsea.analysis.pscore.PhenotypeScorer` computes a phenotype score
for an individual, and the differences in score distributions of genotype groups
can be tested e.g. with Mann-Whitney U test.
An :class:`~gpsea.analysis.temporal.Endpoint` computes a survival
as time until death, disease onset, or onset of an HPO term.

Any GPSEA analysis needs a genotype clasifier and a phenotype classifier/scorer/endpoint.
The following sections show how to prepare the classifiers, scorers, and endpoints
for your analysis. We also show a gallery with examples.


.. toctree::
  :maxdepth: 1
  :caption: Contents:

  genotype/index
  phenotype/index
  gallery
