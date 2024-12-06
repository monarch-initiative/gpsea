.. _counting:


#########################
Counting phenotype scorer
#########################

The :class:`~gpsea.analysis.pscore.CountingPhenotypeScorer` assigns a phenotype score that is equivalent
to the count of observed phenotypes (HPO terms) that are either an exact match to the `query` terms or their descendants.
Typically, the `query` terms will comprise abnormalities in different organ systems.
For instance, we may want to count whether an individual has brain, liver, kidney, and skin abnormalities.
In the case, the query would include the corresponding terms (e.g., Abnormal brain morphology HP:0012443).
An individual can then have between 0 and 4 phenotype group abnormalities.  The scorer does not double count if the individual has multiple
observed abnormalities in one of the organ systems (i.e., multiple descendents of one of the query terms). Each individual can thus have a score 
of between 0 (no relevant abnormalities) to the number of categories (if the individual has an abnormality in each of the categories). 
The genotype groups are then compared with respect to the distribution of counts using the Mann Whitney U test.


*******
Example
*******

Here we use :class:`~gpsea.analysis.pscore.CountingPhenotypeScorer` for scoring
the individuals based on the number of structural defects
from the following 5 categories:

* Brain anomalies
* Eye anomalies
* Congenital heart defects
* Renal anomalies
* Sensorineural hearing loss

For example, an individual with a congenital heart defect would be assigned a score of `1`,
an individual with congenital heart defect and a renal anomaly would be assigned a score of `2`,
and so on. If an individual had two heart defects (e.g., atrial septal defect and ventricular septal defect), 
a score of 1 (not 2) would be assigned for the heart defect category.

The :class:`~gpsea.analysis.pscore.CountingPhenotypeScorer` automatizes this scoring method
by encoding the categories into HPO terms:

>>> structural_defects = (
...     'HP:0012443',  # Abnormal brain morphology (Brain anomalies)
...     'HP:0012372',  # Abnormal eye morphology (Eye anomalies)
...     'HP:0001627',  # Abnormal heart morphology (Congenital heart defects)
...     'HP:0012210',  # Abnormal renal morphology (Renal anomalies)
...     'HP:0000407',  # Sensorineural hearing impairment (Sensorineural hearing loss)
... )


and then tests the individuals for presence of at least one HPO term
that corresponds to the structural defect
(e.g. `Abnormal brain morphology <https://hpo.jax.org/browse/term/HP:0012443>`_, exact match)
or that is its descendant
(e.g. `Cerebellar atrophy <https://hpo.jax.org/browse/term/HP:0001272>`_).

We construct the scorer with
:func:`~gpsea.analysis.pscore.CountingPhenotypeScorer.from_query_curies` function:

>>> from gpsea.analysis.pscore import CountingPhenotypeScorer
>>> pheno_scorer = CountingPhenotypeScorer.from_query_curies(
...     hpo=hpo,
...     query=structural_defects,
... )
>>> pheno_scorer.description
'Assign a phenotype score that is equivalent to the count of present phenotypes that are either an exact match to the query terms or their descendants'
