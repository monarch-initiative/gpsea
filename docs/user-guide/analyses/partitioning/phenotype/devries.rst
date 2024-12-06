.. _devries-scorer:

##############
De Vries Score
##############


The De Vries score was developed as a relatively simple phenotypic severity
score for individuals with intellectual disability in which points are given
for (severity of) intellectual disability,
growth abnormalities (prenatal and postnatal), facial dysmorphisms,
nonfacial dysmorphisms, and other congenital anomalies
(`Dingemans et al. (2022) <https://pubmed.ncbi.nlm.nih.gov/36182950/>`_).
Statistical significance of a difference in the De Vries score between groups can be
determined using the Mann-Whitney-U test.

We refer to `Feenstra et al. (2011) <https://pubmed.ncbi.nlm.nih.gov/21712853/>`_ for
the original description of the adjusted De Vries score. Here we offer an adapted version of the
score that leverages the structure of the Human Phenotype Ontology to assess the phenotype.


*************************
The sections of the score
*************************

The De Vries score has several sections, each of which is scored on a point system. The
final score is obtained as the sum of the points of each of the sections.

Developmental delay
~~~~~~~~~~~~~~~~~~~

The original score assigns one point for mild or moderate developmental delay
and two points for severe developmental delay.

+----------------------------------------------------------------------------------------------------------+-----------+
| HPO term                                                                                                 | Score     |
+==========================================================================================================+===========+
| `Mild global developmental delay (HP:0011342) <https://hpo.jax.org/browse/term/HP:0011342>`_             | 1         |
+----------------------------------------------------------------------------------------------------------+-----------+
| `Moderate global developmental delay (HP:0011343)  <https://hpo.jax.org/browse/term/HP:0011343>`_        | 1         |
+----------------------------------------------------------------------------------------------------------+-----------+
| `Severe global developmental delay (HP:0011344)  <https://hpo.jax.org/browse/term/HP:0011344>`_          | 2         |
+----------------------------------------------------------------------------------------------------------+-----------+
| `Profound global developmental delay (HP:0011344)  <https://hpo.jax.org/browse/term/HP:0012736>`_        | 2         |
+----------------------------------------------------------------------------------------------------------+-----------+
| `Global developmental delay (HP:0001263)  <https://hpo.jax.org/browse/term/HP:0012736>`_                 | 1         |
+----------------------------------------------------------------------------------------------------------+-----------+

Note that one point is assigned for the term Global developmental delay (HP:0001263), which is the
parent of the other terms, because no information was provided about the degree of delay.

If none of the above terms is found, then the scorer assigns terms based on the Intellectual Disability terms.

+----------------------------------------------------------------------------------------------------------+-----------+
| HPO term                                                                                                 | Score     |
+==========================================================================================================+===========+
| `Intellectual disability, borderline (HP:0006889) <https://hpo.jax.org/browse/term/HP:0006889>`_         | 0.5       |
+----------------------------------------------------------------------------------------------------------+-----------+
| `Intellectual disability, mild (HP:0001256)  <https://hpo.jax.org/browse/term/HP:0001256>`_              | 1         |
+----------------------------------------------------------------------------------------------------------+-----------+
| `Intellectual disability, moderate (HP:0002342)  <https://hpo.jax.org/browse/term/HP:0002342>`_          | 1         |
+----------------------------------------------------------------------------------------------------------+-----------+
| `Intellectual disability, severe (HP:0010864)  <https://hpo.jax.org/browse/term/HP:0010864>`_            | 2         |
+----------------------------------------------------------------------------------------------------------+-----------+
| `Intellectual disability, profound (HP:0002187)  <https://hpo.jax.org/browse/term/HP:0002187>`_          | 2         |
+----------------------------------------------------------------------------------------------------------+-----------+
| `Intellectual disability (HP:0001249)  <https://hpo.jax.org/browse/term/HP:0001249>`_                    | 1         |
+----------------------------------------------------------------------------------------------------------+-----------+

If none of these terms is found, a score of zero is assigned for this section.


Prenatal-onset growth retardation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
In the original score, two points are assigned if Prenatal-onset growth retardation is present. In our implementation,
we assign two points if either of the following terms is present (the score is thus either zero or two).

+----------------------------------------------------------------------------------------------------------+-----------+
| HPO term                                                                                                 | Score     |
+==========================================================================================================+===========+
| `Small for gestational age (HP:0001518) <https://hpo.jax.org/browse/term/HP:0001518>`_                   | 2         |
+----------------------------------------------------------------------------------------------------------+-----------+
| `Intrauterine growth retardation (HP:0001511)  <https://hpo.jax.org/browse/term/HP:0001511>`_            | 2         |
+----------------------------------------------------------------------------------------------------------+-----------+




Postnatal growth abnormalities
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the original score, one point is assigned for each of Microcephaly, Short stature, Macrocephaly, and Tall stature,
with the maximum score for the section being limited to 2 points. We implement this as adding one point for either the
corresponding HPO terms or any of their descendents (for instance, `Disproportionate short stature (HP:0003498) <https://hpo.jax.org/browse/term/HP:0003498>`_ would
be counted for `Short stature (HP:0004322) <https://hpo.jax.org/browse/term/HP:0004322>`_).

+----------------------------------------------------------------------------------------------------------+-----------+
| HPO term                                                                                                 | Score     |
+==========================================================================================================+===========+
| `Microcephaly (HP:0000252) <https://hpo.jax.org/browse/term/HP:0000252>`_                                | 1         |
+----------------------------------------------------------------------------------------------------------+-----------+
| `Short stature (HP:0004322) <https://hpo.jax.org/browse/term/HP:0004322>`_                               | 1         |
+----------------------------------------------------------------------------------------------------------+-----------+
| `Macrocephaly (HP:0000256)  <https://hpo.jax.org/browse/term/HP:0000256>`_                               | 1         |
+----------------------------------------------------------------------------------------------------------+-----------+
| `Tall stature (HP:0000098)  <https://hpo.jax.org/browse/term/HP:0010864>`_                               | 1         |
+----------------------------------------------------------------------------------------------------------+-----------+

If none of these terms is found, a score of zero is assigned for this section. Logically, the maximum score obtainable
is 2 because the same individual cannot have both tall and short stature or both micro- and macrocephaly.


Facial dysmorphic features
~~~~~~~~~~~~~~~~~~~~~~~~~~

This section assigns two points if two or more facial dysmorphisms are identified. In contrast to the list of anomalies described
in the original 2011 publication of the DeVries score, we leverage the structure of the HPO to include many more HPO terms that 
denote various kinds of facial dysmorphism (e.g., `Abnormality of globe location <https://hpo.jax.org/browse/term/HP:0100886>`_ instead of just
`Hypertelorism (HP:0000316) <https://hpo.jax.org/browse/term/HP:0000316>`_).

Our implementation of this feature counts the total number of terms or descendents of the following HPO terms. Up to one point is given
for each of the categories.

+----------------------------------------------------------------------------------------------------------+-----------+
| HPO term                                                                                                 | Score     |
+==========================================================================================================+===========+
| `Abnormality of globe location (HP:0000316) <https://hpo.jax.org/browse/term/HP:0100886>`_               | 0 or 1    |
+----------------------------------------------------------------------------------------------------------+-----------+
| `Abnormal lip morphology (HP:0000159) <https://hpo.jax.org/browse/term/HP:0000159>`_                     |  0 or 1   |
+----------------------------------------------------------------------------------------------------------+-----------+
| `Abnormal facial shape (HP:0001999) <https://hpo.jax.org/browse/term/HP:0001999>`_                       |  0 or 1   |
+----------------------------------------------------------------------------------------------------------+-----------+
| `Abnormal midface morphology (HP:0000309) <https://hpo.jax.org/browse/term/HP:0000309>`_                 |  0 or 1   |
+----------------------------------------------------------------------------------------------------------+-----------+
| `Abnormal forehead morphology (HP:0000290) <https://hpo.jax.org/browse/term/HP:0000290>`_                | 0 or 1    |
+----------------------------------------------------------------------------------------------------------+-----------+
| `Abnormal chin morphology (HP:0000306) <https://hpo.jax.org/browse/term/HP:0000306>`_                    |  0 or 1   |
+----------------------------------------------------------------------------------------------------------+-----------+
| `Abnormal external nose morphology (HP:0010938) <https://hpo.jax.org/browse/term/HP:0010938>`_           |  0 or 1   |
+----------------------------------------------------------------------------------------------------------+-----------+
| `Abnormal pinna morphology (HP:0000377)  <https://hpo.jax.org/browse/term/HP:0000377>`_                  |  0 or 1   |
+----------------------------------------------------------------------------------------------------------+-----------+

If items from two or more categories are found, the score is 2, otherwise a score of zero is assigned.


Non-facial dysmorphism and congenital abnormalities
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
One point is assigned for either the corresponding HPO terms or any of their descendents up to a maximum of two points.
A maximum of one point is assigned for each of the following categories.

+----------------------------------------------------------------------------------------------------------+-----------+
| HPO term                                                                                                 | Score     |
+==========================================================================================================+===========+
| `Abnormal hand morphology (HP:0005922) <https://hpo.jax.org/browse/term/HP:0005922>`_                    | 0 or 1    |
+----------------------------------------------------------------------------------------------------------+-----------+
| `Abnormal heart morphology (HP:0001627) <https://hpo.jax.org/browse/term/HP:0001627>`_                   |  0 or 1   |
+----------------------------------------------------------------------------------------------------------+-----------+
| `Abnormal external genitalia morphology (HP:0000811) <https://hpo.jax.org/browse/term/HP:0000811>`_      |  0 or 1   |
+----------------------------------------------------------------------------------------------------------+-----------+

The score for this section can thus be 0, 1, or 2.


Final score
~~~~~~~~~~~

The final score is obtained by summing the scores from each of the sections. The final score ranges from 0 to 10, with
higher scores being considered a proxy for higher clinical severity.


*********************************
Using the De Vries scorer in code
*********************************

GPSEA implements the score in :class:`~gpsea.analysis.pscore.DeVriesPhenotypeScorer` that can be used
as a part of the :ref:`Phenotype score <phenotype-score-stats>` analysis, where it is used
as a :ref:`phenotype scorer <phenotype-score>`.

A De Vries scorer uses HPO hierarchy as a prerequisite.
We can load HPO using HPO toolkit:

>>> import hpotk
>>> store = hpotk.configure_ontology_store()
>>> hpo = store.load_minimal_hpo(release='v2024-07-01')

and use it to create :class:`~gpsea.analysis.pscore.DeVriesPhenotypeScorer`

>>> from gpsea.analysis.pscore import DeVriesPhenotypeScorer
>>> pheno_scorer = DeVriesPhenotypeScorer(hpo)
>>> pheno_scorer.description
'A phenotypic severity score for individuals with intellectual disability'

which we can use as a phenotype scorer.
