.. _devries:

==============
De Vries Score
==============


The De Vries score was developed as a relatively simple phenotypic severity
score for individuals with intellectual disability in which points are given
for (severity of) intellectual disability,
growth abnormalities (prenatal and postnatal), facial dysmorphisms,
nonfacial dysmorphisms, and other congenital anomalies
(`Dingemans et al. (2022) <https://pubmed.ncbi.nlm.nih.gov/36182950/>`_).
Statistical significance of a difference in the De Vries score between groups can be
determined using the Mann–Whitney-U test.

We refer to `Feenstra et al. (2011) <https://pubmed.ncbi.nlm.nih.gov/21712853/>`_ for
the original description of the adjusted De Vries score. Here we offer a version of the
score that leverages the structure of the Human Phenotype Ontology to assess the phenotype.

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

This section assigns two points if two or more anomalies are identified in the following
categories: hypertelorism, nasal anomalies and ear anomalies. Our implementation of this feature counts the total
number of terms or descendents of the following HPO terms.

+----------------------------------------------------------------------------------------------------------+-----------+
| HPO term                                                                                                 | Score     |
+==========================================================================================================+===========+
| `Hypertelorism (HP:0000316) <https://hpo.jax.org/browse/term/HP:0000316>`_                               | 1         |
+----------------------------------------------------------------------------------------------------------+-----------+
| `Abnormal external nose morphology (HP:0010938) <https://hpo.jax.org/browse/term/HP:0010938>`_           | 1 each    |
+----------------------------------------------------------------------------------------------------------+-----------+
| `Abnormal pinna morphology (HP:0000377)  <https://hpo.jax.org/browse/term/HP:0000377>`_                  | 1 each    |
+----------------------------------------------------------------------------------------------------------+-----------+

If two or more terms are found, the score is 2, otherwise a score of zero is assigned.


Non-facial dysmorphism and congenital abnormalities
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
One point is assigned for either the
corresponding HPO terms or any of their descendents up to a maximum of two points.

+----------------------------------------------------------------------------------------------------------+-----------+
| HPO term                                                                                                 | Score     |
+==========================================================================================================+===========+
| `Abnormal hand morphology (HP:0005922) <https://hpo.jax.org/browse/term/HP:0005922>`_                    | 1 each    |
+----------------------------------------------------------------------------------------------------------+-----------+
| `Abnormal heart morphology (HP:0001627) <https://hpo.jax.org/browse/term/HP:0001627>`_                   | 1 each    |
+----------------------------------------------------------------------------------------------------------+-----------+
| `Hypospadias (HP:0000047) <https://hpo.jax.org/browse/term/HP:0000047>`_                                 | 1         |
+----------------------------------------------------------------------------------------------------------+-----------+

Final score
~~~~~~~~~~~

The final score is obtained by summing the scores from each of the sections. The final score ranges from 0 to 10, with
higher scores being considered a proxy for higher clinical severity.


Using the De Vries Scorer in code
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TODO


