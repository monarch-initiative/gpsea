.. _variant-category:

=========================
Group by variant category
=========================

Sometimes we want to compare the individuals who have the same allele count (`AC`) of variant categories :math:`A` and :math:`B`.
For example, in the context of an autosomal dominant disease,
we may want to compare the individuals with :math:`AC_{A}=1` (where :math:`A` is e.g. a predicted loss-of-function mutation)
with those harboring :math:`AC_{B}=1` (where :math:`B` is e.g. a missense mutation).
Similarly, in an autosomal recessive disease, we may be interested in comparing the individuals
with :math:`AC_{A} \ge 1` with those with :math:`AC_{A} = 0`.
In both analyses, we compare two variant categories :math:`A` and :math:`B`
which are described by a :class:`~gpsea.analysis.predicate.VariantPredicate`
(see :ref:`variant-predicates` section),
while ensuring the allele count sum of both variant categories is :math:`k`.

:math:`k = \sum_{i \in \{A, B\}} AC_{i}`

GPSEA provides two genotype classifiers:


.. table::

    +-----------+------------------------+------------------------------------------------------------------+
    | :math:`k` | Name                   | Function                                                         |
    +===========+========================+==================================================================+
    | 1         | Monoallelic classifier | :func:`~gpsea.analysis.clf.monoallelic_classifier`               |
    +-----------+------------------------+------------------------------------------------------------------+
    | 2         | Biallelic classifier   | :func:`~gpsea.analysis.clf.biallelic_classifier`                 |
    +-----------+------------------------+------------------------------------------------------------------+


.. _monoallelic-classifier:

**********************
Monoallelic classifier
**********************

Monoallelic classifier compares the individuals who have *one* allele of a variants of interest.
The classifier uses a pair of variant predicates, `A` and `B`,
to compute the allele counts :math:`AC_{A}` and :math:`AC_{B}`,
in order to assign an individual into a genotype class.

.. table:: Monoallelic genotype classes

    +-----------------+-------------------+-------------------+
    | Genotype class  | :math:`AC_{A}`    | :math:`AC_{B}`    |
    +=================+===================+===================+
    | A               | 1                 | 0                 |
    +-----------------+-------------------+-------------------+
    | B               | 0                 | 1                 |
    +-----------------+-------------------+-------------------+
    | ``None``        | other             | other             |
    +-----------------+-------------------+-------------------+

The individuals with :math:`\sum_{i \in \{A, B\}} AC_{i} \neq 1`
are omitted from the analysis.

.. figure:: img/monoallelic-classifier.png
   :alt: Monoallelic classifier
   :align: center
   :width: 600px


Example
=======

Let's create a classifier to categorize the individuals
to those having one missense allele or to those having
one frameshift allele with respect to fictional transcript ``NM_1234.5``.

>>> tx_id = "NM_1234.5"

We start by creating the variant predicates for missense (`A`)
and frameshift (`B`) variants:

>>> from gpsea.model import VariantEffect
>>> from gpsea.analysis.predicate import variant_effect
>>> is_missense = variant_effect(effect=VariantEffect.MISSENSE_VARIANT, tx_id=tx_id)
>>> is_frameshift = variant_effect(effect=VariantEffect.FRAMESHIFT_VARIANT, tx_id=tx_id)

Monoallelic classifier lets us customize the category names.
Let's use `Missense` and `Frameshift` instead of the defaults:

>>> a_label = "Missense"
>>> b_label = "Frameshift"

Now we have all we need to create the predicate:

>>> from gpsea.analysis.clf import monoallelic_classifier
>>> gt_clf = monoallelic_classifier(
...     a_predicate=is_missense,
...     b_predicate=is_frameshift,
...     a_label=a_label, b_label=b_label,
... )
>>> gt_clf.class_labels
('Missense', 'Frameshift')



.. _biallelic-classifier:

********************
Biallelic classifier
********************

Biallelic classifier compares the individuals with *two* alleles of the variants of interest.
The functionality is very similar to that of monoallelic classifier, with two differences:
(1) genotype classes and (2) partitions.


Genotype classes
================

Biallelic locus can be present in one of three genotypes, allowing an individual
to be assigned into one of the three genotype classes:

.. _biallelic-gt-classes:

.. table:: Biallelic genotype classes

    +-------+----------------+-------------------+-------------------+
    | Index | Genotype class | :math:`AC_{A}`    | :math:`AC_{B}`    |
    +=======+================+===================+===================+
    | 0     | A/A            | 2                 | 0                 |
    +-------+----------------+-------------------+-------------------+
    | 1     | A/B            | 1                 | 1                 |
    +-------+----------------+-------------------+-------------------+
    | 2     | B/B            | 0                 | 2                 |
    +-------+----------------+-------------------+-------------------+
    |       | ``None``       | other             | other             |
    +-------+----------------+-------------------+-------------------+

Note that :math:`\sum_{i \in \{A, B\}} AC_{i} = 2` and the individuals
with a different allele count sum are omitted from the analysis.


.. figure:: img/biallelic-classifier.png
   :alt: Biallelic classifier
   :align: center
   :width: 600px


Example
-------

Let `A` and `B` correspond to *MISSENSE* and *FRAMESHIFT* variants,
and let's reuse the variant predicates ``is_missense`` and ``is_frameshift`` from the previous section,
to compare missense and frameshift variants in the context of an autosomal recessive disease.

>>> from gpsea.analysis.clf import biallelic_classifier
>>> gt_clf = biallelic_classifier(
...     a_predicate=is_missense,
...     b_predicate=is_frameshift,
...     a_label="Missense", b_label="Frameshift",
... )
>>> gt_clf.class_labels
('Missense/Missense', 'Missense/Frameshift', 'Frameshift/Frameshift')


The classifier assigns an individual into one of three genotype classes:

* `Missense/Missense` - two missense alleles
* `Missense/Frameshift` - one missense and one frameshift allele
* `Frameshift/Frameshift` - two frameshift alleles

    
Partitions
==========

Sometimes we are interested in lumping several genotype classes into a partition
and then comparing the partitions.
For instance, in the context of an autosomal recessive disease,
we may want to compare individuals with two "mild" mutations with the individuals
with at least one "severe" mutation.
This comparison can be implemented using the `partitions` option.

A partition is a set of one or more genotype class indices
(see :ref:`biallelic-gt-classes` table).
Then, two (or more) partitions are provided to biallelic classifier
via the `partitions` option.

For example, we can compare the individuals with two missense alleles with those harboring
one frameshift and one missense alleles, or two frameshift alleles.

Let `A` and `B` correspond to *MISSENSE* and *FRAMESHIFT* variant.
According to :ref:`biallelic-gt-classes` table,
the `A/A` genotype class corresponds to index `0`,
and the `A/B` and `B/B` genotype class correspond to indices `1` and `2`, respectively.
We form the partitions accordingly:

>>> partitions = (0, {1, 2})

With the ``partitions``, the biallelic classifier splits the individuals
into two classes:

* two missense alleles
* one missense alelele and one frameshift allele or two frameshift alleles

.. figure:: img/biallelic-classifier-w-partitions.png
   :alt: Biallelic classifier with partitions
   :align: center
   :width: 600px


Using the ``partitions`` in code is a no-brainer:

>>> gt_clf = biallelic_classifier(
...     a_predicate=is_missense,
...     b_predicate=is_frameshift,
...     a_label="Missense", b_label="Frameshift",
...     partitions=partitions,
... )
>>> gt_clf.class_labels
('Missense/Missense', 'Missense/Frameshift OR Frameshift/Frameshift')
