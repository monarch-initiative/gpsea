.. _mode_of_inheritance_predicate:

==============================
Mode of Inheritance Predicates
==============================

There are five basic modes of inheritance for single-gene diseases: autosomal dominant, autosomal recessive, X-linked dominant, X-linked recessive, and mitochondrial (See
`Understanding Genetics, Appendix B <https://www.ncbi.nlm.nih.gov/books/NBK132145/>`_).


The :class:`~gpsea.analysis.predicate.genotype.ModeOfInheritancePredicate`
assigns the individual into a group based on the number of alleles
that match a condition specified by a :class:`~gpsea.analysis.predicate.genotype.VariantPredicate`.
The :class:`~gpsea.analysis.predicate.genotype.ModeOfInheritancePredicate` supports
the following Mendelian modes of inheritance (MoI):


+-----------------------+-----------------------------------+------------------+------------------------+
|  Mode of inheritance  | Sex                               |   Allele count   |  Genotype category     |
+=======================+===================================+==================+========================+
|  Autosomal dominant   | `*`                               |   0              |  `HOM_REF`             |
+                       +-----------------------------------+------------------+------------------------+
|                       | `*`                               |   1              |  `HET`                 |
+                       +-----------------------------------+------------------+------------------------+
|                       | `*`                               |   :math:`\ge 2`  |  ``None``              |
+-----------------------+-----------------------------------+------------------+------------------------+
|  Autosomal recessive  | `*`                               |   0              |  `HOM_REF`             |
+                       +-----------------------------------+------------------+------------------------+
|                       | `*`                               |   1              |  `HET`                 |
+                       +-----------------------------------+------------------+------------------------+
|                       | `*`                               |   2              |  `BIALLELIC_ALT`       |
+                       +-----------------------------------+------------------+------------------------+
|                       | `*`                               |   :math:`\ge 3`  |  ``None``              |
+-----------------------+-----------------------------------+------------------+------------------------+
|  X-linked dominant    | `*`                               |   0              |  `HOM_REF`             |
+                       +-----------------------------------+------------------+------------------------+
|                       | `*`                               |   1              |  `HET`                 |
+                       +-----------------------------------+------------------+------------------------+
|                       | `*`                               |   :math:`\ge 2`  |  ``None``              |
+-----------------------+-----------------------------------+------------------+------------------------+
|  X-linked recessive   | `*`                               |   0              |  `HOM_REF`             |
+                       +-----------------------------------+------------------+------------------------+
|                       | :class:`~gpsea.model.Sex.FEMALE`  |   1              |  `HET`                 |
+                       +                                   +------------------+------------------------+
|                       |                                   |   2              |  `BIALLELIC_ALT`       |
+                       +                                   +------------------+------------------------+
|                       |                                   |   :math:`\ge 3`  |  ``None``              |
+                       +-----------------------------------+------------------+------------------------+
|                       | :class:`~gpsea.model.Sex.MALE`    |   1              |  `HEMI`                |
+                       +                                   +------------------+------------------------+
|                       |                                   |   :math:`\ge 2`  |  ``None``              |
+-----------------------+-----------------------------------+------------------+------------------------+

.. note::

    `BIALLELIC_ALT` includes both homozygous and compound heterozygous genotypes.

Clinical judgment should be used to choose the MoI for the cohort analysis.
Then a predicate for the desired MoI can be created by one of 
:class:`~gpsea.analysis.predicate.genotype.ModeOfInheritancePredicate` static constructors:

* :func:`~gpsea.analysis.predicate.genotype.ModeOfInheritancePredicate.autosomal_dominant`
* :func:`~gpsea.analysis.predicate.genotype.ModeOfInheritancePredicate.autosomal_recessive`
* :func:`~gpsea.analysis.predicate.genotype.ModeOfInheritancePredicate.x_dominant`
* :func:`~gpsea.analysis.predicate.genotype.ModeOfInheritancePredicate.x_recessive`

All constructors take an instance
of :class:`~gpsea.analysis.predicate.genotype.VariantPredicate` as an argument.


Example
-------

Here we show seting up a predicate for grouping individuals based on
having a variant that leads to a frameshift or to a stop gain to a fictional transcript ``NM_1234.5``
to test differences between the genotypes of a disease with an autosomal recessive MoI.

First, we set up a :class:`~gpsea.analysis.predicate.genotype.VariantPredicate`
for testing if a variant meets the condition:

>>> from gpsea.model import VariantEffect
>>> from gpsea.analysis.predicate.genotype import VariantPredicates
>>> tx_id = 'NM_1234.5'
>>> is_frameshift_or_stop_gain = VariantPredicates.variant_effect(VariantEffect.FRAMESHIFT_VARIANT, tx_id) \
...     | VariantPredicates.variant_effect(VariantEffect.STOP_GAINED, tx_id)
>>> is_frameshift_or_stop_gain.get_question()
'(FRAMESHIFT_VARIANT on NM_1234.5 OR STOP_GAINED on NM_1234.5)'

Next, we use :class:`~gpsea.analysis.predicate.genotype.ModeOfInheritancePredicate.autosomal_recessive`
for assigning a patient into a genotype group:

>>> from gpsea.analysis.predicate.genotype import ModeOfInheritancePredicate
>>> gt_predicate = ModeOfInheritancePredicate.autosomal_recessive(is_frameshift_or_stop_gain)
>>> gt_predicate.display_question()
'What is the genotype group?: HOM_REF, HET, BIALLELIC_ALT'

The `gt_predicate` can be used in downstream analysis, such as in :class:
