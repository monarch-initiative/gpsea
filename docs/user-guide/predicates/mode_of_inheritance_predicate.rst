.. _mode-of-inheritance-predicate:

==============================
Mode of Inheritance Predicates
==============================

There are five basic modes of inheritance for single-gene diseases: autosomal dominant, 
autosomal recessive, X-linked dominant, X-linked recessive, and mitochondrial 
(See `Understanding Genetics, Appendix B <https://www.ncbi.nlm.nih.gov/books/NBK132145/>`_).


The :class:`~gpsea.analysis.predicate.genotype.ModeOfInheritancePredicate`
assigns the individual into a group based on the number of alleles
that match a condition specified by a :class:`~gpsea.analysis.predicate.genotype.VariantPredicate`.
The :class:`~gpsea.analysis.predicate.genotype.ModeOfInheritancePredicate` supports
the following Mendelian modes of inheritance (MoI):


+-----------------------+------------------+------------------------+
|  Mode of inheritance  |   Allele count   |  Genotype category     |
+=======================+==================+========================+
|  Autosomal dominant   |   0              |  `HOM_REF`             |
+                       +------------------+------------------------+
|                       |   1              |  `HET`                 |
+                       +------------------+------------------------+
|                       |   :math:`\ge 2`  |  ``None``              |
+-----------------------+------------------+------------------------+
|  Autosomal recessive  |   0              |  `HOM_REF`             |
+                       +------------------+------------------------+
|                       |   1              |  `HET`                 |
+                       +------------------+------------------------+
|                       |   2              |  `BIALLELIC_ALT`       |
+                       +------------------+------------------------+
|                       |   :math:`\ge 3`  |  ``None``              |
+-----------------------+------------------+------------------------+


.. note::

    `BIALLELIC_ALT` includes both homozygous and compound heterozygous genotypes.

Clinical judgment should be used to choose the MoI for the cohort analysis.
Then a predicate for the desired MoI can be created by one of 
:class:`~gpsea.analysis.predicate.genotype.ModeOfInheritancePredicate` static constructors:

* :func:`~gpsea.analysis.predicate.genotype.ModeOfInheritancePredicate.autosomal_dominant`
* :func:`~gpsea.analysis.predicate.genotype.ModeOfInheritancePredicate.autosomal_recessive`

By default, the MoI predicates will use *all* variants recorded in the individual.
However, a :class:`~gpsea.analysis.predicate.genotype.VariantPredicate`
can be provided to select a variant subset, if necessary.


Assign individuals into genotype groups
---------------------------------------

Here we show seting up a predicate for grouping individuals for differences
between genotypes of a disease with an autosomal recessive MoI.

We use :class:`~gpsea.analysis.predicate.genotype.ModeOfInheritancePredicate.autosomal_recessive`
to create the predicate:

>>> from gpsea.analysis.predicate.genotype import ModeOfInheritancePredicate
>>> gt_predicate = ModeOfInheritancePredicate.autosomal_recessive()
>>> gt_predicate.display_question()
'What is the genotype group: HOM_REF, HET, BIALLELIC_ALT'

The predicate will use *all* recorded variants to determine if the individual belongs into
homozygous reference (`HOM_REF`), heterozygous (`HET`), or biallelic alternative (`BIALLELIC_ALT`)
category.


Use a subset of variants for choosing the genotype group
--------------------------------------------------------

To select specific variants, a :class:`~gpsea.analysis.predicate.genotype.VariantPredicate`
can be registered with the MoI predicate.
For instance, the following can be done to only consider the variants that lead
to a missense change on a fictional transcript ``NM_1234.5``
when assigning the genotype group. We set up the variant predicate:

>>> from gpsea.model import VariantEffect
>>> from gpsea.analysis.predicate.genotype import VariantPredicates
>>> tx_id = 'NM_1234.5'
>>> is_missense = VariantPredicates.variant_effect(VariantEffect.MISSENSE_VARIANT, tx_id)
>>> is_missense.get_question()
'MISSENSE_VARIANT on NM_1234.5'

and we use it to create the MoI predicate:

>>> gt_predicate = ModeOfInheritancePredicate.autosomal_recessive(is_missense)
>>> gt_predicate.display_question()
'What is the genotype group: HOM_REF, HET, BIALLELIC_ALT'
