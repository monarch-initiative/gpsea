.. _mode-of-inheritance-predicate:

==============================
Mode of Inheritance Predicates
==============================

There are five basic modes of inheritance (MoI) for single-gene diseases: autosomal dominant, 
autosomal recessive, X-linked dominant, X-linked recessive, and mitochondrial 
(See `Understanding Genetics, Appendix B <https://www.ncbi.nlm.nih.gov/books/NBK132145/>`_).

GPSEA enables testing aligned with an expected MoI, and to categorize individuals
based on the number of variant alleles. Clinical judgment should be used to choose the MoI
for the cohort analysis.

.. 

    By default, MoI predicates use *all* variants recorded in the individual.
    However, a :class:`~gpsea.analysis.predicate.genotype.VariantPredicate`
    can be provided to select the variants of interest. For instance, missense variants
    or the variants that overlap with a specific protein domain.


.. _autosomal-dominant-moi:

******************
Autosomal dominant
******************

In case of diseases with the autosomal dominant MoI, we typically assemble a cohort of individuals
harboring one pathogenic variant allele and we analyze the phenotypic differences between variant qualities
(e.g. missense vs. loss-of-function). As pointed out in the :ref:`monoallelic-predicate` section,
we use :class:`~gpsea.analysis.predicate.genotype.monoallelic_predicate`.
However, we can also use the predicate to compare the phenotype of the individuals with an allele
with those without the allele, and assign the individuals into the respective groups.


Example
=======

We may want to compare the individuals with a missense variant on a fictional transcript ``NM_1234.5``
which are predicted to affect one of the ten aminoacid residues :math:`[31, 40]` of a protein domain of interest:

We start by creating the variant predicate, which needs a variant to be both a missense variant
and overlap with the selected aminoacid residues:

>>> from gpsea.model import VariantEffect
>>> from gpsea.analysis.predicate.genotype import VariantPredicates
>>> is_missense = VariantPredicates.variant_effect(effect=VariantEffect.MISSENSE_VARIANT, tx_id="NM_1234.5")
>>> overlaps_domain = VariantPredicates.region(region=(31, 40), tx_id="NM_1234.5")
>>> var_predicate = is_missense & overlaps_domain
>>> var_predicate.get_question()
'(MISSENSE_VARIANT on NM_1234.5 AND variant affects aminoacid(s) between 31 and 40 on protein encoded by transcript NM_1234.5)'


and we use the variant predicate to create a genotype predicate for comparing the individuals with the
missense variant in the domain (`Missense in domain`) with those with no such variant (`Other`):

>>> from gpsea.analysis.predicate.genotype import monoallelic_predicate
>>> gt_predicate = monoallelic_predicate(
...     var_predicate,
...     names=("Missense in domain", "Other"),
... )
>>> gt_predicate.display_question()
'Allele group: Missense in domain, Other'



.. _autosomal-recessive-predicate:

*****************************
Autosomal recessive predicate
*****************************

The autosomal recessive predicate uses the allele count to assign
an individual into one of the genotype categories:

.. table:: Autosomal recessive predicate categories

    +------------------+-------------------+----------------+
    |   Allele count   |  Category         | Category index |
    +==================+===================+================+
    |   0              |  `No allele`      | 0              |
    +------------------+-------------------+----------------+
    |   1              |  `Monoallelic`    | 1              |
    +------------------+-------------------+----------------+
    |   2              |  `Biallelic`      | 2              |
    +------------------+-------------------+----------------+
    |   :math:`\ge 3`  |  ``None``         |                |
    +------------------+-------------------+----------------+

.. note::

    `Biallelic` includes both homozygous and compound heterozygous genotypes.


Partitions
==========

Sometimes we are interested in lumping several genotype categories into a group or and then comparing the groups.
For instance, to compare phenotype of the individuals with *at least one* frameshift allele
with those with *no* frameshift allele. Alternatively, we may only want to analyze a subset of the genotype categories,
such as `Monoallelic` vs. `Biallelic`.

The `partitions` option of the :func:`~gpsea.analysis.predicate.genotype.autosomal_recessive` function
lets us do this.
The option needs a set of sets of category indices (see table above).
The set is a partition of a set with the following properties:

* no subset is empty
* the intersection of any subsets is empty

These rules are very much alike the properties of the `partitions of a set <https://en.wikipedia.org/wiki/Partition_of_a_set>`_,
with the exception that adherence to *the union of the subsets include all group indices* rule is *NOT* required.


Examples
========


Use all variants
----------------

We create the predicate
with the :func:`~gpsea.analysis.predicate.genotype.autosomal_recessive` function:

>>> from gpsea.analysis.predicate.genotype import autosomal_recessive
>>> gt_predicate = autosomal_recessive()
>>> gt_predicate.display_question()
'What is the genotype group: No allele, Monoallelic, Biallelic'


Use a variant subset
--------------------

Same as in the autosomal dominant version,
we can use a :class:`~gpsea.analysis.predicate.genotype.VariantPredicate`
to only count the alleles of the variants of interest, such as the missense variants
of a fictional transcript ``NM_1234.5``:

>>> from gpsea.model import VariantEffect
>>> from gpsea.analysis.predicate.genotype import VariantPredicates
>>> is_missense = VariantPredicates.variant_effect(effect=VariantEffect.MISSENSE_VARIANT, tx_id="NM_1234.5")
>>> is_missense.get_question()
'MISSENSE_VARIANT on NM_1234.5'

and then use it to create the autosomal recessive predicate:

>>> gt_predicate = autosomal_recessive(is_missense)
>>> gt_predicate.display_question()
'What is the genotype group: No allele, Monoallelic, Biallelic'

This predicate will assign the individuals into one of the listed genotype categories
based on the allele counts of the missense variants.


Compare `Monoallelic` vs. `Biallelic`
-------------------------------------

We can provide ``partitions`` to only compare the heterozygotes with those carrying
biallelic alt mutations (homozygous alternate or compound heterozygous):

We consult the *Autosomal recessive predicate categories* table for the category indices
and we create the genotype group partitions:

>>> # `1` for `Monoallelic` and `2` for `Biallelic`
>>> partitions = ({1,}, {2,})

which we use to create the autosomal recessive predicate:

>>> gt_predicate = autosomal_recessive(
...     partitions=partitions,    
... )
>>> gt_predicate.display_question()
'What is the genotype group: Monoallelic, Biallelic'
