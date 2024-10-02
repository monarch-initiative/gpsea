.. _allele-count-predicates:

=======================
Allele count predicates
=======================

Sometimes we want to analyze the individuals who have the same allele count of the variants of interest.
Most of the time, the allele count is `1` to analyze individuals with an autosomal dominant disease,
or `2` for the individuals with a disease segregating with the autosomal recessive mode of inheritance.

GPSEA provides :func:`~gpsea.analysis.predicate.genotype.monoallelic_predicate`
and :func:`~gpsea.analysis.predicate.genotype.biallelic_predicate` to investigate
the individuals with `1` or `2` alleles of the variants of interest.

Both predicates take two variant predicates (:class:`~gpsea.analysis.predicate.genotype.VariantPredicate`)
that are used to determine the allele count of the target variants. 


.. _monoallelic-predicate:

*********************
Monoallelic predicate
*********************

Monoallelic predicate compares individuals who have *one* allele of a variant of interest.
The predicate requires two variant predicates `A` and `B` for choosing
the variants of interest (e.g. variants overlapping with a specific protein domain).
The allele counts :math:`A_{AC}` and :math:`B_{AC}` for the variants of interest
are computed and the individual is assigned into a category
based on the following table:

.. table:: Monoallelic predicate genotype categories

    +-------------------+-------------------+-----------+
    | :math:`A_{AC}`    | :math:`B_{AC}`    | Category  |
    +===================+===================+===========+
    | 1                 | 0                 | A         |
    +-------------------+-------------------+-----------+
    | 0                 | 1                 | B         |
    +-------------------+-------------------+-----------+
    | other             | other             | ``None``  |
    +-------------------+-------------------+-----------+

The individuals with different allele counts
(e.g. :math:`A_{AC} = 0` and :math:`B_{AC} = 2`)
are assigned into the ``None`` group and, thus, omitted from the analysis.


Example
=======

Let's create a predicate to categorize the individuals
to those having one missense allele or to those having
one frameshift allele with respect to fictional transcript ``NM_1234.5``.

We start by creating the variant predicates for missense (`A`)
and frameshift (`B`) variants:

>>> from gpsea.model import VariantEffect
>>> from gpsea.analysis.predicate.genotype import VariantPredicates
>>> tx_id = "NM_1234.5"
>>> is_missense = VariantPredicates.variant_effect(effect=VariantEffect.MISSENSE_VARIANT, tx_id=tx_id)
>>> is_frameshift = VariantPredicates.variant_effect(effect=VariantEffect.FRAMESHIFT_VARIANT, tx_id=tx_id)

Monoallelic predicate lets us customize the category names.
Let's use `Missense` and `Frameshift` instead of the defaults `A` and `B`:

>>> names = ("Missense", "Frameshift")

Now we have all we need to create the monoallelic predicate:

>>> from gpsea.analysis.predicate.genotype import monoallelic_predicate
>>> gt_predicate = monoallelic_predicate(
...     a_predicate=is_missense,
...     b_predicate=is_frameshift,
...     names=names,
... )
>>> gt_predicate.display_question()
'Allele group: Missense, Frameshift'

`gt_predicate` can be used in any GPSEA analysis, to compare the individuals who harbor
one allele of a missense variant with the individuals with one allele that leads to a frameshift.


.. _biallelic-predicate:

*******************
Biallelic predicate
*******************

Biallelic predicate compares the individuals with *two* alleles of the variants of interest.
The functionality is very similar to that of monoallelic predicate, with two differences.


Categories
----------

A biallelic locus can be present in one of three genotypes.
For example, using the predicates `A` and `B` for the missense vs. frameshift comparison as in the monoallelic case,
we will assign the individuals into one of the following genotype categories:

.. table:: Biallelic predicate genotype categories

    +-------------------+-------------------+-----------+----------------+
    | :math:`A_{AC}`    | :math:`B_{AC}`    | Category  | Category index |
    +===================+===================+===========+================+
    | 2                 | 0                 | A/A       | 0              |
    +-------------------+-------------------+-----------+----------------+
    | 1                 | 1                 | A/B       | 1              |
    +-------------------+-------------------+-----------+----------------+
    | 0                 | 2                 | B/B       | 2              |
    +-------------------+-------------------+-----------+----------------+
    | other             | other             | ``None``  |                |
    +-------------------+-------------------+-----------+----------------+

The individuals with a different allele count combination
(e.g. :math:`A_{AC} = 1` and :math:`B_{AC} = 2`)
are assigned into the ``None`` group and will be, thus, omitted from the analysis.

    
Partitions
----------

Sometimes we are interested in lumping several genotype categories into a group and then comparing the groups.
For instance, we may want to to compare phenotype of the individuals with *at least one* frameshift allele
with those with *no* frameshift allele.

The `partitions` option dictates form of the genotype categories.
The option needs a set of sets of category indices (see table above).
The set is a partition of a set with the following properties:

* no subset is empty
* the union of the subsets include all group indices
* the intersection of any subsets is empty

as in `partitions of a set <https://en.wikipedia.org/wiki/Partition_of_a_set>`_.


Examples
^^^^^^^^

Let `A` and `B` correspond to the variant predicates that select *MISSENSE* and *FRAMESHIFT* variants,
and let's reuse the variant predicates ``is_missense`` and ``is_frameshift`` from the previous section.


Compare missense vs. frameshift
-------------------------------

We compare missense and frameshift variants in the context of an autosomal recessive disease,
and we will use the same allele group names as before - `Missense` and `Frameshift`.


>>> from gpsea.analysis.predicate.genotype import biallelic_predicate
>>> gt_predicate = biallelic_predicate(
...     a_predicate=is_missense,
...     b_predicate=is_frameshift,
...     names=names,
... )
>>> gt_predicate.display_question()
'Allele group: Missense/Missense, Missense/Frameshift, Frameshift/Frameshift'

The predicate will assign the individuals into one of three genotype groups:

* `Missense/Missense` - individual with two missense alleles
* `Missense/Frameshift` - individual with one missense allele and one frameshift allele
* `Frameshift/Frameshift` - individual with two frameshift alleles

.. note::

    This corresponds to the partitions::

        partitions = ({0,}, {1,}, {2,})


Compare missense vs. frameshift
-------------------------------

Here we compare individuals with at least one frameshift allele with those with only missense variants.
We can achieve this by providing the `partitions` option:

>>> partitions = ({0,}, {1, 2})

The partitions includes subsets of the category indices to lump together.
With reference to the *Biallelic predicate genotype categories* table,
the first subset represents the individuals with `Missense/Missense` genotype
and the second subset represents those with either `Missense/Frameshift` or `Frameshift/Frameshift`.

We provide `partitions` as an extra argument
to the :func:`~gpsea.analysis.predicate.genotype.biallelic_predicate` function:

>>> gt_predicate = biallelic_predicate(
...     a_predicate=is_missense,
...     b_predicate=is_frameshift,
...     names=names,
...     partitions=partitions,
... )
>>> gt_predicate.display_question()
'Allele group: Missense/Missense, Missense/Frameshift OR Frameshift/Frameshift'

