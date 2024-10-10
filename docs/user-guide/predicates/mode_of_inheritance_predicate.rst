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

Note, the MoI analysis differs from the analysis based on :ref:`allele-count-predicates`;
the former compares phenotypes of individuals with different counts of the same variant category,
while the latter compares the individuals with the same allele count (e.g. 1, 2) of variants of different categories.


.. _autosomal-dominant-moi:

******************
Autosomal dominant
******************

Sometimes we want to compare phenotypes of individuals with one variant allele
with those harboring no such allele.
In GPSEA, this can be done using
:func:`~gpsea.analysis.predicate.genotype.monoallelic_predicate`


Example
=======

For instance, we may want to compare the individuals with one variant allele that overlaps with
the 6th exon of the fictional transcript `NM_1234.5` with those with no such allele.
We construct the corresponding genotype predicate in two steps.
First, we prepare a variant predicate for choosing the variants affecting exon 6:

>>> from gpsea.analysis.predicate.genotype import VariantPredicates
>>> in_exon_6 = VariantPredicates.exon(exon=6, tx_id="NM_1234.5")
>>> in_exon_6.get_question()
'variant affects exon 6 on NM_1234.5'

Second, by wrapping the `in_exon_6` predicate in :func:`~gpsea.analysis.predicate.genotype.monoallelic_predicate`,
we will partition the individuals into those with one variant allele (`Monoallelic`)
and those with no alleles (`No allele`), while omitting the individuals with 2 or more alleles:

>>> from gpsea.analysis.predicate.genotype import monoallelic_predicate
>>> gt_predicate = monoallelic_predicate(
...     a_predicate=in_exon_6,
...     names=("Monoallelic", "No allele"),
... )
>>> gt_predicate.display_question()
'Allele group: Monoallelic, No allele'

.. note::
    
    Under the hood, `monoallelic_predicate` creates the inverse of `in_exon_6`,
    to capture all variants outside exon 6.


.. _autosomal-recessive-moi:

*******************
Autosomal recessive
*******************

In case of autosomal recessive diseases, where biallelic mutations are needed for disease manifestation,
we can either compare the individuals with biallelic variants of different categories
(described in :ref:`biallelic-predicate`),
or the individuals with different allele counts of the same variant category (described here).

In theory, a biallelic locus can exist in one of the three states ``{"A/A", "A/B", "B/B"}``,
where `A` and `B` denote variant categories and, 
for instance, `A/A` represents a locus with two alleles of the `A` category,
`A/B` represents a locus with one allele each category.
In contrast to the analysis described in :ref:`biallelic-predicate` section,
the focus of autosomal recessive analysis is to compare the individuals
with different allele counts of the same variant category (e.g. `B`).
Despite the possibility of three states, most of the time, we are only interested in
comparing `A/B` with `B/B`.

GPSEA implements the allele counting with
:func:`~gpsea.analysis.predicate.genotype.allele_counting` predicate
which assigns individuals into one of the following genotype categories:


.. table:: Autosomal recessive categories

    +--------------------------+--------------------+
    |   Allele count           | Genotype category  |
    +==========================+====================+
    |   1                      |  `Monoallelic`     |
    +--------------------------+--------------------+
    |   2                      |  `Biallelic`       |
    +--------------------------+--------------------+
    |   :math:`\notin {1, 2}`  |  ``None``          |
    +--------------------------+--------------------+

.. note::

    `Biallelic` includes both homozygous and compound heterozygous genotypes.


Last, a :class:`~gpsea.model.Patient` can include variants unrelated to the "current" analysis.
These may include variants in a different gene, or variants that are otherwise unlikely to be relevant
to the analysis, and we may not want to include these in the allele counting. One way to solve this issue
is to ensure the phenopackets/individuals only include the relevant variants to start with.
However, as a convenience, :func:`~gpsea.analysis.predicate.genotype.allele_counting` takes an optional
:class:`~gpsea.analysis.predicate.genotype.VariantPredicate` to select a subset of variants of interest.


Examples
========


Use all variants
----------------

We can create a predicate to group individuals based on one or two variant alleles,
while considering all variants:

>>> from gpsea.analysis.predicate.genotype import allele_counting
>>> gt_predicate = allele_counting()
>>> gt_predicate.display_question()
'What is the genotype group: Monoallelic, Biallelic'


Use a variant subset
--------------------

If the cohort members include variants that are likely unrelated to the analysis, we can subset the variants
with a :class:`~gpsea.analysis.predicate.genotype.VariantPredicate`.
Here we will include the variants that have a functional annotation to a fictional transcript `NM_1234.5`:

First, we create the variant predicate:

>>> from gpsea.analysis.predicate.genotype import VariantPredicates
>>> affects_tx = VariantPredicates.transcript(tx_id="NM_1234.5")
>>> affects_tx.get_question()
'variant affects transcript NM_1234.5'

and then we wrap `affects_tx` with allele count predicate, to partition individuals
based on presence of one or two alleles that affect `NM_1234.5`:

>>> from gpsea.analysis.predicate.genotype import allele_counting
>>> gt_predicate = allele_counting(variant_predicate=affects_tx)
>>> gt_predicate.display_question()
'What is the genotype group: Monoallelic, Biallelic'
