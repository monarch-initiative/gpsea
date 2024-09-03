.. _predicates:

==========
Predicates
==========

Searching for genotype-phenotype associations usually requires to partition
the individuals into several discrete groups to allow testing for the inter-group differences.
GPSEA reflects these requirements with its predicate API.
Perhaps unsurprisingly, a predicate must be capable of partitioning the individuals into two or more groups.
The groups must be *exclusive* - each individual must be assigned at most into one group.
Moreover, the groups should be *exhaustive* and cover maximum of the possible states.
However, the predicate is allowed to return `None` if the individual cannot be assigned.
In result, the individual will be omitted from the downstream analysis.

Predicates serve both *genotype* and *phenotype* prongs of the analysis.
Genotype predicates (:class:`~gpsea.analysis.predicate.genotype.GenotypePolyPredicate`)
assign the :class:`~gpsea.model.Patient`
into a group (mostly) based on the variant information, while the
phenotype predicates (:class:`~gpsea.analysis.predicate.phenotype.PhenotypePolyPredicate`)
use the HPO terms to assign a group.

It is literally impossible to use GPSEA without the predicates
because all analyses need at least one predicate (typically a *genotype* predicate).
Luckily, the purpose of this guide is to show all that is to know about predicates.
We will first discuss the genotype predicates and end with phenotype predicates.

.. _genotype-predicates:

*******************
Genotype predicates
*******************

A genotype predicate seeks to divide the individuals along an axis that is orthogonal to phenotypes.
Typically, this includes using the genotype data, such as presence of a missense variant
in a heterozygous genotype. However, other categorical variables,
such as diagnoses (TODO - add link to disease predicate) or cluster ids can also be used.

The genotype predicates test the individual for a presence of variants that meet certain inclusion criteria.
The testing is done in two steps. First, we count the alleles
of the matching variants and then we interpret the count, possibly including factors
such as the expected mode of inheritance and sex, to assign the individual into a group.
Finding the matching variants is what
the :class:`~gpsea.analysis.predicate.genotype.VariantPredicate` is all about.


TODO: wordsmith
We must first create the variant predicate and then wrap it in genotype predicate.


Variant predicates
==================

GPSEA uses the :class:`~gpsea.analysis.predicate.genotype.VariantPredicate` class
to test if a :class:`~gpsea.model.Variant` meets the inclusion criteria.
The variant predicate can leverage multiple primary data:

+------------------------+-------------------------------------------------------------------------------------------------+
| Primary data source    |   Example                                                                                       |
+========================+=================================================================================================+
| Allele                 | the variant being a deletion or a single nucleotide variant (SNV)                               |
+------------------------+-------------------------------------------------------------------------------------------------+
| Genome                 | overlaps of a target genomic region                                                             |
+------------------------+-------------------------------------------------------------------------------------------------+
| Functional annotation  | variant is predicted to lead to a missense change or affect an exon of certain transcript       |
+------------------------+-------------------------------------------------------------------------------------------------+
| Protein data           | variant is located in a region encoding a protein domain, protein feature type                  |
+------------------------+-------------------------------------------------------------------------------------------------+

which demands a considerable amount of flexibility for creating the predicate.

As a rule of thumb, the predicates for testing basic conditions are available off the shelf,
and they can be used as building block for testing for more complex conditions,
such as testing if the variant is "a missense or synonymous variant located in exon 6 of transcript `NM_013275.6`".

Let's demonstrate this on few examples.
We will load a cohort of 19 subjects with variants in *RERE*,
leading to `Holt-Oram syndrome MIM:142900 <https://omim.org/entry/142900>`_.
The the clinical signs and symptoms of the subjects were encoded into HPO terms
along with the pathogenic *RERE* variant.

Let's load the phenopackets, as previously described in greater detail the :ref:`input-data` section.
Briefly, we first load HPO:

>>> import hpotk
>>> store = hpotk.configure_ontology_store()
>>> hpo = store.load_minimal_hpo(release='v2024-07-01')

then, we configure the cohort creator:

>>> from gpsea.preprocessing import configure_caching_cohort_creator
>>> cohort_creator = configure_caching_cohort_creator(hpo)

which we use to create a :class:`~gpsea.model.Cohort` from a bunch of phenopackets
from the release `0.1.18` of `Phenopacket Store <https://github.com/monarch-initiative/phenopacket-store>`_.
We will load 19 individuals with mutations in *RERE* gene:

>>> from ppktstore.registry import configure_phenopacket_registry
>>> registry = configure_phenopacket_registry()
>>> with registry.open_phenopacket_store(release='0.1.18') as ps:
...     phenopackets = tuple(ps.iter_cohort_phenopackets('RERE'))
>>> len(phenopackets)
19

and we will convert the phenopacket into a :class:`~gpsea.model.Cohort`:

>>> from gpsea.preprocessing import load_phenopackets
>>> cohort, _ = load_phenopackets(phenopackets, cohort_creator)  # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
Patients Created: ...

To demonstrate the predicate API, we will use the variant ``1_8358231_8358231_T_C`` that corresponds 
to a pathogenic variant `VCV000522858.5 <https://www.ncbi.nlm.nih.gov/clinvar/variation/522858/>`_ 
that replaces the histidine encoded by the 1435th codon of `NM_001042681.2` with arginine: ``NM_001042681.2(RERE):c.4304A>G (p.His1435Arg)``.

>>> variant_key_of_interest = '1_8358231_8358231_T_C'
>>> variant = cohort.get_variant_by_key(variant_key_of_interest)

Building blocks
---------------

We can check that the variant overlaps with *RERE*:

>>> from gpsea.analysis.predicate.genotype import VariantPredicates
>>> gene = VariantPredicates.gene('RERE')
>>> gene.test(variant)
True

it overlaps with the *MANE* transcript:

>>> rere_mane_tx_id = 'NM_001042681.2'
>>> tx = VariantPredicates.transcript(rere_mane_tx_id)
>>> tx.test(variant)
True

it in fact overlaps with the exon 20:

>>> exon20 = VariantPredicates.exon(exon=20, tx_id=rere_mane_tx_id)
>>> exon20.test(variant)
True

and leads to a missense mutation with respect to the MANE transcript:

>>> from gpsea.model import VariantEffect
>>> missense = VariantPredicates.variant_effect(VariantEffect.MISSENSE_VARIANT, tx_id=rere_mane_tx_id)
>>> missense.test(variant)
True

See :class:`~gpsea.analysis.predicate.genotype.VariantPredicates` 
for more info on the predicates available off the shelf.


Complex conditions
------------------

We can combine the building blocks to test for more elaborate conditions.
For instance, we can test if the variant meets *any* or several conditions:

>>> nonsense = VariantPredicates.variant_effect(VariantEffect.STOP_GAINED, tx_id=rere_mane_tx_id)
>>> missense_or_nonsense = missense | nonsense
>>> missense_or_nonsense.test(variant)
True

or *all* conditions:

>>> missense_and_exon20 = missense & exon20
>>> missense_and_exon20.test(variant)
True

The `VariantPredicate` overloads Python ``&`` (AND) and ``|`` (OR) operators to build a compound predicate from lower level building blocks.

Therefore, there is nothing that prevents us to combine the predicates into multi-level tests, 
such as testing if the variant is a *"chromosomal deletion" or a deletion which removes at least 50 bp*:

>>> from gpsea.model import VariantClass
>>> chromosomal_deletion = "SO:1000029"
>>> predicate = VariantPredicates.structural_type(chromosomal_deletion) | (VariantPredicates.variant_class(VariantClass.DEL) & VariantPredicates.change_length("<=", -50))
>>> predicate.get_question()
'(structural type is SO:1000029 OR (variant class is DEL AND change length <= -50))'


Inverting conditions
--------------------

Sometimes we may want to test the variant for a condition that must *not* be met.
For instance, we may want to test if the variant is a deletion 
that is *not* predicted to shift the transcript reading frame.
One of doing this would be to build a compound predicates 
for all variant effects except of :class:`~gpsea.model.VariantEffect.FRAMESHIFT_VARIANT`:

>>> non_frameshift_effects = (
...   VariantEffect.SYNONYMOUS_VARIANT, VariantEffect.MISSENSE_VARIANT, VariantEffect.INTRON_VARIANT,
...   # and many more effects..
... )
>>> non_frameshift_predicate = VariantPredicates.all(VariantPredicates.variant_effect(eff, tx_id=rere_mane_tx_id) for eff in non_frameshift_effects)

However, this is clearly tedious and it would be much better implemented 
by a simple logical not of a predicate for a frameshift variant effect.

To support this, `VariantPredicate` implements *logical inversion* 
which corresponds to Python's ``~`` operator (tilde), to wrap
the underlying predicate and to invert its test result.

This is how we can use the predicate inversion to build the predicate for non-frameshift deletions:

>>> non_frameshift_del = ~VariantPredicates.variant_effect(VariantEffect.FRAMESHIFT_VARIANT, tx_id=rere_mane_tx_id) & VariantPredicates.variant_class(VariantClass.DEL)
>>> non_frameshift_del.get_question()
'(NOT FRAMESHIFT_VARIANT on NM_001042681.2 AND variant class is DEL)'

Note the presence of a tilde ``~`` before the variant effect predicate and resulting ``NOT`` in the predicate question.

The variant predicate offers a flexible API for testing if variants meet a condition.
However, the genotype phenotype correlations are done on the individual level
and the variant predicates are used as a component of the genotype predicate.
The next sections show how to use variant predicates to assign individuals into groups.


Mode of inheritance predicate
=============================

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
'Which genotype group does the patient fit in: HOM_REF, HET, BIALLELIC_ALT'

The `gt_predicate` can be used in downstream analysis, such as in :class:


.. _filtering-predicate:

Filtering predicate
===================

Sometimes a predicate can bin individuals into more genotype groups than necessary and there may be need
to consider only a subset of the groups. A `GenotypePolyPredicate`
created by :class:`~gpsea.analysis.predicate.genotype.filtering_predicate` can retain only a subset
of the target categorizations of interest.

Example
-------

Let's suppose we want test the genotype-phenotype association between variants
that lead to frameshift or a stop gain in a fictional transcript `NM_1234.5`,
and we are specifically interested in comparing the heterozygous variants
in a biallelic alternative allele genotypes (homozygous alternate and compound heterozygous).

First, we set up a :class:`~gpsea.analysis.predicate.genotype.VariantPredicate`
for testing if a variant introduces a premature stop codon or leads to the shift of the reading frame:

>>> from gpsea.model import VariantEffect
>>> from gpsea.analysis.predicate.genotype import VariantPredicates
>>> tx_id = 'NM_1234.5'
>>> is_frameshift_or_stop_gain = VariantPredicates.variant_effect(VariantEffect.FRAMESHIFT_VARIANT, tx_id) \
...     | VariantPredicates.variant_effect(VariantEffect.STOP_GAINED, tx_id)
>>> is_frameshift_or_stop_gain.get_question()
'(FRAMESHIFT_VARIANT on NM_1234.5 OR STOP_GAINED on NM_1234.5)'

Then, we create :class:`~gpsea.analysis.predicate.genotype.ModeOfInheritancePredicate.autosomal_recessive`
to bin according to a genotype group:

>>> from gpsea.analysis.predicate.genotype import ModeOfInheritancePredicate
>>> gt_predicate = ModeOfInheritancePredicate.autosomal_recessive(is_frameshift_or_stop_gain)
>>> gt_predicate.display_question()
'Which genotype group does the patient fit in: HOM_REF, HET, BIALLELIC_ALT'

We see that the `gt_predicate` bins the patients into three groups:

>>> cats = gt_predicate.get_categorizations()
>>> cats
(Categorization(category=HOM_REF), Categorization(category=HET), Categorization(category=BIALLELIC_ALT))

We wrap the categorizations of interest along with the `gt_predicate` by the `filtering_predicate` function,
and we will get a :class:`~gpsea.analysis.predicate.genotype.GenotypePolyPredicate`
that includes only the categories of interest:

>>> from gpsea.analysis.predicate.genotype import filtering_predicate
>>> fgt_predicate = filtering_predicate(
...     predicate=gt_predicate,
...     targets=(cats[1], cats[2]),
... )
>>> fgt_predicate.display_question()
'Which genotype group does the patient fit in: HET, BIALLELIC_ALT'


.. _groups-predicate:

Groups predicate
================

Sometimes, all we want is to compare if there is a difference between individuals
who include one or more alleles of variant $X$ vs. individuals with variants $Y$,
vs. individuals with variants $Z$, where $X$, $Y$ and $Z$ are variant predicates.
We can do this with a *groups* predicate.

The :func:`~gpsea.analysis.predicate.genotype.groups_predicate`
takes *n* variant predicates and *n* group labels, and it will assign the patients
into the respective groups if one or more matching allele is found.
However, only one predicate is allowed to return a non-zero allele count.
Otherwise, the patient is assigned with ``None`` and excluded from the analysis.

Example
-------

Here we show how to build a :class:`~gpsea.analysis.predicate.genotype.GenotypePolyPredicate`
for testing if the individual has at least one missense vs. frameshift vs. synonymous variant.

>>> from gpsea.model import VariantEffect
>>> from gpsea.analysis.predicate.genotype import VariantPredicates, groups_predicate
>>> tx_id = 'NM_1234.5'
>>> gt_predicate = groups_predicate(
...     predicates=(
...         VariantPredicates.variant_effect(VariantEffect.MISSENSE_VARIANT, tx_id),
...         VariantPredicates.variant_effect(VariantEffect.FRAMESHIFT_VARIANT, tx_id),
...         VariantPredicates.variant_effect(VariantEffect.SYNONYMOUS_VARIANT, tx_id),
...     ),
...     group_names=('Missense', 'Frameshift', 'Synonymous'),
... )
>>> gt_predicate.display_question()
'Genotype group: Missense, Frameshift, Synonymous'


.. _phenotype-predicates:

********************
Phenotype predicates
********************

The phenotype predicate assigns the individual into a group with respect to tested phenotype.
Typically, the phenotype corresponds to a clinical sign or symptom encoded into an HPO term.


Propagating phenotype predicate
===============================

When testing for presence or absence of an HPO term, the propagating phenotype predicate
leverages the :ref:`true-path-rule` to take advantage of the HPO hierarchy.
In result, an individual annotated with a term is implicitly annotated with all its ancestors.
For instance, an individual annotated with `Ectopia lentis <https://hpo.jax.org/browse/term/HP:0001083>`_
is also annotated with `Abnormal lens morphology <https://hpo.jax.org/browse/term/HP:0000517>`_,
`Abnormal anterior eye segment morphology <https://hpo.jax.org/browse/term/HP:0004328>`_,
`Abnormal eye morphology <https://hpo.jax.org/browse/term/HP:0012372>`_, ...

Similarly, all descendants of a term, whose presence was specifically excluded in an individual,
are implicitly excluded.

:class:`~gpsea.analysis.predicate.phenotype.PropagatingPhenotypePredicate` implements this logic.

Example
-------

Here we show how to set up :class:`~gpsea.analysis.predicate.phenotype.PropagatingPhenotypePredicate`
to test for a presence of `Abnormal lens morphology <https://hpo.jax.org/browse/term/HP:0000517>`_.


>>> from gpsea.analysis.predicate.phenotype import PropagatingPhenotypePredicate
>>> query = hpotk.TermId.from_curie('HP:0000517')
>>> pheno_predicate = PropagatingPhenotypePredicate(
...     hpo=hpo,
...     query=query,
... )
>>> pheno_predicate.display_question()
'Is Abnormal lens morphology present in the patient: Yes, No'


TODO: explain ``missing_implies_phenotype_excluded``


Predicates for all cohort phenotypes
====================================

Constructing phenotype predicates for all HPO terms of a cohort sounds a bit tedious.
The :func:`~gpsea.analysis.predicate.phenotype.prepare_predicates_for_terms_of_interest`
function cuts down the tedium:

>>> from gpsea.analysis.predicate.phenotype import prepare_predicates_for_terms_of_interest
>>> pheno_predicates = prepare_predicates_for_terms_of_interest(
...     cohort=cohort,
...     hpo=hpo,
... )
>>> len(pheno_predicates)
301

and prepares predicates for testing 301 HPO terms of the *RERE* cohort.


*******
Gallery
*******

Here we show examples of predicates used in some of our analyses.

TODO


**********
Need more?
**********

Please see :class:`~gpsea.analysis.predicate.genotype.VariantPredicates` 
and :class:`~gpsea.analysis.predicate.genotype.ProteinPredicates` 
for a list of the predicates available off the shelf.

However, feel free to open an issue on our `GitHub tracker <https://github.com/monarch-initiative/gpsea/issues>`_
if a predicate seems to be missing.
