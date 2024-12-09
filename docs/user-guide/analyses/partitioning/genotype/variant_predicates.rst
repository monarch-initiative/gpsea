.. _variant-predicates:


==================
Variant Predicates
==================


Variant predicate is a core component to partition a cohort using the genomic variants identified in the cohort members.
A :class:`~gpsea.analysis.predicate.genotype.VariantPredicate`
tests if a :class:`~gpsea.model.Variant` meets the inclusion criteria.
For instance, a predicate can test if a variant is a deletion,
leads to a missense change, or overlaps with a protein domain.

An array of variant predicates is available as static methods
of the :class:`~gpsea.analysis.predicate.VariantPredicates` class.

The predicates operate on several lines of information:

+------------------------+-------------------------------------------------------------------------------------------------+
| Information            | Example                                                                                         |
+========================+=================================================================================================+
| Allele                 | variant is e.g. a deletion of >50 bases                                                         |
+------------------------+-------------------------------------------------------------------------------------------------+
| Functional annotation  | variant leads to a missense change or affects *n*-th exon of a transcript                       |
+------------------------+-------------------------------------------------------------------------------------------------+
| Protein data           | variant is located in a region encoding a protein domain, protein feature type                  |
+------------------------+-------------------------------------------------------------------------------------------------+
| Genome                 | overlap with a genomic region of interest                                                       |
+------------------------+-------------------------------------------------------------------------------------------------+


The scope of the builtin predicates is fairly narrow
and likely insufficient for real-life analyses.
The predicates can, however, be chained into a compound predicate to test conditions,
such as "variant is a missense or synonymous variant located in exon 6 of `NM_013275.6`".


********
Examples
********

Here we show examples of several simple variant predicates and 
and how to combine them to test a complex condition.


Load cohort
-----------

Let's start by loading a :class:`~gpsea.model.Cohort`
of 19 individuals with mutations in *RERE* leading to Holt-Oram syndrome.
The cohort was prepared from phenopackets as described in :ref:`create-a-cohort` section,
and then serialized as
a `JSON file <https://github.com/monarch-initiative/gpsea/tree/main/docs/cohort-data/RERE.0.1.20.json>`_
following the instructions in :ref:`cohort-persistence` section.

.. 
   Prepare the JSON file by running the tests in `tests/tests/test_generate_doc_cohorts.py`.

>>> import json
>>> from gpsea.io import GpseaJSONDecoder
>>> fpath_cohort_json = 'docs/cohort-data/RERE.0.1.20.json'
>>> with open(fpath_cohort_json) as fh:
...     cohort = json.load(fh, cls=GpseaJSONDecoder)
>>> len(cohort)
19


Some individuals were found to harbor the variant ``1_8358231_8358231_T_C`` that corresponds 
to a pathogenic mutation `VCV000522858.5 <https://www.ncbi.nlm.nih.gov/clinvar/variation/522858/>`_ 
that replaces the histidine encoded by the 1435th codon of `NM_001042681.2` with arginine: ``NM_001042681.2(RERE):c.4304A>G (p.His1435Arg)``.
We can retrieve the variant by querying the cohort by the variant key:

>>> variant_key_of_interest = '1_8358231_8358231_T_C'
>>> variant = cohort.get_variant_by_key(variant_key_of_interest)


Builtin predicates
------------------

Let's use builtin predicates to verify the properties of the variant ``1_8358231_8358231_T_C``.

We can check that the variant overlaps with *RERE*

>>> from gpsea.analysis.predicate.genotype import VariantPredicates
>>> gene = VariantPredicates.gene('RERE')
>>> gene.test(variant)
True

it overlaps with the *MANE* transcript

>>> rere_mane_tx_id = 'NM_001042681.2'
>>> tx = VariantPredicates.transcript(rere_mane_tx_id)
>>> tx.test(variant)
True

it in fact overlaps with the exon 20,

>>> exon20 = VariantPredicates.exon(exon=20, tx_id=rere_mane_tx_id)
>>> exon20.test(variant)
True

and leads to a missense mutation with respect to the MANE transcript

>>> from gpsea.model import VariantEffect
>>> missense = VariantPredicates.variant_effect(VariantEffect.MISSENSE_VARIANT, tx_id=rere_mane_tx_id)
>>> missense.test(variant)
True

See :class:`~gpsea.analysis.predicate.genotype.VariantPredicates`
for a complete list of the builtin predicates.


Predicate chain
---------------

Using the builtin predicates, we can build a logical chain to test complex conditions.
For instance, we can test if the variant meets any of several conditions:

>>> nonsense = VariantPredicates.variant_effect(VariantEffect.STOP_GAINED, tx_id=rere_mane_tx_id)
>>> missense_or_nonsense = missense | nonsense
>>> missense_or_nonsense.test(variant)
True

or *all* conditions:

>>> missense_and_exon20 = missense & exon20
>>> missense_and_exon20.test(variant)
True

All variant predicates overload Python ``&`` (AND) and ``|`` (OR) operators, to allow chaining.

Therefore, there is nothing that prevents us to combine the predicates into multi-level tests, 
e.g. to test if the variant is a *"chromosomal deletion" or a deletion which removes at least 50 bp*:

>>> from gpsea.model import VariantClass
>>> chromosomal_deletion = "SO:1000029"
>>> predicate = VariantPredicates.structural_type(chromosomal_deletion) | (VariantPredicates.variant_class(VariantClass.DEL) & VariantPredicates.change_length("<=", -50))
>>> predicate.description
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

However, this is clearly much better implemented by a logical *not* of a "is frameshift" predicate.

Therefore, all variant predicates implement *logical inversion* 
which corresponds to Python's ``~`` operator (tilde),
and results in an inverted predicate.

This is how we can use the predicate inversion to build the predicate for non-frameshift deletions:

>>> non_frameshift_del = ~VariantPredicates.variant_effect(VariantEffect.FRAMESHIFT_VARIANT, tx_id=rere_mane_tx_id) & VariantPredicates.variant_class(VariantClass.DEL)
>>> non_frameshift_del.description
'(NOT FRAMESHIFT_VARIANT on NM_001042681.2 AND variant class is DEL)'

Note the presence of a tilde ``~`` before the variant effect predicate and resulting ``NOT`` in the predicate question.


**********
Need more?
**********

The builtin predicates should cover majority of use cases.
However, if a predicate seems to be missing,
feel free to submit an issue in our
`GitHub tracker <https://github.com/monarch-initiative/gpsea/issues>`_,
or to implement a custom predicate
by extending the :class:`~gpsea.analysis.predicate.genotype.VariantPredicate` class 😎.



The variant predicate offers a flexible API for testing if variants meet a condition.
However, the genotype phenotype correlations are done on the individual level
and the variant predicates are used as a component of the genotype predicate.
The next sections show how to use variant predicates to assign individuals into groups.
