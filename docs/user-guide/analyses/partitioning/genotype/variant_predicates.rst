.. _variant-predicates:


==================
Variant Predicates
==================


Variant predicate is a core component to classify a cohort using the genomic variants identified in the cohort members.
A :class:`~gpsea.analysis.predicate.VariantPredicate`
tests if a :class:`~gpsea.model.Variant` meets the inclusion criteria.
For instance, a predicate can test if a variant is a deletion,
leads to a missense change, or overlaps with a protein domain.

An array of builtin variant predicates is available as functions
of the :mod:`gpsea.analysis.predicate` module.

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


The scope of the builtin predicates is fairly narrow
and likely insufficient for real-life analyses.
However, several predicates can be "chained" into a compound predicate using a boolean logic,
to achive more expressivity for testing complex conditions,
such as "variant is a missense or synonymous variant located in exon 6 of `NM_013275.6`".


********
Examples
********

Here we show how to use the builtin predicates for simple tests
and how to build a compound predicate from the builtin predicates,
for testing complex conditions.


Load cohort
===========

Let's start by loading a :class:`~gpsea.model.Cohort`
of 19 individuals with mutations in *RERE* leading to Holt-Oram syndrome.
The cohort was prepared from phenopackets as described in :ref:`create-a-cohort` section,
and then serialized as
a `JSON file <https://github.com/P2GX/gpsea/tree/main/docs/cohort-data/RERE.0.1.20.json>`_
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
==================

Let's use builtin predicates to verify the properties of the variant ``1_8358231_8358231_T_C``.

We can check that the variant overlaps with *RERE*

>>> import gpsea.analysis.predicate as vp
>>> gene = vp.gene('RERE')
>>> gene.test(variant)
True

it overlaps with the *MANE* transcript

>>> rere_mane_tx_id = 'NM_001042681.2'
>>> tx = vp.transcript(rere_mane_tx_id)
>>> tx.test(variant)
True

it in fact overlaps with the exon 20,

>>> exon20 = vp.exon(exon=20, tx_id=rere_mane_tx_id)
>>> exon20.test(variant)
True

and leads to a missense mutation with respect to the MANE transcript

>>> from gpsea.model import VariantEffect
>>> missense = vp.variant_effect(VariantEffect.MISSENSE_VARIANT, tx_id=rere_mane_tx_id)
>>> missense.test(variant)
True

See the :mod:`gpsea.analysis.predicate` module
for a complete list of the builtin predicates.


Compound predicates
===================

A compound predicate for testing complex conditions can be built from two or more predicates.
For instance, we can test if the variant meets any of several conditions:

>>> import gpsea.analysis.predicate as vp
>>> nonsense = vp.variant_effect(VariantEffect.STOP_GAINED, tx_id=rere_mane_tx_id)
>>> missense_or_nonsense = missense | nonsense
>>> missense_or_nonsense.test(variant)
True

or *all* conditions:

>>> missense_and_exon20 = missense & exon20
>>> missense_and_exon20.test(variant)
True

All variant predicates overload Python ``&`` (AND) and ``|`` (OR) operators,
to combine a predicate pair into a compound predicate.

.. note::

  Combining three or or more predicates can be achieved with :func:`~gpsea.analysis.allof`
  and :func:`~gpsea.analysis.anyof` functions.

Therefore, there is nothing that prevents us to combine the predicates into multi-level tests, 
e.g. to test if the variant is a *"chromosomal deletion" or a deletion which removes at least 50 bp*:

>>> from gpsea.model import VariantClass
>>> chromosomal_deletion = "SO:1000029"
>>> predicate = vp.structural_type(chromosomal_deletion) | (vp.variant_class(VariantClass.DEL) & vp.change_length("<=", -50))
>>> predicate.description
'(structural type is SO:1000029 OR (variant class is DEL AND change length <= -50))'


Inverting conditions
====================

Sometimes we may want to test the variant for a condition that must *not* be met.
For instance, we may want to test if the variant is a deletion 
that is *not* predicted to shift the transcript reading frame.
One of doing this would be to build a compound predicates 
for all variant effects except of :class:`~gpsea.model.VariantEffect.FRAMESHIFT_VARIANT`:

>>> non_frameshift_effects = (
...   VariantEffect.SYNONYMOUS_VARIANT, VariantEffect.MISSENSE_VARIANT, VariantEffect.INTRON_VARIANT,
...   # and many more effects..
... )
>>> non_frameshift_predicate = vp.allof(vp.variant_effect(eff, tx_id=rere_mane_tx_id) for eff in non_frameshift_effects)

However, this is clearly much better implemented by a logical *not* of a "is frameshift" predicate.

Therefore, all variant predicates implement *logical inversion* 
which corresponds to Python's ``~`` operator (tilde),
and results in an inverted predicate.

This is how we can use the predicate inversion to build the predicate for non-frameshift deletions:

>>> non_frameshift_del = ~vp.variant_effect(VariantEffect.FRAMESHIFT_VARIANT, tx_id=rere_mane_tx_id) & vp.variant_class(VariantClass.DEL)
>>> non_frameshift_del.description
'(NOT FRAMESHIFT_VARIANT on NM_001042681.2 AND variant class is DEL)'

Note the presence of a tilde ``~`` before the variant effect predicate and resulting ``NOT`` in the predicate question.


**********
Need more?
**********

The builtin predicates should cover majority of use cases.
However, if a predicate seems to be missing,
feel free to submit an issue in our
`GitHub tracker <https://github.com/P2GX/gpsea/issues>`_,
or implement your own predicate by following the :ref:`custom-variant-predicate`
guide.


****
Next
****

The variant predicate offers a flexible API for testing if variants meet a condition.
However, the genotype phenotype correlations are studied on the level of individuals.
As described in :ref:`genotype-classifiers`, GPSEA uses the :class:`~gpsea.analysis.clf.GenotypeClassifier` API
to assign individuals into non-overlapping classes. Variant predicates are essential for creating such classifier.
We explain the details in the following sections.
