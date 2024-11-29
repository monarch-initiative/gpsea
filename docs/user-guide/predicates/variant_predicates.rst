.. _variant_predicates:

==================
Variant Predicates
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


As a rule of thumb, the predicates for testing basic conditions are available off the shelf,
and they can be used as building block for testing for more complex conditions,
such as testing if the variant is "a missense or synonymous variant located in exon 6 of transcript `NM_013275.6`".

Let's demonstrate the variant predicate usage on a few examples.


Load cohort
-----------

For the purpose of this example, we will load a :class:`~gpsea.model.Cohort`
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

However, this is clearly tedious and it would be much better implemented 
by a simple logical not of a predicate for a frameshift variant effect.

To support this, `VariantPredicate` implements *logical inversion* 
which corresponds to Python's ``~`` operator (tilde), to wrap
the underlying predicate and to invert its test result.

This is how we can use the predicate inversion to build the predicate for non-frameshift deletions:

>>> non_frameshift_del = ~VariantPredicates.variant_effect(VariantEffect.FRAMESHIFT_VARIANT, tx_id=rere_mane_tx_id) & VariantPredicates.variant_class(VariantClass.DEL)
>>> non_frameshift_del.description
'(NOT FRAMESHIFT_VARIANT on NM_001042681.2 AND variant class is DEL)'

Note the presence of a tilde ``~`` before the variant effect predicate and resulting ``NOT`` in the predicate question.

The variant predicate offers a flexible API for testing if variants meet a condition.
However, the genotype phenotype correlations are done on the individual level
and the variant predicates are used as a component of the genotype predicate.
The next sections show how to use variant predicates to assign individuals into groups.
