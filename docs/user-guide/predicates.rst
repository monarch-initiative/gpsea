.. _predicates:

==========
Predicates
==========

Genophenocorr uses predicates to test if variant, patient, or any tested item 
meets a condition. Based on the test results, the items are assigned into groups.

As described in the :class:`genophenocorr.analysis.predicate.PolyPredicate` API, 
the groups must be *exclusive* - the item can be assigned with one and only one group,
and *exhaustive* - the groups must cover all possible scenarios.

However, if the item cannot be assigned into any meaningful category, 
the predicate can return `None`, and the item will be omitted from the analysis.

The predicates can be chained to test for more complex conditions. 
For instance, "test if a patient has a missense or synonymous variant located in exon 6 of transcript `NM_013275.6`".

Let's demonstrate this on an example with a :class:`genophenocorr.analysis.predicate.genotype.VariantPredicate`.
We will load a cohort of 5 subjects with variants in *ANKRD11*, leading to KBG syndrome. 
The the clinical signs and symptoms of the subjects were encoded into HPO terms 
along with the pathogenic *ANKRD11* variant.

Let's load the phenopackets, as previously described in greater detail the :ref:`input-data` section.

First, we load HPO using HPO toolkit:

>>> import hpotk
>>> store = hpotk.configure_ontology_store()
>>> hpo = store.load_minimal_hpo(release='v2024-03-06')

then, we will configure the cohort creator:

>>> from genophenocorr.preprocessing import configure_caching_cohort_creator, load_phenopacket_folder
>>> cohort_creator = configure_caching_cohort_creator(hpo)

last, we will load the cohort from a directory with phenopackets:

>>> import os
>>> cohort_dir = os.path.join('docs', 'data', 'simple_cohort')
>>> cohort = load_phenopacket_folder(cohort_dir, cohort_creator) # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
Patients Created...
>>> len(cohort)
5

We loaded a cohort of 5 patients.

Let's use the variant ``16_89281397_89281397_G_C`` that corresponds 
to a pathogenic variant `VCV001029215.1 <https://www.ncbi.nlm.nih.gov/clinvar/variation/1029215/>`_ 
that replaces the tyrosine encoded by the 1715th codon of `NM_013275.6` with a premature stop codon: ``NM_013275.6(ANKRD11):c.5145C>G (p.Tyr1715Ter)``.

>>> variant_key_of_interest = '16_89281397_89281397_G_C'
>>> variant = cohort.get_variant_by_key(variant_key_of_interest)

We will use the variant to exemplify the predicate API.

Simple predicates
*****************

We can check that the variant overlaps with *ANKRD11*:

>>> from genophenocorr.analysis.predicate.genotype import VariantPredicates
>>> gene = VariantPredicates.gene('ANKRD11')
>>> gene.test(variant)
True

it overlaps with the *MANE* transcript:

>>> ankrd11_mane_tx_id = 'NM_013275.6'
>>> tx = VariantPredicates.transcript(ankrd11_mane_tx_id)
>>> tx.test(variant)
True

it in fact overlaps with the exon 9:

>>> exon9 = VariantPredicates.exon(exon=9, tx_id=ankrd11_mane_tx_id)
>>> exon9.test(variant)
True

and it is predicted to introduce a premature termination codon to the MANE transcript:

>>> from genophenocorr.model import VariantEffect
>>> nonsense = VariantPredicates.variant_effect(VariantEffect.STOP_GAINED, tx_id=ankrd11_mane_tx_id)
>>> nonsense.test(variant)
True

See :class:`genophenocorr.analysis.predicate.genotype.VariantPredicates` 
for more info on the predicates available off the shelf.


Compound predicates
*******************

The simple predicates can be combined to test for more elaborate conditions.
For instance, we can test if the variant meets *any* or several conditions:

>>> missense = VariantPredicates.variant_effect(VariantEffect.MISSENSE_VARIANT, tx_id=ankrd11_mane_tx_id)
>>> missense_or_nonsense = missense | nonsense
>>> missense_or_nonsense.test(variant)
True

or *all* conditions:

>>> nonsense_and_exon9 = nonsense & exon9
>>> nonsense_and_exon9.test(variant)
True

The `VariantPredicate` overloads Python ``&`` (AND) and ``|`` (OR) operators to build a compound predicate from lower level building blocks.

Therefore, there is nothing that prevents us to combine the predicates into multi-level tests, 
such as testing if the variant is a *"chromosomal deletion" or a deletion which removes at least 50 bp*:

>>> from genophenocorr.model import VariantClass
>>> chromosomal_deletion = "SO:1000029"
>>> predicate = VariantPredicates.structural_type(chromosomal_deletion) | (VariantPredicates.variant_class(VariantClass.DEL) & VariantPredicates.change_length("<=", -50))
>>> predicate.get_question()
'(structural type is SO:1000029 OR (variant class is DEL AND change length is <=-50))'


That's it for predicates. Please see :class:`genophenocorr.analysis.predicate.genotype.VariantPredicates` 
and :class:`genophenocorr.analysis.predicate.genotype.ProteinPredicates` 
for a comprehensive list of the predicates available off the shelf.

Please open an issue on our `GitHub tracker <https://github.com/monarch-initiative/genophenocorr/issues>`_ if a predicate seems to be missing.
