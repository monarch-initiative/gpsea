.. _predicates:

==========
Predicates
==========

GPSEA uses predicates to test if variant, patient, or any tested item 
meets a condition. Based on the test results, the items are assigned into groups.

As described in the :class:`~gpsea.analysis.predicate.PolyPredicate` API, 
the groups must be *exclusive* - the item can be assigned with one and only one group,
and *exhaustive* - the groups must cover all possible scenarios.

However, if the item cannot be assigned into any meaningful category, 
the predicate can return `None`, and the item will be omitted from the analysis.

The predicates can be chained to test for more complex conditions. 
For instance, "test if a patient has a missense or synonymous variant located in exon 6 of transcript `NM_013275.6`".

Let's demonstrate this on an example with a :class:`~gpsea.analysis.predicate.genotype.VariantPredicate`.
We will load a cohort of 5 subjects with variants in *ANKRD11*, leading to KBG syndrome. 
The the clinical signs and symptoms of the subjects were encoded into HPO terms 
along with the pathogenic *ANKRD11* variant.

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
This time, however, we will load 19 individuals with mutations in *RERE* gene:

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

Simple predicates
*****************

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


Compound predicates
*******************

The simple predicates can be combined to test for more elaborate conditions.
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


Inverted predicate
******************

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


That's it for predicates. Please see :class:`~gpsea.analysis.predicate.genotype.VariantPredicates` 
and :class:`~gpsea.analysis.predicate.genotype.ProteinPredicates` 
for a comprehensive list of the predicates available off the shelf.

Please open an issue on our `GitHub tracker <https://github.com/monarch-initiative/gpsea/issues>`_ if a predicate seems to be missing.
