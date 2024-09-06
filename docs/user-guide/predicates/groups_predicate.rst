.. _groups-predicate:

================
Groups Predicate
================



Sometimes, all we want is to compare if there is a difference between individuals
who include one or more alleles of variant `X` vs. individuals with variants `Y`,
vs. individuals with variants `Z`, where `X`, `Y` and `Z` are variant predicates.
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



