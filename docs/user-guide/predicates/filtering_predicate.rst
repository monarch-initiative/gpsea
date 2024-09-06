.. _filtering-predicate:


===================
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
'What is the genotype group?: HOM_REF, HET, BIALLELIC_ALT'

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
'What is the genotype group?: HET, BIALLELIC_ALT'