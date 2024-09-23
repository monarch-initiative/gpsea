.. _genotype-predicates:

===================
Genotype Predicates
===================


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


.. toctree::
  :maxdepth: 1
  :caption: Contents:

  variant_predicates
  protein_domains
  mode_of_inheritance_predicate
  male_female_predicate
  diagnosis_predicate
  groups_predicate



