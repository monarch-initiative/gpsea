.. _genotype-predicates:

###################
Genotype predicates
###################


A genotype predicate partitions the individuals based on their genotype.
In GPSEA, genotype predicates leverage information from one of the three areas:

* Sex
* Disease diagnosis
* Presence of variant(s) that meet certain inclusion criteria (e.g. a missense variant in heterozygous genotype)

Partitioning based on sex or disease diagnosis is relatively straightforward - the individuals
are assigned by the biological sex or presence of a specific diagnosis.
See :ref:`group-by-sex` and :ref:`group-by-diagnosis` for more details. 

Partitioning based on variants is, however, much more flexible,
to support the analysis of the broad spectrum of pathomechanisms
that have been shown to lead to genetic diseases.
In general, we first create one or more variant predicates to select the variants of interest (:ref:`variant-predicates`),
then count and interpret the number of alleles in the cohort members,
including factors such as the expected mode of inheritance and sex,
to assign the tested individual into a genotype group
(:ref:`variant-category` or :ref:`allele-count`).

.. toctree::
  :maxdepth: 1
  :caption: Contents:

  sex
  diagnosis
  variant_predicates
  variant_category
  allele_count
