.. _genotype-predicates:

###################
Genotype predicates
###################


A genotype predicate partitions the individuals based on their genotype.
In GPSEA, a predicate leverages information from one of the three areas:

* sex (Male vs. Female)
* disease diagnosis
* presence of variant(s) that meet certain inclusion criteria (e.g. a missense variant in hetezygous genotype)

Partitioning based on sex or disease diagnosis is relatively straightforward
(see :ref:`group-by-sex` and :ref:`group-by-diagnosis` for more details). 

Partitioning based on variants needs, however, more thought.
The process is split into two steps.
First, we count the alleles of the matching variants
(:ref:`variant-predicates`)
and then we interpret the count,
possibly including factors such as the expected mode of inheritance and sex,
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
