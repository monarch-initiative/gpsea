.. _male_female_predicate:

Partition by the sex of the individual
======================================

It is easy to investigate the phenotypic differences between females and males.
The :func:`~gpsea.analysis.predicate.genotype.sex_predicate` provides a predicate
for partitioning based on the sex of the individual:

>>> from gpsea.analysis.predicate.genotype import sex_predicate
>>> gt_predicate = sex_predicate()
>>> gt_predicate.display_question()
'Sex of the individual: FEMALE, MALE'

The individuals with :class:`~gpsea.model.Sex.UNKNOWN_SEX` will be omitted from the analysis.

Note that we have implemented this predicate as a genotype predicate, because it is used in 
place of other genotype predicates. Currently, it is not possible to compare the distribution of genotypes across sexes.



