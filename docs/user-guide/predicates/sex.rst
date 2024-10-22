.. _group-by-sex:

============
Group by sex
============

It is easy to investigate the differences between males and females.
The :func:`~gpsea.analysis.predicate.genotype.sex_predicate` partitions
the individuals based on their :class:`~gpsea.model.Sex`:

>>> from gpsea.analysis.predicate.genotype import sex_predicate
>>> gt_predicate = sex_predicate()
>>> gt_predicate.display_question()
'Sex of the individual: FEMALE, MALE'

The individuals with :class:`~gpsea.model.Sex.UNKNOWN_SEX` will be omitted from the analysis.

Note that we implemented this predicate as a genotype predicate.
Currently, it is not possible to compare the distribution of genotypes across sexes.

