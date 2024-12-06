.. _group-by-diagnosis:

==================
Group by diagnosis
==================

We can easily compare individuals diagnosed with different diseases.
:func:`~gpsea.analysis.predicate.genotype.diagnosis_predicate` groups the individuals
based on the disease diagnoses present in the individuals.

>>> from gpsea.analysis.predicate.genotype import diagnosis_predicate
>>> gt_predicate = diagnosis_predicate(
...     diagnoses=('OMIM:154700', 'OMIM:129600'),
...     labels=('Marfan syndrome', 'Ectopia lentis, familial'),
... )
>>> gt_predicate.group_labels
('OMIM:154700', 'OMIM:129600')

The predicate takes two or more disease identifiers (`diagnoses`) as well as their names (`labels`),
and it assigns the individuals based on their diagnoses.

Note, the assignment must be unambiguous; any individual labeled with two or more target diagnoses
(e.g. an individual diagnosed with both *Marfan syndrome* and *Ectopia lentis, familial* in the example above)
will be *omitted* from the analysis.
