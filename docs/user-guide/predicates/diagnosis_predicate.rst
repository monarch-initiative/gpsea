.. _diagnosis-predicate:

========================
Partition by a diagnosis
========================

It is also possible to bin the individuals based on a diagnosis.
The :func:`~gpsea.analysis.predicate.genotype.diagnosis_predicate` 
prepares a genotype predicate for assigning an individual into a diagnosis group:

>>> from gpsea.analysis.predicate.genotype import diagnosis_predicate
>>> gt_predicate = diagnosis_predicate(
...     diagnoses=('OMIM:154700', 'OMIM:129600'),
...     labels=('Marfan syndrome', 'Ectopia lentis, familial'),
... )
>>> gt_predicate.display_question()
'What disease was diagnosed: OMIM:154700, OMIM:129600'

Note, an individual must match only one diagnosis group. Any individuals labeled with two or more diagnoses
(e.g. an individual with both *Marfan syndrome* and *Ectopia lentis, familial*)
will be automatically omitted from the analysis.