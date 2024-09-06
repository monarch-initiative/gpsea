.. _hpo_predicate:


Propagating phenotype predicate
===============================

When testing for presence or absence of an HPO term, the propagating phenotype predicate
leverages the :ref:`true-path-rule` to take advantage of the HPO hierarchy.
In result, an individual annotated with a term is implicitly annotated with all its ancestors.
For instance, an individual annotated with `Ectopia lentis <https://hpo.jax.org/browse/term/HP:0001083>`_
is also annotated with `Abnormal lens morphology <https://hpo.jax.org/browse/term/HP:0000517>`_,
`Abnormal anterior eye segment morphology <https://hpo.jax.org/browse/term/HP:0004328>`_,
`Abnormal eye morphology <https://hpo.jax.org/browse/term/HP:0012372>`_, ...

Similarly, all descendants of a term, whose presence was specifically excluded in an individual,
are implicitly excluded.

:class:`~gpsea.analysis.predicate.phenotype.PropagatingPhenotypePredicate` implements this logic.

Example
-------

Here we show how to set up :class:`~gpsea.analysis.predicate.phenotype.PropagatingPhenotypePredicate`
to test for a presence of `Abnormal lens morphology <https://hpo.jax.org/browse/term/HP:0000517>`_.


>>> from gpsea.analysis.predicate.phenotype import HpoPredicate
>>> query = hpotk.TermId.from_curie('HP:0000517')
>>> pheno_predicate = HpoPredicate(
...     hpo=hpo,
...     query=query,
... )
>>> pheno_predicate.display_question()
'Is Abnormal lens morphology present in the patient: Yes, No'



missing_implies_phenotype_excluded
----------------------------------

In many cases, published reports of clinical data about individuals with rare diseases describes phenotypic features that were observed, but do not 
provide a comprehensive list of features that were explicitly excluded. By default, GPSEA will only include features that are recorded as observed or excluded in a phenopacket.
Setting this argument to True will cause "n/a" entries to be set to "excluded". We provide this option for exploration but do not recommend its use for the 
final analysis unless the assumption behind it is known to be true.