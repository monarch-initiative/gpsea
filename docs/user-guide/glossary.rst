.. _glossary:

========
Glossary
========

The glossary summarizes several frequently used concepts.

.. _true-path-rule:

True path rule
~~~~~~~~~~~~~~

The true path rule of ontologies states that an item (e.g. an individual) annotated with an ontology term
(e.g. `Focal-onset seizure <https://hpo.jax.org/browse/term/HP:0007359>`_)
is implicitly annotated with all its *ancestor* terms
(`Seizure <https://hpo.jax.org/browse/term/HP:0001250>`_,
`Abnormal nervous system physiology <https://hpo.jax.org/browse/term/HP:0012638>`_, ...).
Conversely, exclusion of a term (e.g. `Abnormal ventricular septum morphology <https://hpo.jax.org/browse/term/HP:0010438>`_)
implies exclusion of all its *descendants*
(`Ventricular septal defect <https://hpo.jax.org/browse/term/HP:0001629>`_,
`Ventricular septal aneurysm <https://hpo.jax.org/browse/term/HP:0030957>`_, ...).
