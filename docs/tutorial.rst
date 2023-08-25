.. _tutorial:

========
Tutorial
========


Set up analysis
^^^^^^^^^^^^^^^

Genophenocorr needs HPO to do the analysis. Let's load the ontology:

.. doctest:: tutorial

  >>> import hpotk
  >>> hpo = hpotk.ontology.load.obographs.load_minimal_ontology('data/hp.toy.json')

.. tip::

  Use the latest HPO which you can get at `http://purl.obolibrary.org/obo/hp.json`

TODO - move the code from `workflow` and the notebook here.

Prepare samples
^^^^^^^^^^^^^^^

Now we need some samples. To keep things simple in this tutorial, we will use a toy cohort that is shipped
with the package:

.. doctest:: tutorial

  >>> from genophenocorr.data import get_toy_cohort
  >>> cohort = get_toy_cohort()

.. seealso::

  See :ref:`input-data` section to learn about preparing your data for the analysis.


TODO - move the code from `workflow` and the notebook here.

