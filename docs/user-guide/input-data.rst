.. _input-data:

==========
Input data
==========

The `gpsea` analysis needs to be provided with a standardized form of genotype and phenotype data.
The analyses require an instance of :class:`~gpsea.model.Cohort` that consists
of individuals in form of a :class:`~gpsea.model.Patient` class.

.. seealso::

  See the :ref:`cohort-exploratory` section for more info.

The first step of the `gpsea` analysis involves standardization of the genotype and phenotype data
and performing functional annotation of the variants. Here we describe how to prepare a `Cohort`
for the exploratory and downstream analysis.


***************************************
Create a cohort from GA4GH phenopackets
***************************************

The easiest way to input data into `gpsea` is to use the
`GA4GH Phenopacket Schema <https://phenopacket-schema.readthedocs.io/en/latest>`_ phenopackets.
`gpsea` provides an out-of-the-box solution for loading a cohort from a folder of phenopacket JSON files.


Create cohort creator
=====================

Next, let's prepare a :class:`~gpsea.preprocessing.CohortCreator` that will turn a phenopacket collection
into a :class:`~gpsea.model.Cohort`. The cohort creator also performs an input validation.
The validation needs Human Phenotype Ontology data.
Let's start with loading Human Phenotype Ontology, a requisite for the input Q/C steps. We'll use the amazing
`hpo-toolkit <https://github.com/TheJacksonLaboratory/hpo-toolkit>`_ library which is installed along with
the standard `gpsea` installation:

>>> import hpotk
>>> store = hpotk.configure_ontology_store()
>>> hpo = store.load_minimal_hpo(release='v2024-07-01')

The easiest way to get the `CohortCreator` is to use the
:func:`~gpsea.preprocessing.configure_caching_cohort_creator` convenience method:

.. doctest:: input-data 

  >>> from gpsea.preprocessing import configure_caching_cohort_creator

  >>> cohort_creator = configure_caching_cohort_creator(hpo) 

.. note::

  The default `CohortCreator` will call Variant Effect Predictor and Uniprot APIs
  to perform the functional annotation and protein annotation, 
  and the responses will be cached in the current working directory to reduce the network bandwidth.
  See the :func:`~gpsea.preprocessing.configure_caching_cohort_creator` pydoc for more options.


Load phenopackets
=================

We can create a cohort starting from a collection of `Phenopacket` objects
provided by Python  `Phenopackets <https://pypi.org/project/phenopackets>`_ library.
For the purpose of this example, we will load a cohort of patients with pathogenic mutations in *RERE* gene
included in the release `0.1.18` of `Phenopacket Store <https://github.com/monarch-initiative/phenopacket-store>`_.
We use `Phenopacket Store Toolkit <https://github.com/monarch-initiative/phenopacket-store-toolkit>`_
(``ppktstore`` in the code) to reduce the boilerplate code needed to load the phenopackets:

>>> from ppktstore.registry import configure_phenopacket_registry
>>> registry = configure_phenopacket_registry()
>>> with registry.open_phenopacket_store(release='0.1.18') as ps:
...     phenopackets = tuple(ps.iter_cohort_phenopackets('RERE'))
>>> len(phenopackets)
19

We loaded 19 phenopackets. Now we can turn the phenopackets into a `Cohort`
using the `cohort_creator` and the :func:`~gpsea.preprocessing.load_phenopackets`
loader function:

>>> from gpsea.preprocessing import load_phenopackets
>>> cohort, qc_results = load_phenopackets(phenopackets, cohort_creator)  # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
Individuals Processed: ...
>>> len(cohort)
19

The cohort includes all 19 individuals. 
On top of the ``cohort``, the loader function also provides Q/C results  ``qc_results``.
We call :meth:`~gpsea.preprocessing.PreprocessingValidationResult.summarize`
to display the Q/C summary:

>>> qc_results.summarize()  # doctest: +SKIP
Validated under none policy
No errors or warnings were found


Alternative phenopacket sources
===============================

In case you do not already have a `Phenopacket` collection at your fingertips,
GPSEA provides a few other convenience functions for loading phenopackets from JSON files.

The :func:`~gpsea.preprocessing.load_phenopacket_files` function can be used to load
a bunch of phenopacket JSON files:

>>> from gpsea.preprocessing import load_phenopacket_files
>>> pp_files = ('path/to/phenopacket1.json', 'path/to/phenopacket2.json')
>>> cohort, qc_results = load_phenopacket_files(pp_files, cohort_creator)  # doctest: +SKIP

or you can load an entire directory of JSON files with :func:`~gpsea.preprocessing.load_phenopacket_folder`:

>>> from gpsea.preprocessing import load_phenopacket_folder
>>> pp_dir = 'path/to/folder/with/many/phenopacket/json/files'
>>> cohort, qc_results = load_phenopacket_folder(pp_dir, cohort_creator)  # doctest: +SKIP

