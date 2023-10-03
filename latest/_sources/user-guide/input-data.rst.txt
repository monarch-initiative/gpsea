.. _input-data:

==========
Input data
==========

The `genophenocorr` analysis needs to be provided with a standardized form of genotype and phenotype data.
The analyses require an instance of :class:`genophenocorr.model.Cohort` that consists
of :class:`genophenocorr.model.Patient`\ s - the cohort members. The cohort and the members
hold the standardized data and provide convenience functions for dataset exploration.

.. seealso::

  See the :ref:`cohort-exploratory` section for more info.

The first step of the `genophenocorr` analysis involves standardization of the genotype and phenotype data
and performing functional annotation of the variants. Here we describe how to prepare a `Cohort`
for the exploratory and downstream analysis.

Create a cohort from GA4GH phenopackets
---------------------------------------

The easiest way to input data into `genophenocorr` is to use the
`GA4GH Phenopacket Schema <https://phenopacket-schema.readthedocs.io/en/latest>`_ phenopackets.
`genophenocorr` provides :class:`genophenocorr.preprocessing.PhenopacketPatientCreator`,
an out-of-the-box class for loading phenopackets.


Let's start with loading Human Phenotype Ontology, a requisite for the input Q/C steps. We'll use the amazing
`hpo-toolkit <https://github.com/TheJacksonLaboratory/hpo-toolkit>`_ library which is installed along with
the standard `genophenocorr` installation:

.. doctest:: input-data

  >>> import hpotk

  >>> hpo = hpotk.load_minimal_ontology('http://purl.obolibrary.org/obo/hp.json')

Next, let's create the `PhenopacketPatientCreator`. We use a convenience method
:func:`genophenocorr.preprocessing.configure_caching_patient_creator`:

.. doctest:: input-data

  >>> from genophenocorr.preprocessing import configure_caching_patient_creator

  >>> patient_creator = configure_caching_patient_creator(hpo)

.. note::

  By default, the method creates the patient creator that will call Variant Effect Predictor
  and Uniprot APIs to perform the functional annotation and protein annotation and cache the responses
  in the current working directory to save the bandwidth.
  See the :func:`genophenocorr.preprocessing.configure_caching_patient_creator` for more configuration options.

We can create a cohort starting from a folder with phenopackets stored as JSON files.
For the purpose of this example, we will use a folder `simple_cohort` with 5 example phenopackets located in
`genophenocorr's repository <https://github.com/monarch-initiative/genophenocorr/tree/main/docs/data/simple_cohort>`_

.. doctest:: input-data

  >>> import os
  >>> simple_cohort_path = os.path.join(os.getcwd(), 'data', 'simple_cohort')

Here we walk the file system, load all phenopacket JSON files, and transform the phenopackets into instances of
:class:`genophenocorr.model.Patient`:

.. doctest:: input-data

  >>> import os
  >>> from phenopackets import Phenopacket
  >>> from google.protobuf.json_format import Parse

  >>> patients = []
  >>> for dirpath, _, filenames in os.walk(simple_cohort_path):
  ...   for filename in filenames:
  ...     if filename.endswith('.json'):
  ...       pp_path = os.path.join(dirpath, filename)
  ...       with open(pp_path) as fh:
  ...         pp = Parse(fh.read(), Phenopacket())
  ...       patient = patient_creator.create_patient(pp)
  ...       patients.append(patient)


  >>> f'Loaded {len(patients)} phenopackets'
  'Loaded 5 phenopackets'

Now we can construct a `Cohort`:

.. doctest:: input-data

  >>> from genophenocorr.model import Cohort

  >>> cohort = Cohort.from_patients(patients)
  >>> f'Created a cohort with {len(cohort)} members'
  'Created a cohort with 5 members'


Create a cohort from other data
-------------------------------

TODO - describe how to construct a Patient from raw HPO terms and variant coordinates.

