.. _user-guide:

##########
User guide
##########

The user guide is the primary reference for all GPSEA functionality
and it aims to connect all the information you need to build
your own G/P association analysis workflows.

A typical workflow starts with ingest, Q/C, and exploration of the cohort data.
Then, one or more G/P hypotheses need to be conceived and shaped as an analysis.
The shaping includes selecting the ways for partitioning the cohort along the genotype and phenotype axes,
as well as the statistical test and multiple-testing correction.
The resultant tables and figures can include nominal and corrected p values
or a figure with phenotype score distributions or the survival curves.

The guide is broken into short self-contained sections,
each discussing a single step of a G/P association workflow.
We encourage the users to run the code examples
in an interactive Python environment, such as Jupyter notebook.


.. toctree::
  :maxdepth: 1
  :caption: Contents:

  input-data
  exploratory
  predicates/index
  analyses/index
  mtc
  glossary
