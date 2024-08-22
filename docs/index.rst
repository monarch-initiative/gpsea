=====
GPSEA
=====


The concept of phenotype denote the observable attributes of an individual, but in 
medical contexts, the word "phenotype" is used to refer to some deviation from normal morphology, physiology, or behavior 
(c.f. `Deep phenotyping for precision medicine <https://pubmed.ncbi.nlm.nih.gov/22504886/>`_).
A key question in biology and human genetics concerns the relationships between phenotypic abnormalities and genotype. In Mendelian
genetics, the focus is generally placed on the study of whether specific disease-causing alleles are associated with specific phenotypic 
manifestations of the disease. 

`GPSEA`  (genotypes and phenotypes - study and evaluation of associations) is a Python package designed to support genotype-phenotype correlation analysis.
We pronounce GPSEA as "G"-"P"-"C". The input to `GPSEA` is a collection of `Global Alliance for Genomics and Health (GA4GH) Phenopackets <https://pubmed.ncbi.nlm.nih.gov/35705716/>`_.
`gpsea` ingests data from these phenopackets and performs analysis of the correlation of specific variants,
variant types (e.g., missense vs. premature termination codon), or variant location in protein motifs or other features.
The phenotypic abnormalities are represented by `Human Phenotype Ontology (HPO) <https://hpo.jax.org/app/>`_ terms.
Statistical analysis is performed using a `Fisher Exact Test <https://en.wikipedia.org/wiki/Fisher%27s_exact_test>`_,
and results are reported for each tested HPO term.

We recommend that users create a Jupyter notebook for each cohort of patients they would like to test.

This documentation includes installation instructions, a brief tutorial, and a comprehensive user guide.
The technical information is available in API reference.


Literature
^^^^^^^^^^

We provide recommended reading for background on the study of genotype-phenotype correlations.

- `Orgogozo V, et al. (2015) <https://pubmed.ncbi.nlm.nih.gov/26042146/>`_ The differential view of genotype-phenotype relationships.

--------
Feedback
--------

The best place to leave feedback, ask questions, and report bugs is the `GPSEA Issue Tracker <https://github.com/monarch-initiative/genophenocorr/issues>`_.

.. toctree::
    :caption: Installation & Tutorial
    :name: tutorial
    :maxdepth: 1
    :hidden:

    installation
    tutorial
    user-guide/index
    apidocs/modules

..  Comment out the below for now

..  .. toctree::
    :caption: Project Info
    :name: project-info
    :maxdepth: 1
    :hidden:

..    contributing
    authors
    history
    LICENSE
