#####
GPSEA
#####


The concept of phenotype denotes the observable attributes of an individual, but in 
medical contexts, the word "phenotype" is used to refer to some deviation from normal morphology, physiology, or behavior 
(c.f. `Deep phenotyping for precision medicine <https://pubmed.ncbi.nlm.nih.gov/22504886/>`_).
A key question in biology and human genetics concerns the relationships between phenotypic abnormalities and genotype.
In Mendelian genetics, the focus is generally placed on the study of whether specific disease-causing alleles
are associated with specific phenotypic manifestations of the disease. 

GPSEA (Genotypes and Phenotypes - Statistical Evaluation of Associations)
is a Python package for finding genotype-phenotype associations.
The input to GPSEA is a collection of `Global Alliance for Genomics and Health (GA4GH) Phenopackets <https://pubmed.ncbi.nlm.nih.gov/35705716/>`_.
GPSEA ingests the phenopackets and analyzes the genotype-phenotype associations.
The genotype can include variant types (e.g., missense vs. premature termination codon),
or variant location in protein motifs or other features.
Phenotype can be represented by `Human Phenotype Ontology (HPO) <https://hpo.jax.org/app/>`_ terms, but using other phenotypes is possible.
Statistical analysis is performed using e.g `Fisher Exact Test <https://en.wikipedia.org/wiki/Fisher%27s_exact_test>`_.

GPSEA integrates with modern interactive computing environments, such as Jupyter notebook.
Therefore, we recommend to install GPSEA into a Jupyter kernel and to perform the analyses
in a Jupyter notebooks.
The documentation is structured such that each section can be executed in a notebook,
and we encourage running the examples while reading the documentation.

The documentation includes :ref:`setup` instructions,
a :ref:`tutorial` with an end-to-end example, and a comprehensive :ref:`user-guide`.
The technical information is available in the API reference.

**********
Literature
**********

We provide recommended reading for background on the study of genotype-phenotype correlations.

- `Orgogozo V, et al. (2015) <https://pubmed.ncbi.nlm.nih.gov/26042146/>`_ The differential view of genotype-phenotype relationships.

********
Feedback
********

The best place to leave feedback, ask questions, and report bugs is the GPSEA `Issue Tracker <https://github.com/P2GX/gpsea/issues>`_.

.. toctree::
    :caption: Table of Contents
    :name: index-toc
    :maxdepth: 1
    :hidden:

    setup
    tutorial
    user-guide/index
    apidocs/modules
