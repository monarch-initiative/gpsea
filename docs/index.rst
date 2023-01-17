=============
genophenocorr
=============

Todo intro


Installation
^^^^^^^^^^^^
See :ref:`installation` for instructions on how to set up the package. 

PubMed abstracts
^^^^^^^^^^^^^^^^

To obtain PubMed abstracts, following the instructions
of the NCBI `Download PubMed Data <https://pubmed.ncbi.nlm.nih.gov/download/>`_ website. 

Replacement of biomedical concepts
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We use the `marea <https://github.com/TheJacksonLaboratory/marea>`_ package to replace 
single- or multi-word biomedical concepts with concept identifiers.

See :ref:`rst_marea` for instructions in how to run marea to perform biomedical concept replacement.


WordNet-based synonym Replacement
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The XYZ script from the wn2vec package is used to replace non-biomedical synonyms with WordNet identifiers.

See :ref:`wnreplacement` for details.


word2vec
^^^^^^^^

We perform word2vec embedding using TensorFlow2 with the XYZ script from the wn2vec package. In order to 
compare the results of embedding with and without WordNet replacment, embedding is performed on both datasets.

See :ref:`w2v` for details.


Evaluating concept sets 
^^^^^^^^^^^^^^^^^^^^^^^

Our hypothesis is that non-biomedical concept replacement will improve embeddings as judged by a smaller distance of related 
concepts to each other. TO assess this, we defined XYZ concept sets representing genetic and genomic functions (gene sets) and biomedical concepts taken from MeSH.
The XYZ script from wn2vec is used to assess the mean intracluster cosine distances of the corresponding concepts.

See :ref:`conceptseteval` for details.

--------
Feedback
--------

The best place to leave feedback, ask questions, and report bugs is the `WN2vec Issue Tracker <https://github.com/TheJacksonLaboratory/wn2vec/issues>`_.

.. toctree::
    :caption: Installation & Tutorial
    :name: tutorial
    :maxdepth: 1
    :hidden:

    install
    marea
    wordnetreplacement
    word2vec
    conceptset_evaluation

.. toctree::
    :caption: Project Info
    :name: project-info
    :maxdepth: 1
    :hidden:

    contributing
    authors
    history
    LICENSE
  #  release_howto