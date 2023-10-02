.. _tutorial:

========
Tutorial
========

The tutorial demonstrates how to load an example Phenopacket cohort and perform genotype-phenotype analysis.

Set up analysis
^^^^^^^^^^^^^^^

`genophenocorr` needs HPO to do the analysis. Let's load the ontology:

.. doctest:: tutorial

  >>> import hpotk
  >>> hpo = hpotk.load_minimal_ontology('data/hp.toy.json')

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

We can then view the data using the list commands. 

.. doctest:: tutorial
  
  >>> sorted(cohort.list_all_patients())
  ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']
  >>> sorted(cohort.list_all_phenotypes())
  [('HP:0001166', 14), ('HP:0001250', 20), ('HP:0001257', 17)]
  >>> sorted(cohort.list_all_variants())
  [('1_281_A/G', 16), ('1_361_TTC/T', 13)]
  >>> sorted(cohort.list_all_proteins())
  [('NP_09876.5', 26)]
  >>> tx_dict = cohort.list_data_by_tx('NM_1234.5')
  >>> sorted(tx_dict['NM_1234.5'].items())
  [('frameshift_variant', 1), ('missense_variant', 1)]

Using the counts, we can choose and run what analyses we want.
For instance, we can partition the patients into two groups based on presence/absence of a *frameshift* variant:

.. doctest:: tutorial

  >>> from genophenocorr.analysis import CohortAnalysis
  >>> from genophenocorr.constants import VariantEffect
  >>> cohort_analysis = CohortAnalysis(cohort, 'NM_1234.5', hpo, include_unmeasured=False)
  >>> frameshift = cohort_analysis.compare_by_variant_type(VariantEffect.FRAMESHIFT_VARIANT)
  >>> frameshift # doctest: +NORMALIZE_WHITESPACE
                              With frameshift_variant         Without frameshift_variant
                                                Count Percent                      Count Percent  p-value Corrected p-values
  HP:0001166 (Arachnodactyly)                       4  30.77%                         10  76.92%  0.04718            0.14154
  HP:0001250 (Seizure)                             11  84.62%                          9  69.23%  0.64472            1.00000
  HP:0001257 (Spasticity)                           8  61.54%                          9  69.23%  1.00000            1.00000


Or perform similar partitioning based on presence/absence of a *missense* variant:

.. doctest:: tutorial

  >>> missense = cohort_analysis.compare_by_variant_type(VariantEffect.MISSENSE_VARIANT)
  >>> missense # doctest: +NORMALIZE_WHITESPACE
                              With missense_variant         Without missense_variant
                                              Count Percent                    Count Percent   p-value Corrected p-values
  HP:0001166 (Arachnodactyly)                    13  81.25%                        1  10.00%  0.000781           0.002342
  HP:0001257 (Spasticity)                        11  68.75%                        6  60.00%  0.692449           1.000000
  HP:0001250 (Seizure)                           12  75.00%                        8  80.00%  1.000000           1.000000


The tables present the HPO terms that annotate the cohort members and report their counts and p values
for each genotype group. The rows are sorted by the p value in ascending order.
