.. _genotype-phenotype-classes:


======================================
Compare genotype and phenotype classes
======================================

.. doctest::
  :hide:

  >>> from gpsea import _overwrite

In this section, we show how to test the association between genotype and phenotype classes.
We assume a cohort was preprocessed following the :ref:`input-data` section,
and we use classifiers described in the :ref:`partitioning` to assign each cohort member
into a class along the genotype and phenotype axes.
We use Fisher exact test (FET) to test for differences between the classes
and we apply multiple testing correction to mitigate finding significant associations by chance. 


.. _fisher-exact-test:

*****************
Fisher exact test
*****************

The Fisher exact test calculates the exact probability value
for the relationship between two dichotomous variables.
In our implementation, the two dichotomous variables are the genotype and the phenotype.
For instance, the individuals of the cohort may be divided
according to whether or not they have a nonsense variant
and according to whether or not they have Strabismus (`HP:0000486 <https://hpo.jax.org/browse/term/HP:0000486>`_).


The results of FET are expressed in terms of an exact probability (P-value), varying within 0 and 1.
Two groups are considered statistically significant if the P-value is less
than the chosen significance level (usually :math:`\alpha = 0.05`).

The following graphic shows an example contingency table that is used to conduct a Fisher exact test.
We are comparing the frequency of *strabismus* in individuals with missense and nonsense variants:

.. image:: /img/fisher.png
   :alt: Fisher exact text example
   :align: center
   :width: 600px

To perform the corresponding test in Python, we would use the following code.

>>> import scipy.stats as stats
>>> contingency_table = [
... #   Missense,    Nonsense
...    [17,          6       ],  # Strabismus observed
...    [1,          19       ],  # Strabismus excluded
... ]
>>> oddsratio, p_value = stats.fisher_exact(contingency_table)
>>> float(oddsratio)
53.833333333333336
>>> float(p_value)
5.432292015291845e-06

The ``p_value`` evaluates to `5.432292015291845e-06`, meaning there is a significant difference between the groups.

The Fisher exact test evaluates whether the observed frequencies in a contingency table significantly
deviate from the frequencies we would expect if there were no association between the variables.
We want to test whether the frequency of `HP:0000486`` is significantly higher or lower in
one genotype class compared to what would be expected if there were no association.
Note that by default, the *two-tailed* Fisher exact test is performed, meaning we have no
hypothesis as to whether there is a higher or lower frequency in one of the genotype classes.

However, we are typically interested in testing the associations between the genotype and multiple phenotypic features at once.
GPSEA takes advantage of the HPO structure and simplifies the testing for all HPO terms encoded in the cohort.


****************
Example analysis
****************

Let's illustrate this in a real-life example of the analysis of the association between frameshift variants in *TBX5* gene
and congenital heart defects in the dataset of 156 individuals with mutations in *TBX5* whose signs and symptoms were
encoded into HPO terms, stored as phenopackets of the `GA4GH Phenopacket Schema <https://pubmed.ncbi.nlm.nih.gov/35705716>`_,
and deposited in `Phenopacket Store <https://github.com/monarch-initiative/phenopacket-store>`_.

.. note::

   The shorter version of the same analysis has been presented in the :ref:`tutorial`.


Load cohort
===========

For the purpose of this analysis, we will load the :class:`~gpsea.model.Cohort`
from a `JSON file <https://github.com/P2GX/gpsea/tree/main/docs/cohort-data/TBX5.0.1.20.json>`_.
The cohort was prepared from phenopackets as described in :ref:`create-a-cohort` section,
and then serialized as a JSON file following the instructions in :ref:`cohort-persistence` section.

.. 
   Prepare the JSON file by running the tests in `tests/tests/test_generate_doc_cohorts.py`.

>>> import json
>>> from gpsea.io import GpseaJSONDecoder
>>> fpath_cohort_json = 'docs/cohort-data/TBX5.0.1.20.json'
>>> with open(fpath_cohort_json) as fh:
...     cohort = json.load(fh, cls=GpseaJSONDecoder)
>>> len(cohort)
156


Configure analysis
==================

We want to test the association between frameshift *TBX5* variants and phenotypic abnormalities.
GPSEA exposes a flexible classifier API that lets us create genotype and phenotype classifiers
to assign the cohort members into genotype and phenotype categories based on the variants
and the HPO terms.
We need to create one genotype classifier and one or more phenotype classifiers.


Genotype classifier
-------------------

We want to separate the patients into two classes: a class *with* a frameshift variant
and a class *without* a frameshift variant (i.e. any other heterozygous variant).
We will use the *MANE* transcript for the analysis:

>>> tx_id = 'NM_181486.4'

Building a genotype classifier is a two step process. 
First, we create a :class:`~gpsea.analysis.predicate.VariantPredicate`
to test if the variant is predicted to lead to a frameshift in `NM_181486.4`:

>>> from gpsea.model import VariantEffect
>>> from gpsea.analysis.predicate import variant_effect
>>> is_frameshift = variant_effect(VariantEffect.FRAMESHIFT_VARIANT, tx_id)
>>> is_frameshift.description
'FRAMESHIFT_VARIANT on NM_181486.4'

.. note::

   The :mod:`gpsea.analysis.predicate` documentation lists all available variant predicates
   and :ref:`variant-predicates` exemplifies their usage.

To build a genotype classifier, we wrap `is_frameshift`
in a Monoallelic classifier (:class:`~gpsea.analysis.clf.monoallelic_classifier`),
to classify each *TBX5* cohort member either as an individual with one *frameshift* allele (`Frameshift`)
or as an individual with one *non-frameshift* allele (`Other`):

>>> from gpsea.analysis.clf import monoallelic_classifier
>>> gt_clf = monoallelic_classifier(
...     a_predicate=is_frameshift,
...     a_label="Frameshift",
...     b_label="Other",
... )
>>> gt_clf.class_labels
('Frameshift', 'Other')

.. note::

   See the :ref:`genotype-classifiers` for other genotype classifier examples.

In the subsequent analysis, `gt_clf` assigns an individual into a genotype class.
Note, any individual with :math:`0` or :math:`\ge 2` alleles will be *omitted* from the analysis.


Phenotype classifiers
---------------------

We recommend testing the genotype phenotype association for all HPO terms that annotate the cohort members,
while taking advantage of the HPO graph structure and of the :ref:`true-path-rule`.
We will use the :func:`~gpsea.analysis.clf.prepare_classifiers_for_terms_of_interest`
utility function to generate phenotype classifiers for all HPO terms.

The function needs HPO to prepare classifiers, hence we need to load HPO:

>>> import hpotk
>>> store = hpotk.configure_ontology_store()
>>> hpo = store.load_minimal_hpo(release='v2024-07-01')


and then we can create the classifiers

>>> from gpsea.analysis.clf import prepare_classifiers_for_terms_of_interest
>>> pheno_clfs = prepare_classifiers_for_terms_of_interest(
...     cohort=cohort,
...     hpo=hpo,
... )
>>> len(pheno_clfs)
369

The function finds 369 HPO terms that annotate at least one individual,
including the *indirect* annotations whose presence is implied by the :ref:`true-path-rule`.


.. _phenotype-classes-statistical-analysis:


Statistical analysis
--------------------

We will use :ref:`fisher-exact-test` to test the association
between genotype and phenotype classes, as described previously.

In the case of this cohort, we can test association between having a frameshift variant and one of 369 HPO terms.
However, testing multiple hypotheses on the same dataset increases the risk of finding
a significant association solely by chance.
GPSEA uses a two-pronged strategy to reduce the number of tests and, therefore, mitigate this risk:
a phenotype multiple testing (MT) filter and multiple testing correction (MTC).

Phenotype MT filter selects a (sub)set of HPO terms for testing,
for instance only the user-selected terms (see :class:`~gpsea.analysis.mtc_filter.SpecifiedTermsMtcFilter`)
or the terms selected by :class:`~gpsea.analysis.mtc_filter.IfHpoFilter`.

MTC then adjusts the nominal p values for the increased risk
of false positive G/P associations.
The available MTC procedures are listed in the :ref:`mtc-correction-procedures` section.

We must choose a phenotype MT filter as well as a MTC procedure to perform genotype-phenotype analysis.


.. _default-hpo-analysis:

Default analysis
^^^^^^^^^^^^^^^^

We recommend using Independent filtering for HPO (:class:`~gpsea.analysis.mtc_filter.IfHpoFilter`)
and Benjamini-Hochberg MT correction.
The default analysis can be configured with :func:`~gpsea.analysis.pcats.configure_hpo_term_analysis` convenience method.

>>> from gpsea.analysis.pcats import configure_hpo_term_analysis
>>> analysis = configure_hpo_term_analysis(hpo)

At this point, the ``analysis`` configured to test
a cohort for G/P associations.


.. _custom-hpo-analysis:

Custom analysis
^^^^^^^^^^^^^^^

If the default selection of phenotype MT filter and multiple testing correction is not an option,
we can configure the analysis manually.

First, we choose a phenotype MT filter (e.g. :class:`~gpsea.analysis.mtc_filter.IfHpoFilter`):

>>> from gpsea.analysis.mtc_filter import IfHpoFilter
>>> mtc_filter = IfHpoFilter.default_filter(hpo)

.. note::

   See the :ref:`mtc-filters` section for info regarding other phenotype MT filters.

then a statistical test (e.g. Fisher Exact test):

>>> from gpsea.analysis.pcats.stats import FisherExactTest
>>> count_statistic = FisherExactTest()

.. note::

   See the :mod:`gpsea.analysis.pcats.stats` module for the available multiple testing procedures
   (TL;DR, just Fisher Exact test at this time).

and we finalize the setup by choosing a MTC procedure
(e.g. `fdr_bh` for Benjamini-Hochberg) along with the MTC alpha:

>>> mtc_correction = 'fdr_bh'
>>> mtc_alpha = 0.05

.. note::

   See the :ref:`mtc-correction-procedures` section for a list of available MTC procedure codes.

The final :class:`~gpsea.analysis.pcats.HpoTermAnalysis` is created as:

>>> from gpsea.analysis.pcats import HpoTermAnalysis
>>> analysis = HpoTermAnalysis(
...     count_statistic=count_statistic,
...     mtc_filter=mtc_filter,
...     mtc_correction='fdr_bh',
...     mtc_alpha=0.05,
... )

The ``analysis`` is identical to the one configured in the :ref:`default-hpo-analysis` section.


Analysis
========

We can now test associations between the genotype classes and the HPO terms:

>>> result = analysis.compare_genotype_vs_phenotypes(
...     cohort=cohort,
...     gt_clf=gt_clf,
...     pheno_clfs=pheno_clfs,
... )
>>> len(result.phenotypes)
369
>>> result.total_tests
32


We tested the ``cohort`` for association between the genotype classes (``gt_clf``)
and HPO terms (``pheno_clfs``).
Thanks to phenotype MT filter, we only tested 32 out of 369 terms.
The MT filter report shows the filtering details:

>>> from gpsea.view import MtcStatsViewer
>>> mtc_viewer = MtcStatsViewer()
>>> mtc_report = mtc_viewer.process(result)
>>> mtc_report  # doctest: +SKIP

.. raw:: html
  :file: report/tbx5_frameshift.mtc_report.html

.. doctest:: phenotype-classes
   :hide:

   >>> if _overwrite: mtc_report.write('docs/user-guide/analyses/report/tbx5_frameshift.mtc_report.html')


Genotype phenotype associations
===============================

Last, let's explore the associations. 

GPSEA displays the associations between genotypes and HPO terms in a table,
one HPO term per row. The rows are ordered by the corrected p value and nominal p value in descending order.

>>> from gpsea.view import summarize_hpo_analysis
>>> summary_df = summarize_hpo_analysis(hpo, result)
>>> summary_df  # doctest: +SKIP

.. csv-table:: *TBX5* frameshift vs rest
   :file: report/tbx5_frameshift.csv
   :header-rows: 1

.. doctest:: phenotype-classes
   :hide:

   >>> if _overwrite: summary_df.to_csv('docs/user-guide/analyses/report/tbx5_frameshift.csv')


The table shows that several HPO terms are significantly associated
with presence of a heterozygous (`Frameshift`) frameshift variant in *TBX5*.
For example, `Ventricular septal defect <https://hpo.jax.org/browse/term/HP:0001629>`_
was observed in 42/71 (59%) patients with no frameshift allele (`Other`)
but it was observed in 19/19 (100%) patients with a frameshift allele (`Frameshift`).
Fisher exact test computed a p value of `~0.000242`
and the p value corrected by Benjamini-Hochberg procedure
is `~0.00774`.

