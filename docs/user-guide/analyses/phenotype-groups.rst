.. _genotype-phenotype-groups:


=====================================
Compare genotype and phenotype groups
=====================================


.. _fisher-exact-test:

***********************
Fisher exact test (FET)
***********************

The Fisher exact test (FET) calculates the exact probability value
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
one genotype group compared to what would be expected if there were no association.
Note that by default, the *two-tailed* Fisher exact test is performed, meaning we have no
hypothesis as to whether there is a higher or lower frequency in one of the genotype groups.

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
from a `JSON file <https://github.com/monarch-initiative/gpsea/tree/main/docs/cohort-data/TBX5.0.1.20.json>`_.
The cohort was prepared from phenopackets as described in :ref:`create-cohort-from-phenopackets` section,
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
GPSEA exposes a flexible predicate API that lets us create genotype and phenotype predicates
to assign the cohort members into genotype and phenotype categories based on the variants
and the HPO terms. We need to create one genotype predicate and one or more phenotype predicates.


Genotype predicate
------------------

We want to separate the patients into two groups: a group *with* a frameshift variant
and a group *without* a frameshift variant (i.e. any other heterozygous variant).
We will use the *MANE* transcript for the analysis:

>>> tx_id = 'NM_181486.4'

Building a genotype predicate is a two step process. 
First, we create a :class:`~gpsea.analysis.predicate.genotype.VariantPredicate`
to test if the variant is predicted to lead to a frameshift in `NM_181486.4`:

>>> from gpsea.model import VariantEffect
>>> from gpsea.analysis.predicate.genotype import VariantPredicates
>>> is_frameshift = VariantPredicates.variant_effect(VariantEffect.FRAMESHIFT_VARIANT, tx_id)
>>> is_frameshift.get_question()
'FRAMESHIFT_VARIANT on NM_181486.4'

and then we wrap `is_frameshift` in a :class:`~gpsea.analysis.predicate.genotype.monoallelic_predicate` 
to classify each *TBX5* cohort member either as an individual with one frameshift allele (`Frameshift`)
or as an idividual with one non-frameshift allele (`Other`):

>>> from gpsea.analysis.predicate.genotype import monoallelic_predicate
>>> gt_predicate = monoallelic_predicate(
...     a_predicate=is_frameshift,
...     a_label="Frameshift",
...     b_label="Other",
... )
>>> gt_predicate.display_question()
'Allele group: Frameshift, Other'

In the subsequent analysis, `gt_predicate` will assign a cohort member into the respective group.
Note, any patient with :math:`0` or :math:`\ge 2` alleles will be *omitted* from the analysis.


Phenotype predicates
--------------------

We recommend testing the genotype phenotype association for all HPO terms that annotate the cohort members,
while taking advantage of the HPO graph structure and of the :ref:`true-path-rule`.
We will use the :func:`~gpsea.analysis.predicate.phenotype.prepare_predicates_for_terms_of_interest`
utility function to generate phenotype predicates for all HPO terms.

The function needs HPO to prepare predicates, hence we need to load HPO:

>>> import hpotk
>>> store = hpotk.configure_ontology_store()
>>> hpo = store.load_minimal_hpo(release='v2024-07-01')


and then we can create the predicates

>>> from gpsea.analysis.predicate.phenotype import prepare_predicates_for_terms_of_interest
>>> pheno_predicates = prepare_predicates_for_terms_of_interest(
...     cohort=cohort,
...     hpo=hpo,
... )
>>> len(pheno_predicates)
369

The function finds 369 HPO terms that annotate at least one individual,
including the *indirect* annotations whose presence is implied by the :ref:`true-path-rule`.


Statistical test
----------------

We will use :ref:`fisher-exact-test` to test the association
between genotype and phenotype groups, as described previously.

>>> from gpsea.analysis.pcats.stats import FisherExactTest
>>> count_statistic = FisherExactTest()

FET will compute a p value for each genotype phenotype group.


Multiple testing correction
---------------------------

In the case of this cohort, we could test association between having a frameshift variant and one of 369 HPO terms.
However, testing multiple hypotheses on the same dataset increases the risk of finding a significant association
by chance.
GPSEA uses a two-pronged strategy to mitigate this risk - a phenotype MTC filter and multiple testing correction.

.. note::

   See the :ref:`mtc` section for more info on multiple testing procedures.

Here we will use a combination of the HPO MTC filter (:class:`~gpsea.analysis.mtc_filter.HpoMtcFilter`)
with Benjamini-Hochberg procedure (``mtc_correction='fdr_bh'``)
with a false discovery control level set to `0.05` (``mtc_alpha=0.05``):

>>> from gpsea.analysis.mtc_filter import HpoMtcFilter
>>> mtc_filter = HpoMtcFilter.default_filter(hpo, term_frequency_threshold=0.2)
>>> mtc_correction = 'fdr_bh'
>>> mtc_alpha = 0.05


Final analysis
--------------

We finalize the analysis setup by putting all components together
into :class:`~gpsea.analysis.pcats.HpoTermAnalysis`:

>>> from gpsea.analysis.pcats import HpoTermAnalysis
>>> analysis = HpoTermAnalysis(
...     count_statistic=count_statistic,
...     mtc_filter=mtc_filter,
...     mtc_correction=mtc_correction,
...     mtc_alpha=mtc_alpha,
... )


Analysis
========

We can now execute the analysis:

>>> result = analysis.compare_genotype_vs_phenotypes(
...     cohort=cohort,
...     gt_predicate=gt_predicate,
...     pheno_predicates=pheno_predicates,
... )
>>> len(result.phenotypes)
369
>>> result.total_tests
24


Thanks to phenotype MTC filter, we only tested 24 out of 369 terms.
We can learn more by showing the MTC filter report:

>>> from gpsea.view import MtcStatsViewer
>>> mtc_viewer = MtcStatsViewer()
>>> mtc_report = mtc_viewer.process(result)
>>> mtc_report  # doctest: +SKIP

.. raw:: html
  :file: report/tbx5_frameshift.mtc_report.html

.. doctest:: phenotype-groups
   :hide:

   >>> mtc_report.write('docs/user-guide/analyses/report/tbx5_frameshift.mtc_report.html')  # doctest: +SKIP


Genotype phenotype associations
===============================

Last, let's explore the associations. The results include a table with all tested HPO terms
ordered by the corrected p value (Benjamini-Hochberg FDR):

>>> from gpsea.view import summarize_hpo_analysis
>>> summary_df = summarize_hpo_analysis(hpo, result)
>>> summary_df  # doctest: +SKIP

.. csv-table:: *TBX5* frameshift vs rest
   :file: report/tbx5_frameshift.csv
   :header-rows: 2

.. doctest:: phenotype-groups
   :hide:

   >>> summary_df.to_csv('docs/user-guide/analyses/report/tbx5_frameshift.csv')  # doctest: +SKIP


The table shows that several HPO terms are significantly associated
with presence of a heterozygous (`Frameshift`) frameshift variant in *TBX5*.
For example, `Ventricular septal defect <https://hpo.jax.org/browse/term/HP:0001629>`_
was observed in 42/71 (59%) patients with no frameshift allele (`Other`)
but it was observed in 19/19 (100%) patients with a frameshift allele (`Frameshift`).
Fisher exact test computed a p value of `~0.000242`
and the p value corrected by Benjamini-Hochberg procedure
is `~0.005806`.

