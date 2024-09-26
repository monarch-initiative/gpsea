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
and deposited in `Phenopacket Store <https://github.com/monarch-initiative/phenopacket-store>`_ (version `0.1.18`).

.. note::

   The shorter version of the same analysis has been presented in the :ref:`tutorial`.


Create cohort
=============

We will load and transform the phenopackets into a :class:`~gpsea.model.Cohort`,
as described in :ref:`input-data` section. Briefly, we will load the phenopackets:

>>> from ppktstore.registry import configure_phenopacket_registry
>>> registry = configure_phenopacket_registry()
>>> with registry.open_phenopacket_store(release='0.1.18') as ps:
...     phenopackets = tuple(ps.iter_cohort_phenopackets('TBX5'))
>>> len(phenopackets)
156

followed by loading HPO release `v2024-07-01`:

>>> import hpotk
>>> store = hpotk.configure_ontology_store()
>>> hpo = store.load_minimal_hpo(release='v2024-07-01')

and we will perform Q/C and functional annotations for the mutations
with the default cohort creator:

>>> from gpsea.preprocessing import configure_caching_cohort_creator, load_phenopackets
>>> cohort_creator = configure_caching_cohort_creator(hpo)
>>> cohort, qc_results = load_phenopackets(phenopackets, cohort_creator)  # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
Individuals Processed: ...
>>> qc_results.summarize()  # doctest: +SKIP
Validated under none policy
No errors or warnings were found


Configure analysis
==================

We want to test the association between frameshift *TBX5* variants and phenotypic abnormalities.
GPSEA exposes a flexible predicate API that lets us create genotype and phenotype predicates
to assign the cohort members into genotype and phenotype categories based on the variants
and the HPO terms. We need to create one genotype predicate and one or more phenotype predicates.


Genotype predicate
------------------

We want to separate the patients into two groups: a group *with* a frameshift variant
and a group *without* a frameshift variant, based on the functional annotation.
We will use the *MANE* transcript for the analysis:

Building a genotype predicate is a two step process. 
First, we create a :class:`~gpsea.analysis.predicate.genotype.VariantPredicate`
to test if the variant leads to a frameshift (in this case):

>>> from gpsea.model import VariantEffect
>>> from gpsea.analysis.predicate.genotype import VariantPredicates
>>> tx_id = 'NM_181486.4'
>>> is_frameshift = VariantPredicates.variant_effect(VariantEffect.FRAMESHIFT_VARIANT, tx_id)
>>> is_frameshift.get_question()
'FRAMESHIFT_VARIANT on NM_181486.4'

and then we choose the expected mode of inheritance to test. In case of *TBX5*,
we expect the autosomal dominant mode of inheritance:

>>> from gpsea.analysis.predicate.genotype import autosomal_dominant
>>> gt_predicate = autosomal_dominant(is_frameshift)
>>> gt_predicate.display_question()
'What is the genotype group: HOM_REF, HET'

`gt_predicate` will assign the patients with no frameshift variant allele into `HOM_REF` group
and the patients with one frameshift allele will be assigned into `HET` group.
Note, any patient with 2 or more alleles will be *omitted* from the analysis.

.. note::

   Mode of inheritance testing is not the only way to dissect by a genotype.
   See the :ref:`genotype-predicates` section for more info.


Phenotype predicates
--------------------

We recommend testing the genotype phenotype association for all HPO terms that are present in 2 or more cohort members,
while taking advantage of the HPO graph structure and of the :ref:`true-path-rule`.
We will use the :func:`~gpsea.analysis.predicate.phenotype.prepare_predicates_for_terms_of_interest`
utility function to generate phenotype predicates for all HPO terms:

>>> from gpsea.analysis.predicate.phenotype import prepare_predicates_for_terms_of_interest
>>> pheno_predicates = prepare_predicates_for_terms_of_interest(
...     cohort=cohort,
...     hpo=hpo,
...     min_n_of_patients_with_term=2,
... )
>>> len(pheno_predicates)
260

The function finds all HPO terms that annotate at least *n* (``min_n_of_patients_with_term=2`` above) individuals,
including the *indirect* annotations whose presence is implied by the true path rule.


Statistical test
----------------

We will use :ref:<fisher-exact-test> to test the association
between genotype and phenotype groups, as described previously.

>>> from gpsea.analysis.pcats.stats import FisherExactTest
>>> count_statistic = FisherExactTest()

FET will compute a p value for each genotype phenotype group.


Multiple testing correction
---------------------------

In the case of this cohort, we could test association between having a frameshift variant and one of 260 HPO terms.
However, testing multiple hypotheses on the same dataset increases the risk of finding a significant association
by chance.
GPSEA uses a two-pronged strategy to mitigate this risk - use Phenotype MTC filter and multiple testing correction.

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
260
>>> result.total_tests
16


Thanks to Phenotype MTC filter, we only tested 16 out of 260 terms.
We can learn more by showing the MTC filter report:

>>> from gpsea.view import MtcStatsViewer
>>> mtc_viewer = MtcStatsViewer()
>>> mtc_report = mtc_viewer.process(result)
>>> with open('docs/user-guide/report/tbx5_frameshift.mtc_report.html', 'w') as fh:  # doctest: +SKIP
...     _ = fh.write(mtc_report)


.. raw:: html
  :file: report/tbx5_frameshift.mtc_report.html


Genotype phenotype associations
===============================

Last, let's explore the associations. The results include a table with all tested HPO terms
ordered by the corrected p value (Benjamini-Hochberg FDR):

>>> from gpsea.view import summarize_hpo_analysis
>>> summary_df = summarize_hpo_analysis(hpo, result)
>>> summary_df.to_csv('docs/user-guide/report/tbx5_frameshift.csv')  # doctest: +SKIP

.. csv-table:: *TBX5* frameshift vs rest
   :file: report/tbx5_frameshift.csv
   :header-rows: 2


The table shows that several HPO terms are significantly associated
with presence of a heterozygous (`HET`) frameshift variant in *TBX5*.
For example, `Ventricular septal defect <https://hpo.jax.org/browse/term/HP:0001629>`_
was observed in 31/60 (52%) patients with a missense variant
but it was observed in 19/19 (100%) patients with a frameshift variant.
Fisher exact test computed a p value of `~0.000242`
and the p value corrected by Benjamini-Hochberg procedure
is `~0.00387`.

The table includes all HPO terms of the cohort, including the terms that were not selected for testing
and thus have no associated p value.

