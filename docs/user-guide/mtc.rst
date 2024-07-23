.. _mtc:

===========================
Multiple-testing correction
===========================

Background
~~~~~~~~~~

A p-value is the probability that a test result, under the null hypothesis, assumes the observed or a more extreme value. It is important to realize that if we
perform many tests, we are likely to get a "significant" result by chance alone. For instance, if 
we test a null hypothesis that is true using a significance level of
:math:`\alpha = 0.05`, then there is a probability of :math:`1-\alpha = 0.95` of
arriving at a correct conclusion of non-significance. If we now test
two independent true null hypotheses, the probability that neither
test is significant is :math:`0.95\times 0.95 = 0.90.` If we test 20
independent null hypotheses, the probability that none will be
significant is then :math:`(0.95)^{20}=0.36`. This corresponds to a
probability of :math:`1-0.36=0.64` of getting at least one spurious
significant result, and the expected number of spurious significant
results in 20 tests is :math:`20\times 0.05=1`. If we perform 100 such
tests, the probability that none will be significant is
:math:`(0.95)^{100}=0.01` and there is a 99\% probability of getting at
least one significant result.


Implementation in genophenocorr
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
By default, genephenocorr performs a hypothesis test for each HPO term found at least twice
in the cohort, meaning that we may
perform up to hundreds of tests. Therefore, unless we take into
account the fact that multiple statistical tests are being performed,
it is likely that we will obtain one or more false-positive
results.

The genephenocorr offers two approaches to mitigate this problem: multiple-testing correction procedures
and heuristic methods to reduce the number of terms to be tested. 

Multiple-testing correction (MTC) procedures
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A number of MTC procedures have
been developed to limit the probability of false-positive results. The
MTC procedures differ in complexity, in their assumptions about the
data, and in the type of control they provide.

The genephenocorr package uses the Python package `statsmodels <https://www.statsmodels.org/devel/>`_ to implement
MTC. See the `docuumentation <https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html>`_ for details;
the following table shows allowable options.

+---------------+--------------+
| MTC procedure | abbreviation |
+===============+==============+
| bonferroni    | b            |
+---------------+--------------+
| sidak         | s            |
+---------------+--------------+
|  holm-sidak   |     hs       |
+---------------+--------------+
|     holm      |      h       |
+---------------+--------------+
| simes-hochberg|   sh         |
+---------------+--------------+
|     hommel    |  ho          |
+---------------+--------------+
|     fdr_bh    |              |
+---------------+--------------+
|    fdr_by     |              |
+---------------+--------------+
|     fdr_tsbh  |              |
+---------------+--------------+
|     fdr_tsbky |              |
+---------------+--------------+
|     fdr_gbs   |              |
+---------------+--------------+


The oldest and simplest MTC procedure is the Bonferroni
correction. The Bonferroni procedure thus provides control of the family-wise
error rate (FWER), which is the probability of at least one Type I
error.  The Bonferroni method multiplies the p-value
returned by each test (which is call the *nominal* p-value)
by the number of tests performed (the result is capped at 1.0). This is the **default** method in genephenocorr.
Alternatively, procedures that control the false-discovery rate (FDR),
limit the proportion of significant results that are type I
errors (false discoveries). The Benjamini and Hochberg method (``fdr_bh``) is probably the most commonly used one. 

MTC Heuristics: Choosing which terms to test
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

By choosing which terms to test we can reduce the overall MTC burden. For example, if we choose to test only ten terms out of 450, then
the mutliplication factor of the Bonferroni correction is only 10 instead of 450, and more p-values may "survive" the multiple-testing correction.

The genephenocorr package offers two options, each of which can be combined with any of the MTC methods (in this example, we
choose the Benjamini Hochberg FDR method).

specify terms strategy
~~~~~~~~~~~~~~~~~~~~~~

If you have a specific hypothesis as to which terms may be different between groups, then you can specify these terms using the ``specify_terms_strategy`` function.

For example, let's say that we want to specifically test 
`Abnormal putamen morphology (HP:0031982) <https://hpo.jax.org/browse/term/HP:0031982>`_ and
`Abnormal caudate nucleus morphology (HP:0002339) <https://hpo.jax.org/browse/term/HP:0002339>`_.

We create a set with these two terms and the pass the set as an argument to ``specify_terms_strategy``.



>>> import hpotk
>>> store = hpotk.configure_ontology_store()
>>> hpo = store.load_minimal_hpo(release='v2024-03-06')
>>> from genophenocorr.analysis import configure_cohort_analysis, CohortAnalysisConfiguration
>>> from genophenocorr.model import VariantEffect
>>> from genophenocorr.analysis.predicate import PatientCategories

>>> abn_putamen = hpotk.TermId.from_curie("HP:0031982")
>>> abn_caudate_nucleus = hpotk.TermId.from_curie("HP:0002339")
>>> term_set = {abn_putamen, abn_caudate_nucleus}

>>> config = CohortAnalysisConfiguration()
>>> config.pval_correction = 'fdr_bh'
>>> config.specify_terms_strategy(specified_term_set=term_set)
>>> analysis = configure_cohort_analysis(cohort, hpo, config=config)

This will cause genophenocorr to perform only two hypothesis tests, one for each of the two terms.


heuristic sampler strategy
~~~~~~~~~~~~~~~~~~~~~~~~~~

This method is made up of eight individual heuristics that are explained below. To apply this strategy, enter the following commands.



>>> analysis_config = CohortAnalysisConfiguration()
>>> analysis_config.heuristic_strategy(threshold_HPO_observed_frequency=0.5)
>>> analysis = configure_cohort_analysis(cohort, hpo, config=analysis_config)


This method applies a number of different strategies to skip terms that are unlikely to yield significant or interesting results.

#. Skipping top level terms 
    All the direct children of the root phenotype term `Phenotypic abnormality (HP:0000118) <https://hpo.jax.org/browse/term/HP:0000118>`_ are skipped, because of the assumption that if there is a valid signal, it will derive from one of the more specific descendents. For instance, `Abnormality of the nervous system (HP:0000707) <https://hpo.jax.org/browse/term/HP:0000707>`_ is a child of *Phenotypic abnormality*, and this assumption implies that we think that if there is a signal from the nervous system, it will lead to at least one of the descendents of *Abnormality of the nervous system* being singnificant.

#. Skipping terms that are not descendents of `Phenotypic abnormality (HP:0000118) <https://hpo.jax.org/browse/term/HP:0000118>`_
    The HPO has a number of other branches that describe modes of inheritance, past medical history, and clinical modifiers. We do not think it makes much sense to test for enrichment of these terms, and so they are filtered out.

#. Skip terms with less than threshold observations
    The ``threshold_HPO_observed_frequency`` argument (see above) determines the mininum proportion of individuals. We check each of the genotype groups (e.g., MISSENSE vs. not-MISSENSE), and only retain a term for testing if the proportion of individuals in at least one genotype group is at least ``threshold_HPO_observed_frequency``. This is because of our assumption that even if there is statistical significance, if a term is only seen in (for example) in 7% of individuals in the MISSENSE group and 2% in the not-MISSENSE group, the term is unlikely to be of great interest because it rare.

#. Skip terms if no cell has more than one count 
    In a related heuristic, we skip terms if no genotype group has more than one count. This is not completely redundant with the previous condition, because some terms may have a small number of total observations.

#. Skip terms if genotypes have same HPO proportions
    If both (or all) of the genotype groups have the same proportion of individuals observed to be annotated to an HPO term, e.g., both are 50%, then skip the term, because it is not possible that the Fisher exact test will return a significant result.

#. Skip terms if there are no HPO observations in a group
    If one of the genotype groups has neither observed nor excluded observations for an HPO term, skip it.

#. Skip terms if all counts are identical to counts for a child term.
    Let's say a term such as `Posterior polar cataract (HP:0001115) <https://hpo.jax.org/browse/term/HP:0001115>`_ was observed in 7 of 11 individuals with a MISSENSE variant and in 3 of 8 individuals with NONSENSE variants. If we then test the parent term `Polar cataract HP:0010696 <https://hpo.jax.org/browse/term/HP:0010696>`_ and also find the same counts (7 of 11 and 3 of 8), the we choose not to test it. This is because the more specific an HPO term is, the more information it has (the more interesting the correlation would be if it exists), and the  result of the Fisher Exact Test for *Polar cataract* would of course be exactly the same as for *Posterior polar cataract*.

#. Second-level terms 
    For the same reason that the first level term are judged to be less relevant than more specific terms, the Second-level terms are also assumed by this heuristic to be less interesting. For instance, a typical first-level term in the HPO would be `Abnormal eye morphology HP:0012372 <https://hpo.jax.org/browse/term/HP:0012372>`_, and a typical second-level term would be `Abnormal eye morphology HP:0012372 <https://hpo.jax.org/browse/term/HP:0012372>`_, but an interesting term would be one such as *Posterior polar cataract*. With this heuristic,  if any child of a second-level term  has at least 75% of the terms, then do not test because it is unlikely to add much insight. Only of the second level term has a lot more annotations do we test (this might mean that many different  specific forms of *Abnormal eye morphology* are overrepresented).

    
