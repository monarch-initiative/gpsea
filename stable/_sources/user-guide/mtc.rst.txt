.. _mtc:

===========================
Multiple-testing correction
===========================

**********
Background
**********

A p-value is the probability that a test result, under the null hypothesis, 
assumes the observed or a more extreme value. It is important to realize that if we
perform many tests, we are likely to get a "significant" result by chance alone. 
For instance, if we test a null hypothesis that is true using a significance level 
of :math:`\alpha = 0.05`, then there is a probability of :math:`1-\alpha = 0.95` 
of arriving at a correct conclusion of non-significance. If we now test
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


***********************
Implementation in GPSEA
***********************

By default, GPSEA performs a hypothesis test for each HPO term found at least twice
in the cohort, meaning that we may perform up to hundreds of tests.
Therefore, unless we take into account the fact that multiple statistical tests are being performed,
it is likely that we will obtain one or more false-positive results.

GPSEA offers two approaches to mitigate this problem: multiple-testing correction (MTC) procedures
and MTC filters to choose the terms to be tested.


Multiple-testing correction procedures
======================================

A number of MTC procedures have
been developed to limit the probability of false-positive results. The
MTC procedures differ in complexity, in their assumptions about the
data, and in the type of control they provide.

The GPSEA package uses the Python package `statsmodels <https://www.statsmodels.org/devel/>`_ to implement
MTC. See the `documentation <https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html>`_ for details;
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
correction (``bonferroni``). The Bonferroni procedure thus provides control of the family-wise
error rate (FWER), which is the probability of at least one Type I
error.  The Bonferroni method multiplies the p-value
returned by each test (which is call the *nominal* p-value)
by the number of tests performed (the result is capped at 1.0).

Alternatively, procedures that control the false-discovery rate (FDR),
limit the proportion of significant results that are type I
errors (false discoveries). 
The Benjamini and Hochberg method (``fdr_bh``) is probably the most commonly used one.
This is the *default* method in GPSEA.

To set an alternative MTC procedure, we use the `mtc_correction` option
when creating an instance of :class:`~gpsea.analysis.pcats.HpoTermAnalysis`:

>>> from gpsea.analysis.mtc_filter import UseAllTermsMtcFilter
>>> from gpsea.analysis.pcats import HpoTermAnalysis
>>> from gpsea.analysis.pcats.stats import FisherExactTest
>>> analysis = HpoTermAnalysis(
...     count_statistic=FisherExactTest(),
...     mtc_filter=UseAllTermsMtcFilter(),
...     mtc_correction='bonferroni',  #      <--- The MTC correction setup
... )


.. _mtc-filters:

MTC filters: Choosing which terms to test
=========================================

We can reduce the overall MTC burden by choosing which terms to test. 
For example, if we choose to test only ten terms out of 450, 
then the mutliplication factor of the Bonferroni correction 
is only 10 instead of 450, and more p-values 
may "survive" the multiple-testing correction.

In the context of GPSEA, we represent the concept of phenotype filtering 
by :class:`~gpsea.analysis.mtc_filter.PhenotypeMtcFilter`.
The filter must be chosen before the :class:`~gpsea.analysis.pcats.MultiPhenotypeAnalysis`,
such as :class:`~gpsea.analysis.pcats.HpoTermAnalysis`, is run:

>>> from gpsea.analysis.pcats import HpoTermAnalysis
>>> analysis = HpoTermAnalysis()  # doctest: +ELLIPSIS
Traceback (most recent call last):
  ...
TypeError: HpoTermAnalysis.__init__() missing 2 required positional arguments: 'count_statistic' and 'mtc_filter'

Note the missing `mtc_filter` option.

We describe the three filtering strategies in the following sections.


.. _use-all-terms-strategy:

Test all terms
--------------

The first MTC filtering strategy is the simplest - do not apply any filtering at all.
This will result in testing all terms. We do not recommend this strategy, 
but it can be useful to disable MTC filtering.

The strategy is implemented in :class:`~gpsea.analysis.mtc_filter.UseAllTermsMtcFilter`.

>>> from gpsea.analysis.mtc_filter import UseAllTermsMtcFilter
>>> use_all = UseAllTermsMtcFilter()

.. _specify-terms-strategy:

Specify terms strategy
----------------------

In presence of a specific hypothesis as to which terms may be different between groups, 
then you can specify these terms in :class:`~gpsea.analysis.mtc_filter.SpecifiedTermsMtcFilter`.

For example if we want to specifically test
`Abnormal putamen morphology (HP:0031982) <https://hpo.jax.org/browse/term/HP:0031982>`_ and
`Abnormal caudate nucleus morphology (HP:0002339) <https://hpo.jax.org/browse/term/HP:0002339>`_
we pass an iterable (e.g. a tuple) with these two terms as an argument:

>>> from gpsea.analysis.mtc_filter import SpecifiedTermsMtcFilter
>>> specified_terms = SpecifiedTermsMtcFilter(
...     terms_to_test=(
...         "HP:0031982",  # Abnormal putamen morphology
...         "HP:0002339",  # Abnormal caudate nucleus morphology
...     )
... )
>>> len(specified_terms.terms_to_test)
2


.. _hpo-mtc-filter-strategy:

HPO MTC filter strategy
-----------------------

Last, the HPO MTC strategy involves making several domain judgments to take advantage of the HPO structure.
The strategy needs access to HPO:

>>> import hpotk
>>> store = hpotk.configure_ontology_store()
>>> hpo = store.load_minimal_hpo(release='v2024-07-01')

and it is implemented in the :class:`~gpsea.analysis.mtc_filter.HpoMtcFilter` class:

>>> from gpsea.analysis.mtc_filter import HpoMtcFilter
>>> hpo_mtc = HpoMtcFilter.default_filter(hpo=hpo)


We use static constructor :func:`~gpsea.analysis.mtc_filter.HpoMtcFilter.default_filter`
for creating :class:`~gpsea.analysis.mtc_filter.HpoMtcFilter`.
The constructor takes a ``term_frequency_threshold`` option (40% by default) 
and the method's logic is made up of 8 individual heuristics 
designed to skip testing the HPO terms that are unlikely to yield significant or interesting results.


`HMF01` - Skip terms that occur very rarely
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``term_frequency_threshold`` determines the mininum proportion of individuals
with direct or indirect annotation by the HPO term to test.
We check each of the genotype groups (e.g., MISSENSE vs. not-MISSENSE),
and we only retain a term for testing if the proportion of individuals
in at least one genotype group is greater than
or equal to ``term_frequency_threshold``.
This is because of our assumption that even if there is statistical significance,
if a term is only seen in (for example) 7% of individuals
in the MISSENSE group and 2% in the not-MISSENSE group,
the term is unlikely to be of great interest because it is rare.


`HMF02` - Skip terms if no genotype group has more than one count
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In a related heuristic, we skip terms if no genotype group has more 
than one count. This is not completely redundant with the previous condition,
because some terms may have a small number of total observations.


`HMF03` - Skip terms if all counts are identical to counts for a child term
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Let's say a term such as    
`Posterior polar cataract (HP:0001115) <https://hpo.jax.org/browse/term/HP:0001115>`_
was observed in 7 of 11 individuals with MISSENSE variants
and in 3 of 8 individuals with NONSENSE variants.
If we find the same patient counts (7 of 11 and 3 of 8) in the parent term
`Polar cataract HP:0010696 <https://hpo.jax.org/browse/term/HP:0010696>`_,
then we choose to not test the parent term.
                                                                                         
This is because the more specific an HPO term is,
the more information it has (the more interesting the correlation would be if it exists),
and the result of a test, such as the Fisher Exact test, would be exactly the same
for *Polar cataract* as for *Posterior polar cataract*.


`HMF05` - Skip term if one of the genotype groups has neither observed nor excluded observations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Skip terms if there are no HPO observations in a group.


`HMF06` - Skip term if underpowered for 2x2 or 2x3 analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If the individuals are binned into 2 phenotype groups and 2 genotype groups (2x2)
and the total count of patients in all genotype-phenotype groups is less than 7,
or into 2 phenotype groups and 3 genotype groups (2x3) and the total count of patients
is less than 6, then there is a lack even of the nominal statistical power
and the counts can never be significant.


`HMF07` - Skipping terms that are not descendents of *Phenotypic abnormality*
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The HPO has a number of other branches that describe modes of inheritance,
past medical history, and clinical modifiers.
We do not think it makes much sense to test for enrichment of these terms,
so, all terms that are not descendants of
`Phenotypic abnormality <https://hpo.jax.org/browse/term/HP:0000118>`_ are filtered out.


`HMF08` - Skipping "general" level terms
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

All the direct children of the root phenotype term
`Phenotypic abnormality (HP:0000118) <https://hpo.jax.org/browse/term/HP:0000118>`_
are skipped, because of the assumption that if there is a valid signal, 
it will derive from one of the more specific descendents.

For instance,
`Abnormality of the nervous system <https://hpo.jax.org/browse/term/HP:0000707>`_
(HP:0000707) is a child of *Phenotypic abnormality*, and this assumption implies
that if there is a signal from the nervous system,
it will lead to at least one of the descendents of
*Abnormality of the nervous system* being significant.


`HMF09` - Skipping terms that are rare on the cohort level 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We skip terms that occur in less than a certain percentage of cohort members.
The purpose of this threshold is to omit terms for which we simply
do not have much data overall.

For instance, if the cohort consists of 100 individuals,
and we have explicit observed observations for 20 and excluded for 10 individuals,
then the annotation frequency is `0.3`. 

The threshold is set as ``annotation_frequency_threshold`` option
of the :func:`~gpsea.analysis.mtc_filter.HpoMtcFilter.default_filter` constructor,
with the default value of `0.4` (40%).


See :ref:`general-hpo-terms` section for details.
