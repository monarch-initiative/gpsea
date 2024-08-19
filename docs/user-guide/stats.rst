.. _stats:

=================
Statistical tests
=================

There are many different ways of statistically testing for genotype-phenotype correlations, 
and the appropriate statistical test depends on the question. 
This document provides an overview of the tests offered by the genophenocorr library 
and explanations of how they are implemented by our software.


Fisher exact test (FET)
~~~~~~~~~~~~~~~~~~~~~~~

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


.. _phenotype-score-stats:

Mann-Whitney U Test 
~~~~~~~~~~~~~~~~~~~

We may want to compare the total number of occurences of a specific set of phenotypic features between two different genotype groups.
For instance, `Jordan et al (2018) <https://pubmed.ncbi.nlm.nih.gov/29330883/>`_ found that the total number of structural defects 
of the brain, eye, heart, and kidney and sensorineural hearing loss seen in individuals with point mutations in the Atrophin-1 domain of the RERE gene
is significantly higher than expected based on the number of similar defects seen in individuals with putative loss-of-function variants.
Since there are five potential defects, each individual has a count ranging between 0 and 5. 

We perform a Mann-Whitney U Test (or Wilcoxon Rank-Sum Test) to compare the distribution of such counts between genotype groups.
This is a non-parametric test that compares the medians of the two groups to determine if they come from the same distribution. 

>>> import scipy.stats as stats
>>> group1 = [0, 0, 1, 0, 2, 0, 1, 1, 1, 0, 2, 0, 0, 3, 1, 1, 1, 0]
>>> group2 = [4, 5, 3, 4, 3, 3, 3, 4, 4, 5, 5, 2, 3, 0, 3, 5, 2, 3]
>>> r = stats.mannwhitneyu(x=group1, y=group2, alternative = 'two-sided')
>>> p_value = r.pvalue
>>> float(p_value)
6.348081479150902e-06

``p_value`` evaluates to `6.348081479150901e-06`, meaning there is a significant difference between the groups.


Example
^^^^^^^

Let's now analyze the subjects reported in *Jordan et al*. 
We start by loading the cohort from Phenopacket Store (version `0.1.18`):

>>> from ppktstore.registry import configure_phenopacket_registry
>>> registry = configure_phenopacket_registry()
>>> with registry.open_phenopacket_store(release='0.1.18') as ps:
...     phenopackets = tuple(ps.iter_cohort_phenopackets('RERE'))
>>> len(phenopackets)
19

We loaded 19 phenopackets. 

Now, we need to prepare the phenopackets for using with Genophenocorr.
We will need HPO (version `v2024-07-01`)

>>> import hpotk
>>> store = hpotk.configure_ontology_store()
>>> hpo = store.load_minimal_hpo(release='v2024-07-01')

to create cohort creator

>>> from genophenocorr.preprocessing import configure_caching_cohort_creator
>>> cohort_creator = configure_caching_cohort_creator(hpo)

which we will use to preprocess the cohort

>>> from genophenocorr.preprocessing import load_phenopackets
>>> cohort, _ = load_phenopackets(phenopackets, cohort_creator)  # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
Patients Created: ...
>>> len(cohort)
19

Now we can set up the phenotype and genotype predicates. Jordan et al tests ...

.. todo: improve the text

>>> rere_mane_tx_id = 'NM_001042681.2'

Now let's create a predicate for testing if the variant is a point mutation or a loss of function mutation.
The point mutation predicate is defined as ... 
TODO: improve!

>>> from genophenocorr.model import VariantEffect
>>> from genophenocorr.analysis.predicate.genotype import VariantPredicates
>>> point_mutation_effects = (
...     VariantEffect.MISSENSE_VARIANT,
... )
>>> point_mutation = VariantPredicates.change_length('==', 0) \
...     & VariantPredicates.ref_length('==', 1) \
...     & VariantPredicates.any(VariantPredicates.variant_effect(effect, rere_mane_tx_id) for effect in point_mutation_effects)
>>> point_mutation.get_question()
'((change length == 0 AND ref allele length == 1) AND MISSENSE_VARIANT on NM_001042681.2)'

For the loss of function predicate, these variant effects are considered loss of function:

>>> lof_effects = (
...     VariantEffect.TRANSCRIPT_ABLATION,
...     VariantEffect.FRAMESHIFT_VARIANT,
...     VariantEffect.START_LOST,
...     VariantEffect.STOP_GAINED,
... )
>>> lof_mutation = VariantPredicates.any(VariantPredicates.variant_effect(eff, rere_mane_tx_id) for eff in lof_effects)
>>> lof_mutation.get_question()
'(TRANSCRIPT_ABLATION on NM_001042681.2 OR FRAMESHIFT_VARIANT on NM_001042681.2 OR START_LOST on NM_001042681.2 OR STOP_GAINED on NM_001042681.2)'

The genotype predicate will bin the patient into two groups: a point mutation group or the loss of function group:

>>> from genophenocorr.analysis.predicate.genotype import groups_predicate
>>> gt_predicate = groups_predicate(
...     predicates=(point_mutation, lof_mutation),
...     group_names=('Point', 'LoF'),
... )
>>> gt_predicate.get_question()
'What group does the patient belong to: Point, LoF'

Now phenotype predicate. The authors divide the patients into groups according to the count of structural defects
in these groups:

>>> structural_defects = (
...     'HP:0012443',  # Abnormal brain morphology (Brain anomalies)
...     'HP:0012372',  # Abnormal eye morphology (Eye anomalies)
...     'HP:0001627',  # Abnormal heart morphology (Congenital heart defects)
...     'HP:0012210',  # Abnormal renal morphology (Renal anomalies)
...     'HP:0000407',  # Sensorineural hearing impairment (Sensorineural hearing loss)
... )

Let's run the analysis.

>>> from genophenocorr.analysis import configure_cohort_analysis
>>> analysis = configure_cohort_analysis(
...     cohort, hpo,
... )
>>> result = analysis.compare_genotype_vs_phenotype_group_count(
...     gt_predicate=gt_predicate,   
...     phenotype_group_terms=structural_defects,
... )
>>> round(result.p_value, 9)
0.027066902


We have the counts:

>>> counts = result.genotype_phenotype_scores
>>> counts.head()  # doctest: +NORMALIZE_WHITESPACE
                                     genotype phenotype
patient_id                                             
Subject 10[PMID_27087320_Subject_10]      LoF         0
Subject 1[PMID_27087320_Subject_1]      Point         4
Subject 2[PMID_27087320_Subject_2]       None         4
Subject 2[PMID_29330883_Subject_2]        LoF         1
Subject 3[PMID_27087320_Subject_3]      Point         4


Let's plot the data:

>>> import matplotlib.pyplot as plt
>>> import seaborn as sns
>>> fig, ax = plt.subplots(figsize=(6, 4), dpi=120)
>>> data = counts.loc[counts['genotype'].notna()]
>>> _ = sns.stripplot(
...     hue='genotype', y='phenotype',
...     dodge=True, jitter=.15, palette='dark', alpha=.8,
...     data=data, ax=ax,
... )
>>> _ = ax.grid(axis='y')
>>> _ = ax.set(ylim=(-.5, len(structural_defects) + .5))
>>> fig.savefig('docs/img/phenotype_group_counts.png')  # doctest: +SKIP


.. image:: /img/phenotype_group_counts.png
   :alt: Phenotype group counts
   :align: center
   :width: 600px
