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


.. _phenotype-group-stats:

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