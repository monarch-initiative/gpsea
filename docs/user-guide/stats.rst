.. _stats:

=================
Statistical tests
=================

There are many different ways of statistically testing for genotype-phenotype correlations, and the appropriate statistical test depends on the question. 
This document provides an overview of the tests offered by the genophenocorr library and explanations of how they 
are implemented by our software.


Fisher Exact Test (FET)
~~~~~~~~~~~~~~~~~~~~~~~

The Fisher exact test (FET) calculates the exact probability value for the
relationship between two dichotomous variables. In our implementation, the two dichotomous variables are the genotype and the phenotype.
For instance, the individuals of the cohort may be divided according to whether or not they have a nonsense variant and according to whether
or not they have ataxia.


The results of FET are expressed in terms of an exact probability (P-value), varying within 0 and 1. Two groups are
considered statistically significant if the P-value is less than the chosen
significance level (usually &alpha;=0.05). 





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
>>> pval = r.pvalue
>>> float(pval)
6.348081479150902e-06

``pval`` evaluates to `6.348081479150901e-06`, meaning there is a significant difference between the groups.