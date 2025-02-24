.. _measurement-stat:


====================
Compare measurements
====================


****************
Example analysis
****************

We will analyze 69 individuals reported by
`Xu et al. (2019) <https://pubmed.ncbi.nlm.nih.gov/30968594/>`_.
The reported data includes the *CYP21A2* mutations, the lab measurement,
and other clinical signs and symptoms encoded as HPO terms.

TODO: wordsmith


Load cohort
===========

For the purpose of this analysis, we will load the :class:`~gpsea.model.Cohort`
from a `JSON file <https://github.com/P2GX/gpsea/tree/main/docs/cohort-data/CYP21A2.0.1.20.json>`_.
The cohort was prepared from phenopackets as described in :ref:`create-a-cohort` section,
and then serialized as a JSON file following the instructions in :ref:`cohort-persistence` section.

.. 
   Prepare the JSON file by running the tests in `tests/tests/test_generate_doc_cohorts.py`.

>>> import json
>>> from gpsea.io import GpseaJSONDecoder
>>> fpath_cohort_json = 'docs/cohort-data/CYP21A2.0.1.20.json'
>>> with open(fpath_cohort_json) as fh:
...     cohort = json.load(fh, cls=GpseaJSONDecoder)
>>> len(cohort)
69


Configure analysis
==================

*MANE* transcript of *CYP21A2*. 

>>> tx_id = 'NM_000500.9'


Genotype predicate
------------------

For now, just Missense vs rest.
Nonsense variants + big SVs,
Missense should *NOT* be severe.
TODO - create real predicate.

>>> from gpsea.model import VariantEffect
>>> from gpsea.analysis.predicate import variant_effect
>>> is_missense = variant_effect(VariantEffect.MISSENSE_VARIANT, tx_id=tx_id)
>>> is_missense.description
'MISSENSE_VARIANT on NM_000500.9'

Assuming AR inheritance, we compare missense vs. rest:

>>> from gpsea.analysis.clf import biallelic_classifier
>>> gt_clf = biallelic_classifier(
...     a_predicate=is_missense,
...     b_predicate=~is_missense,
...     a_label="Missense", b_label="Other",
...     partitions=({0,}, {1, 2}),
... )
>>> gt_clf.class_labels
('Missense/Missense', 'Missense/Other OR Other/Other')

Phenotype score
---------------

We investigate an association between a measurement value and a genotype group.

We use the measurement of `Testosterone [Mass/volume] in Serum or Plasma <https://loinc.org/2986-8/>`_
(`LOINC:2986-8`).

>>> from gpsea.analysis.pscore import MeasurementPhenotypeScorer
>>> testosterone = 'LOINC:2986-8'
>>> pheno_scorer = MeasurementPhenotypeScorer.from_measurement_id(
...     term_id=testosterone,
...     label="Testosterone [Mass/volume] in Serum or Plasma",
... )
>>> pheno_scorer.description
'Value of Testosterone [Mass/volume] in Serum or Plasma [LOINC:2986-8]'


Statistical test
----------------

We will use T test to test for differences between scores
of the different genotype groups

>>> from gpsea.analysis.pscore.stats import TTestStatistic
>>> score_statistic = TTestStatistic()


Final analysis
--------------

We will put the final analysis together into :class:`~gpsea.analysis.pscore.PhenotypeScoreAnalysis`.

>>> from gpsea.analysis.pscore import PhenotypeScoreAnalysis
>>> score_analysis = PhenotypeScoreAnalysis(
...     score_statistic=score_statistic,   
... )


Analysis
========

We execute the analysis by running

>>> result = score_analysis.compare_genotype_vs_phenotype_score(
...     cohort=cohort,
...     gt_clf=gt_clf,
...     pheno_scorer=pheno_scorer,
... )

>>> result.pval
0.741216622359659

Show data frame with scores

>>> scores = result.data.sort_index()
>>> scores.head()  # doctest: +NORMALIZE_WHITESPACE
                                            genotype  phenotype
patient_id                                                     
individual 10[PMID_30968594_individual_10]         1      614.0
individual 11[PMID_30968594_individual_11]         1      630.0
individual 12[PMID_30968594_individual_12]         1        NaN
individual 13[PMID_30968594_individual_13]         1      303.0
individual 14[PMID_30968594_individual_14]         1      664.0


Prepare genotype category legend:

>>> gt_id_to_name = {c.category.cat_id: c.category.name for c in gt_clf.get_categorizations()}
>>> gt_id_to_name
{0: 'Missense/Missense', 1: 'Missense/Other OR Other/Other'}

TODO: wordsmith & finish!

