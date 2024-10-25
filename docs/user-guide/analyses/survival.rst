.. _survival:

=================
Survival analysis
=================


****************
Example analysis
****************

We will analyze the time until end stage renal disease in 207 individuals with mutations in *UMOD*.
Specifically, we will test for difference between the onset of the end stage renal disease in the individuals with mutation
in exon 3 of *UMOD* vs. individuals with other *UMOD* mutation.


Load cohort
===========

For the purpose of this analysis, we will load the :class:`~gpsea.model.Cohort`
from a `JSON file <https://github.com/monarch-initiative/gpsea/tree/main/docs/cohort-data/UMOD.0.1.20.json>`_.
The cohort was prepared from phenopackets as described in :ref:`create-cohort-from-phenopackets` section,
and then serialized as a JSON file following the instructions in :ref:`cohort-persistence` section.

.. 
   Prepare the JSON file by running the tests in `tests/tests/test_generate_doc_cohorts.py`.

>>> import json
>>> from gpsea.io import GpseaJSONDecoder
>>> fpath_cohort_json = 'docs/cohort-data/UMOD.0.1.20.json'
>>> with open(fpath_cohort_json) as fh:
...     cohort = json.load(fh, cls=GpseaJSONDecoder)
>>> len(cohort)
207


Configure analysis
==================

*MANE* transcript of *UMOD*.

>>> tx_id = 'NM_003361.4'

Genotype predicate
------------------

One allele of exon 3 vs. one allele of elsewhere.

>>> from gpsea.analysis.predicate.genotype import VariantPredicates
>>> is_in_exon3 = VariantPredicates.exon(exon=3, tx_id=tx_id)
>>> is_in_exon3.get_question()
'variant affects exon 3 on NM_003361.4'

Monoallelic predicate to compare one allele of *UMOD* exon 3 variant
versus one allele of other *UMOD* variant:

>>> from gpsea.analysis.predicate.genotype import monoallelic_predicate
>>> gt_predicate = monoallelic_predicate(
...     a_predicate=is_in_exon3,
...     b_predicate=~is_in_exon3,
...     a_label="Exon 3", b_label="Other",
... )
>>> gt_predicate.group_labels
('Exon 3', 'Other')


Survival endpoint
-----------------

The endpoint of our study is defined as development of end stage renal disease.
In the *UMOD* cohort, this is encoded with
`Stage 5 chronic kidney disease <https://hpo.jax.org/browse/term/HP:0003774>`_
(`HP:0003774`) HPO term.
We need to leverage the HPO hierarchy when computing
the onset of an HPO term. Let's load HPO:

>>> import hpotk
>>> store = hpotk.configure_ontology_store()
>>> hpo = store.load_minimal_hpo(release='v2024-07-01')

and now we can create an :class:`~gpsea.analysis.temporal.Endpoint`
to compute the time until an individual develops end stage renal disease:

>>> from gpsea.analysis.temporal.endpoint import hpo_onset
>>> term_id = "HP:0003774"  # Stage 5 chronic kidney disease
>>> endpoint = hpo_onset(hpo=hpo, term_id=term_id)
>>> endpoint.summary
'Compute time until onset of Stage 5 chronic kidney disease'


Statistical test
----------------

We will use Log rank test to compare the age until the endpoint between
the genotype groups:

>>> from gpsea.analysis.temporal.stats import LogRankTest
>>> survival_statistic = LogRankTest()

Final analysis
--------------

We will put the final analysis together into :class:`~gpsea.analysis.temporal.PhenotypeScoreAnalysis`.

>>> from gpsea.analysis.temporal import SurvivalAnalysis
>>> survival_analysis = SurvivalAnalysis(
...     statistic=survival_statistic,
... )

Analysis
========

We execute the analysis by running

>>> result = survival_analysis.compare_genotype_vs_survival(
...     cohort=cohort,
...     gt_predicate=gt_predicate,
...     endpoint=endpoint,
... )

>>> result.pval
0.06200425830044376


Kaplan-Meier curves
-------------------


We can plot Kaplan-Meier curves:

>>> from gpsea.model import Age
>>> import matplotlib as mpl
>>> import matplotlib.pyplot as plt
>>> fig, ax = plt.subplots(figsize=(6, 4), dpi=120)
>>> result.plot_kaplan_meier_curves(
...     ax=ax,
... )
>>> _ = ax.xaxis.set(
...     # Show X axis in years ...
...     major_formatter=mpl.ticker.FuncFormatter(lambda x, pos: f"{x / Age.DAYS_IN_YEAR:.0f}"),  
...     # ... with a tick for every decade
...     major_locator=mpl.ticker.MultipleLocator(10 * Age.DAYS_IN_YEAR),
... )
>>> _ = ax.set(
...     xlabel=endpoint.name + " [years]",
...     ylabel="Empirical survival",
... )
>>> _ = ax.grid(axis="y")

.. image:: /img/umod_km_curves.png
   :alt: UMOD Kaplan-Meier curves
   :align: center
   :width: 600px

.. doctest:: survival
   :hide:

   >>> fig.savefig('docs/img/umod_km_curves.png')  # doctest: +SKIP
   

Raw data
--------

The `result` includes the survival values for all cohort members:

>>> survivals = result.data.sort_index()
>>> survivals.head()  # doctest: +NORMALIZE_WHITESPACE
                          genotype    phenotype
patient_id                                                                        
AII.1[PMID_22034507_AII_1]       0    Survival(value=18262.5, is_censored=True)
AII.2[PMID_22034507_AII_2]       0    None
AII.3[PMID_22034507_AII_3]       0    Survival(value=16436.25, is_censored=True)
AII.5[PMID_22034507_AII_5]       0    Survival(value=22280.25, is_censored=False)
AIII.4[PMID_22034507_AIII_4]     0    Survival(value=19723.5, is_censored=False)

Each line corresponeds to an individual and the dataframe is indexed by the individual's identifier/label.
The `genotype` column contains the genotype group code,
and `phenotype` column includes a :class:`~gpsea.analysis.temporal.Survival` value
or `None` if computing the survival was impossible (see :func:`~gpsea.analysis.temporal.endpoint.hpo_onset` for details).
The `Survival` reports the number of days until attaining the endpoint,
here defined as end stage renal disease (`is_censored=False`),
or until the individual dropped out of the analysis (`is_censored=True`).
