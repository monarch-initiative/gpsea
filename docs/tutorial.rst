.. _tutorial:

========
Tutorial
========

Here we demonstrate an end-to-end genotype-phenotype analysis with GPSEA.
The tutorial illustrates just one of the many ways GPSEA can be used to characterize genotype-phenotype correlations.
The :ref:`user-guide` contains details about additional methods and functionalities.


The tutorial presents an analysis of a cohort of individuals with pathogenic variants in *TBX5* leading to
`Holt-Oram syndrome MIM:142900 <https://omim.org/entry/142900>`_.

Holt-Oram syndrome is an autosomal dominant disorder characterized by
upper limb defects, congenital heart defects, and arrhythmias (`PMID:38336121 <https://pubmed.ncbi.nlm.nih.gov/38336121/>`_).
It has been observed in the literature that congenital defects of the ventricular and atrial septum are more
common in the truncating than in the missense variants (`PMID:30552424 <https://pubmed.ncbi.nlm.nih.gov/30552424/>`_).
Additionally, upper limb defects are more frequent in patients with protein-truncating variants (`PMID:38336121 <https://pubmed.ncbi.nlm.nih.gov/38336121/>`_).

To perform the analysis, we curated the literature and created on `GA4GH phenopacket <https://pubmed.ncbi.nlm.nih.gov/35705716/>`_ for each
affected individual. The phenopackets are made available in `Phenopacket Store <https://github.com/monarch-initiative/phenopacket-store>`_.



The analysis
~~~~~~~~~~~~

For the analysis, the `MANE <https://www.ncbi.nlm.nih.gov/refseq/MANE/>`_ transcript
(i.e., the "main" biomedically relevant transcript of a gene) should be chosen unless
there is a specific reason not to (which should occur rarely if at all).

In the case of *TBX5* the MANE transcript is `NM_181486.4`. Note that the trascript identifier (`NM_181486`) and the version (`4`) are both required.
A good way to find the MANE transcript is to search on the gene symbol (e.g., *TBX5*) in `ClinVar <https://www.ncbi.nlm.nih.gov/clinvar/>`_ and to
choose a variant that is specifically located in the gene. The MANE transcript will be displayed here (e.g., `NM_181486.4(TBX5):c.1221C>G (p.Tyr407Ter)
<https://www.ncbi.nlm.nih.gov/clinvar/variation/495227/>`_).

We additionally need the corresponding protein identifier.
A good way to find this is to search on the transcript id in `NCBI Nucleotide <https://www.ncbi.nlm.nih.gov/nuccore/>`_.
In our case, search on `NM_181486.4` will bring us to `this page <https://www.ncbi.nlm.nih.gov/nuccore/NM_181486.4>`_.
If we search within this page for `"NP_"`, this will bring us to the
corresponding protein accession `NP_852259.1`.

>>> cohort_name = 'TBX5'
>>> tx_id = 'NM_181486.4'
>>> px_id = 'NP_852259.1'


Load HPO
^^^^^^^^

GPSEA needs HPO to do the analysis.
We use HPO toolkit to load HPO version `v2024-07-01`:

>>> import hpotk
>>> store = hpotk.configure_ontology_store()
>>> hpo = store.load_minimal_hpo(release='v2024-07-01')

.. tip::

  Use the latest HPO release by omitting the `release` option.

Prepare cohort
^^^^^^^^^^^^^^

Now we will load the samples to analyze. We will use the cohort of 156 individuals with mutations in *TBX5*
whose clinical signs and symptoms were encoded into HPO terms
and stored in `Phenopacket Store <https://github.com/monarch-initiative/phenopacket-store>`_.

>>> from ppktstore.registry import configure_phenopacket_registry
>>> phenopacket_registry = configure_phenopacket_registry()
>>> with phenopacket_registry.open_phenopacket_store('0.1.18') as ps:
...     phenopackets = tuple(ps.iter_cohort_phenopackets(cohort_name))
>>> len(phenopackets)
156

We loaded 156 phenopackets which need further preprocessing to prepare for the analysis.
We will compute functional annotations for the mutations and then include the individuals into
a :class:`~gpsea.model.Cohort`:

>>> from gpsea.preprocessing import configure_caching_cohort_creator, load_phenopackets
>>> cohort_creator = configure_caching_cohort_creator(hpo)
>>> cohort, validation = load_phenopackets(  # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
...     phenopackets=phenopackets,
...     cohort_creator=cohort_creator,
... )
Patients Created: ...

and we will check that there are no Q/C issues:

>>> validation.summarize()  # doctest: +SKIP
Validated under none policy
No errors or warnings were found

We loaded the patient data into a `cohort` which is ready for the next steps.

.. seealso::

  Here we show how to create a :class:`~gpsea.model.Cohort` from phenopackets.
  See :ref:`input-data` section to learn how to create a cohort from another inputs.


Explore cohort
^^^^^^^^^^^^^^

We can now explore the cohort to see how many patients are included.

>>> from gpsea.view import CohortViewable
>>> viewer = CohortViewable(hpo)
>>> report = viewer.process(cohort=cohort, transcript_id=tx_id)
>>> with open('docs/report/tbx5_cohort_info.html', 'w') as fh:  # doctest: +SKIP
...     _ = fh.write(report)

.. raw:: html
  :file: report/tbx5_cohort_info.html

.. note::

  The report can also be displayed directly in a Jupyter notebook by running::

    from IPython.display import HTML, display
    display(HTML(report))

Now we can show the distribution of variants with respect to the encoded protein.
We first obtain `tx_coordinates` (:class:`~gpsea.model.TranscriptCoordinates`)
and `protein_meta` (:class:`~gpsea.model.ProteinMetadata`)
with information about the transcript and protein "anatomy":

>>> from gpsea.model.genome import GRCh38
>>> from gpsea.preprocessing import configure_protein_metadata_service, VVMultiCoordinateService
>>> txc_service = VVMultiCoordinateService(genome_build=GRCh38)
>>> pms = configure_protein_metadata_service()
>>> tx_coordinates = txc_service.fetch(tx_id)
>>> protein_meta = pms.annotate(px_id)

and we follow with plotting the diagram of the mutations on the protein:

>>> from gpsea.view import ProteinVisualizer
>>> import matplotlib.pyplot as plt
>>> fig, ax = plt.subplots(figsize=(15, 8))
>>> visualizer = ProteinVisualizer()
>>> visualizer.draw_protein_diagram(
...     tx_coordinates,
...     protein_meta,
...     cohort,
...     ax=ax,
... )
>>> fig.tight_layout()
>>> fig.savefig('docs/img/tutorial/tbx5_protein_diagram.png')  # doctest: +SKIP

.. image:: /img/tutorial/tbx5_protein_diagram.png
   :alt: TBX5 protein diagram
   :align: center
   :width: 600px


Prepare genotype and phenotype predicates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We will create a predicate to bin patients into group
depending on presence of a missense and frameshift variant to test
if there is a difference between frameshift and non-frameshift variants
in the individuals of the *TBX5* cohort.

>>> from gpsea.model import VariantEffect
>>> from gpsea.analysis.predicate.genotype import VariantPredicates, groups_predicate
>>> gt_predicate = groups_predicate(
...     predicates=(
...         VariantPredicates.variant_effect(VariantEffect.MISSENSE_VARIANT, tx_id),
...         VariantPredicates.variant_effect(VariantEffect.FRAMESHIFT_VARIANT, tx_id)
...     ),
...     group_names=('Missense', 'Frameshift'),
... )
>>> gt_predicate.display_question()
'Genotype group: Missense, Frameshift'

.. note::

  There are many other ways to set up a predicate for testing
  for a GP correlation.
  See the :ref:`predicates` section to learn more about building
  a predicate of interest.

The phenotype grouping is based on presence or absence of an HPO term.
We take advantage of the ontology "true path rule" to infer presence
of the ancestor terms for all present HPO terms
(e.g. presence of `Abnormal ventricular septum morphology <https://hpo.jax.org/browse/term/HP:0010438>`_
in an individual annotated with `Ventricular septal defect <https://hpo.jax.org/browse/term/HP:0001629>`_)
and exclusion of the descendant terms for all excluded terms (e.g. exclusion of
`Motor seizure <https://hpo.jax.org/browse/term/HP:0020219>`_
in an individual where `Seizure <https://hpo.jax.org/browse/term/HP:0001250>`_
was excluded):

>>> from gpsea.analysis.predicate.phenotype import prepare_predicates_for_terms_of_interest
>>> pheno_predicates = prepare_predicates_for_terms_of_interest(
...     cohort=cohort,
...     hpo=hpo,
...     min_n_of_patients_with_term=2,
... )

By default, GPSEA will perform one hypothesis test for each HPO term used to annotate two or more individuals in the cohort
(see ``min_n_of_patients_with_term=2`` above).
Testing multiple hypothesis on the same dataset increases the chance of receiving false positive result.
However, GPSEA simplifies the application of an appropriate multiple testing correction.

For general use, we recommend using a combination
of a *Phenotype MTC filter* (:class:`~gpsea.analysis.PhenotypeMtcFilter`) with a *multiple testing correction*.
Phenotype MTC filter chooses the HPO terms to test according to several heuristics, which
reduce the multiple testing burden and focus the analysis
on the most interesting terms (see :ref:`HPO MTC filter <hpo-mtc-filter-strategy>` for more info).
Then the multiple testing correction, such as Bonferroni or Benjamini-Hochberg,
is used to control the family-wise error rate or the false discovery rate.
See :ref:`mtc` for more information.

In this example, we will use a combination of the HPO MTC filter (:class:`~gpsea.analysis.mtc_filter.HpoMtcFilter`)
with Benjamini-Hochberg procedure (``mtc_correction='fdr_bh'``)
with a false discovery control level at (``mtc_alpha=0.05``):

>>> from gpsea.analysis.mtc_filter import HpoMtcFilter
>>> mtc_filter = HpoMtcFilter.default_filter(hpo, term_frequency_threshold=0.2)
>>> mtc_correction = 'fdr_bh'
>>> mtc_alpha = 0.05

Choosing the statistical procedure for assessment of association between genotype and phenotype
groups is the last missing piece of the analysis. We will use Fisher Exact Test:

>>> from gpsea.analysis.pcats.stats import FisherExactTest
>>> count_statistic = FisherExactTest()

and we finalize the analysis setup by putting all components together
into :class:`~gpsea.analysis.pcats.HpoTermAnalysis`:

>>> from gpsea.analysis.pcats import HpoTermAnalysis
>>> analysis = HpoTermAnalysis(
...     count_statistic=count_statistic,
...     mtc_filter=mtc_filter,
...     mtc_correction=mtc_correction,
...     mtc_alpha=mtc_alpha,
... )

Now we can perform the analysis and investigate the results.

>>> result = analysis.compare_genotype_vs_phenotypes(
...     cohort=cohort,
...     gt_predicate=gt_predicate,
...     pheno_predicates=pheno_predicates,
... )
>>> result.total_tests
17

We only tested 1y HPO terms. This is despite the individuals being collectively annotated with
260 direct and indirect HPO terms

>>> len(result.phenotypes)
260

We can show the reasoning behind *not* testing 243 (`260 - 17`) HPO terms
by exploring the phenotype MTC filtering report.

>>> from gpsea.view import MtcStatsViewer
>>> mtc_viewer = MtcStatsViewer()
>>> mtc_report = mtc_viewer.process(result)
>>> with open('docs/report/tbx5_frameshift_vs_missense.mtc_report.html', 'w') as fh:  # doctest: +SKIP
...     _ = fh.write(mtc_report)

.. raw:: html
  :file: report/tbx5_frameshift_vs_missense.mtc_report.html

and these are the top 20 HPO terms ordered by the p value corrected with the Benjamini-Hochberg procedure:

>>> from gpsea.view import summarize_hpo_analysis
>>> summary_df = summarize_hpo_analysis(hpo, result)
>>> summary_df.head(20).to_csv('docs/report/tbx5_frameshift_vs_missense.csv')  # doctest: +SKIP

.. csv-table:: *TBX5* frameshift vs missense
   :file: report/tbx5_frameshift_vs_missense.csv
   :header-rows: 2

We see that several HPO terms are significantly associated
with presence of a frameshift variant in *TBX5*.
For example, `Ventricular septal defect <https://hpo.jax.org/browse/term/HP:0001629>`_
was observed in 31/60 (52%) patients with a missense variant
but it was observed in 19/19 (100%) patients with a frameshift variant.
Fisher exact test computed a p value of `~0.0000562`
and the p value corrected by Benjamini-Hochberg procedure
is `~0.000955`.
