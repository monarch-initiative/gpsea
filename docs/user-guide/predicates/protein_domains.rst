.. _protein_domains:

===============
Protein Domains
===============

Specific domains of a protein may be associated with genotype-phenotype correlations. 
For instance, variants in the pore domain of *PIEZO1* are associated with more severe clinical 
manifestions in dehydrated hereditary stomatocytosis `Andolfo et al.,  2018 <https://pubmed.ncbi.nlm.nih.gov/30187933>`_.

GPSEA uses the protein data in several places: to show distribution of variants with respect to the protein domains
or other features of biological interest, and to group the individuals based on presence of a variant predicted
to affect the protein features.
In all cases, the protein data must be formatted as an instance of :class:`~gpsea.model.ProteinMetadata`
and here we show how to get the data and use it in the analysis.

****************************************
Get the data for the protein of interest
****************************************

The protein data (:class:`~gpsea.model.ProteinMetadata`) can be obtained in several ways,
ordered by the convenience:

* fetched from UniProt REST API
* parsed from a JSON file downloaded from UniProt
* entered manually from a data frame


Fetch data from UniProt REST API
================================

The most convenient way to obtain the protein data is to use a :class:`~gpsea.preprocessing.ProteinMetadataService`.
We recommend using the :func:`~gpsea.preprocessing.configure_default_protein_metadata_service`
to reduce the amount of the associated boiler-plate code:

>>> from gpsea.preprocessing import configure_default_protein_metadata_service
>>> pms = configure_default_protein_metadata_service()


Then, fetching the data for protein accession *NP_852259.1* encoded by the *NM_181486.4* transcript of *TBX5*
is as simple as running:

>>> protein_meta = pms.annotate('NP_852259.1')
>>> protein_meta.protein_id
'NP_852259.1'
>>> protein_meta.protein_length
518
>>> len(protein_meta.protein_features)
2

The `protein_meta` represents the *TBX5* isoform that includes 518 aminoacids and two features of interest,
which we can see on the following screenshot of the UniProt entry for *TBX5*:

.. figure:: img/TBX5_uniprot_features.png
   :alt: *TBX5* (P37173, UniProt entry)
   :align: center
   :width: 800px

   Protein features of *TBX5* (Q99593, UniProt entry)

UniProt shows four protein features:

- the Disordered region (1-46) 
- the Disordered region (250-356) 
- presence of Polar residues (263-299)
- presence of Basic and acidic residues (320-346).


Parse UniProt JSON dump
=======================

In the cases, when the REST API cannot give us the data for a protein of interest,
we can download a JSON file representing the protein features manually,
and load the file into :class:`~gpsea.model.ProteinMetadata`.

To do this, click on the *Download* symbol (see the UniProt screenshot figure above). This will open a dialog
that allows the user to choose the contents of the JSON file. 
Do not change the default option (Features - Domain, Region).
Provided that the file has been saved as `docs/user-guide/data/Q99593.json`,
the ``ProteinMetadata`` can be loaded using :func:`~gpsea.model.ProteinMetadata.from_uniprot_json` function.
Note that you will need to obtain information about the protein name (`label`)
and `protein_length`, but these are shown in the UniProt entry:

>>> from gpsea.model import ProteinMetadata
>>> downloaded_uniprot_json = "docs/user-guide/data/Q99593.json"
>>> protein_meta = ProteinMetadata.from_uniprot_json(
...     protein_id="NP_852259.1",
...     label="transforming growth factor beta receptor 2",
...     uniprot_json=downloaded_uniprot_json,
...     protein_length=518,
... )


Enter features manually
=======================

The information about protein features provided by UniProt entries may not always be complete. 
Here we show how to enter the same information manually, in a custom protein dataframe. 

The frame can be created e.g. by running:

>>> import pandas as pd
>>> domains = [
...    {"region": "Disordered","category": "region", "start": 1, "end": 46, },
...    {"region": "Disordered", "category": "region", "start": 250, "end": 356, },
...    {"region": "Polar residues", "category": "compositional bias", "start": 263, "end": 299, },
...    {"region": "Basic and acidic residues", "category": "compositional bias", "start": 320, "end": 346, },
... ]
>>> df = pd.DataFrame(domains)

The `ProteinMetadata` is then created using :func:`~gpsea.model.ProteinMetadata.from_feature_frame` function:

>>> protein_meta = ProteinMetadata.from_feature_frame(
...     protein_id="NP_852259.1",
...     label="transforming growth factor beta receptor 2",
...     features=df,
...     protein_length=518,
... )


***************************************************
Plot distribution of cohort variants on the protein
***************************************************

Having the protein data on hand, we can plot the distribution of the variants found in the cohort members.
GPSEA leverages Matplotlib to create a diagram with variants and protein features,
to get insights into the cohort and to formulate genotype-phenotype association hypotheses.


Example
=======

Let's plot a distribution of the variants found in *TBX5* cohort of Phenopacket Store.
First, some boiler-plate code is needed to load HPO and the 156 phenopackets

>>> import hpotk
>>> from ppktstore.registry import configure_phenopacket_registry
>>> store = hpotk.configure_ontology_store()
>>> hpo = store.load_minimal_hpo(release="v2024-07-01")
>>> registry = configure_phenopacket_registry()
>>> with registry.open_phenopacket_store("0.1.18") as ps:
...     phenopackets = tuple(ps.iter_cohort_phenopackets("TBX5"))
>>> len(phenopackets)
156

which we load into :class:`~gpsea.model.Cohort`:

>>> from gpsea.preprocessing import configure_caching_cohort_creator, load_phenopackets
>>> creator = configure_caching_cohort_creator(hpo)
>>> cohort, _ = load_phenopackets(  # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
...     phenopackets=phenopackets,
...     cohort_creator=creator,
... )
Individuals Processed: ...

Now we need to fetch the transcript information:

>>> from gpsea.model.genome import GRCh38
>>> from gpsea.preprocessing import VVMultiCoordinateService
>>> txc_service = VVMultiCoordinateService(genome_build=GRCh38)
>>> tx_coordinates = txc_service.fetch("NM_181486.4")


and draw the diagram using :class:`~gpsea.view.ProteinVisualizer`:

>>> import matplotlib.pyplot as plt
>>> from gpsea.view import ProteinVisualizer, ProteinVisualizable
>>> pvis = ProteinVisualizable(tx_coordinates=tx_coordinates, protein_meta=protein_meta, cohort=cohort)
>>> visualizer = ProteinVisualizer()
>>> fig, ax = plt.subplots(figsize=(12, 8), dpi=120)
>>> visualizer.draw_fig(pvis=pvis, ax=ax)
>>> fig.tight_layout()
>>> fig.savefig('docs/user-guide/predicates/img/TBX5_gpsea_with_uniprot_domains.png')  # doctest: +SKIP


.. figure:: img/TBX5_gpsea_with_uniprot_domains.png
   :alt: TBX5 (Q99593, UniProt entry)
   :align: center
   :width: 800px

   GPSEA display of variants and protein features of *TBX5*


.. note::

    We call `fig.savefig` to save the diagram for the purpose of showing it in the user guide.
    You certainly do not need to do it as part of your analysis.


*************************************
Use the protein data to test variants
*************************************

TODO

