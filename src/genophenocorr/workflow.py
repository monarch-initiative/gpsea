# ## Load
# Load the phenopackets from a folder. Later, replace this with a `PhenopacketCohortCreator`
#
import typing

import hpotk
# pyright: reportGeneralTypeIssues=false
from phenopackets import Phenopacket

from genophenocorr.patient import PhenopacketPatientCreator
from genophenocorr.phenotype import PhenotypeCreator
from genophenocorr.protein import UniprotProteinMetadataService, ProteinAnnotationCache, ProtCachingFunctionalAnnotator
from genophenocorr.variant import VarCachingFunctionalAnnotator, VariantAnnotationCache, VepFunctionalAnnotator
from genophenocorr.cohort import PhenopacketCohortCreator

# These are the analysis parameters. We'll expect these parameters to be provided from the user, either in a notebook
# or in CLI.
fpath_hpo = 'path/to/hp.json'
cache_datadir = '/path/to/cache/dir'
fpath_phenopackets = '/path/to/folder/with/phenopackets'
tx_identifier = 'NM_003172.4'  # The ID of the transcript to use for the analysis (e.g. `NM_003172.4` for *SURF1*)
protein_identifier = 'Q15526'  # The ID of the protein to restrict the analysis into (e.g. `Q15526` for *SURF1*)

# Now, let's configure the components for parsing phenopackets, functional annotation, and
# for gathering protein metadata.
# Phenotype Creator
hpo: hpotk.ontology.MinimalOntology = hpotk.ontology.load.obographs.load_minimal_ontology(fpath_hpo)
validators = [
    hpotk.validate.AnnotationPropagationValidator(hpo),
    hpotk.validate.ObsoleteTermIdsValidator(hpo),
    hpotk.validate.PhenotypicAbnormalityValidator(hpo)
]
phenotype_creator = PhenotypeCreator(hpo, hpotk.validate.ValidationRunner(validators))

# Functional annotator
vac = VariantAnnotationCache(cache_datadir)
vep = VepFunctionalAnnotator()
vfa = VarCachingFunctionalAnnotator(vac, vep)

# Protein metadata
pm = UniprotProteinMetadataService()
pac = ProteinAnnotationCache(cache_datadir)
pfa = ProtCachingFunctionalAnnotator(pac, pm)

# Assemble the patient creator
pc = PhenopacketPatientCreator(phenotype_creator, vfa, pfa)
cc = PhenopacketCohortCreator(pc)

# Next, load the patients from phenopackets
##phenopackets: typing.List[Phenopacket] = []  # TODO - implement real parsing based on a path above

##patients = [pc.create_patient(pp) for pp in phenopackets]

# Daniel - Is below okay? Or should we remove cohort all together and have is like you do above? 
patientCohort = cc.create_cohort(fpath_phenopackets, tx_identifier, protein_identifier)

# TODO - implement the rest below

# ## Preview
# function per table


# Functions for tables.
# What comparisons does the user want to do?
# we provide a dict:
# key: variant class
# value: count
#
# - variant effects/types (missense, nonsense, etc)
#   - we do not filter counts
#   - nice column order, start with common variant effects (missense, nonse)
#   - we only show columns with at least one entry
#
# - exon table
#

# a table of variants with at least 1 allele in n or more samples/%

# protein domains
# Ideally, the Uniprot response would return more structured information.
# We can offer a function for basic summary based on protein feature type and on all feature ids/names

# ## Analysis
# (1) predicate is either a function that returns a bool or a named int
# (2) offer a library of predicates
