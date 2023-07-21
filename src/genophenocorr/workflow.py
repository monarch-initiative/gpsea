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
genes_in_focus = {'SURF1', 'FBN1'}
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
# TODO - it may be a good idea to create a `cache_datadir/variant` subfolder for the variants data
vac = VariantAnnotationCache(cache_datadir)
vep = VepFunctionalAnnotator()
vfa = VarCachingFunctionalAnnotator(vac, vep)

# Protein metadata
pm = UniprotProteinMetadataService()
# TODO - it may be a good idea to create a `cache_datadir/protein` subfolder for the proteins data
pac = ProteinAnnotationCache(cache_datadir)
pfa = ProtCachingFunctionalAnnotator(pac, pm)

# Assemble the patient creator
pc = PhenopacketPatientCreator(phenotype_creator, vfa, pfa)
cc = PhenopacketCohortCreator(pc)

# Next, load the patients from phenopackets
##phenopackets: typing.List[Phenopacket] = []  # TODO - implement real parsing based on a path above

##patients = [pc.create_patient(pp) for pp in phenopackets]

# Daniel - Is below okay? Or should we remove cohort all together and have is like you do above?
# In principle yes, we can work with a cohort or a sequence of Patients.
# I would not, however, restrict the annotation to a particular tx and protein IDs, I favor getting all available data.
# I also outlined this in other parts of the code.
patientCohort = cc.create_cohort(fpath_phenopackets)

# TODO - implement the rest below
## remove empty variant/hpo patients
## Remove or identify 


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
