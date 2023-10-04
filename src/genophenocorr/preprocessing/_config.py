import os
import typing

import hpotk

from genophenocorr.model.genome import GRCh37, GRCh38

from ._phenotype import PhenotypeCreator
from ._phenopacket import PhenopacketPatientCreator
from ._api import FunctionalAnnotator
from ._protein import ProteinAnnotationCache, ProtCachingFunctionalAnnotator
from ._uniprot import UniprotProteinMetadataService
from ._variant import VarCachingFunctionalAnnotator, VariantAnnotationCache
from ._vep import VepFunctionalAnnotator


def configure_caching_patient_creator(hpo: hpotk.MinimalOntology,
                                      genome_build: str = 'GRCh38.p13',
                                      validation_runner: typing.Optional[hpotk.validate.ValidationRunner] = None,
                                      cache_dir: typing.Optional[str] = None,
                                      variant_fallback: str = 'VEP',
                                      protein_fallback: str = 'UNIPROT') -> PhenopacketPatientCreator:
    """
    A convenience function for configuring a caching :class:`genophenocorr.preprocessing.PhenopacketPatientCreator`.

    To create the patient creator, we need hpo-toolkit's representation of HPO, the validator

    :param hpo: a HPO instance.
    :param genome_build: name of the genome build to use, choose from `{'GRCh37.p13', 'GRCh38.p13'}`.
    :param validation_runner: an instance of the validation runner.
    :param cache_dir: path to the folder where we will cache the results fetched from the remote APIs or `None`
     if the data should be cached in `.cache` folder in the current working directory.
     In any case, the directory will be created if it does not exist (including non-existing parents).
    :param variant_fallback: the fallback variant annotator to use if we cannot find the annotation locally.
     Choose from ``{'VEP'}`` (just one fallback implementation is available at the moment).
    :param protein_fallback: the fallback protein metadata annotator to use if we cannot find the annotation locally.
     Choose from ``{'UNIPROT'}`` (just one fallback implementation is available at the moment).
    """
    if cache_dir is None:
        cache_dir = os.path.join(os.getcwd(), '.genophenocorr_cache')

    if genome_build == 'GRCh38.p13':
        build = GRCh38
    elif genome_build == 'GRCh37.p13':
        build = GRCh37
    else:
        raise ValueError(f'Unknown build {genome_build}. Choose from [\'GRCh37.p13\', \'GRCh38.p13\']')

    phenotype_creator = _setup_phenotype_creator(hpo, validation_runner)
    functional_annotator = _configure_functional_annotator(cache_dir, variant_fallback, protein_fallback)
    return PhenopacketPatientCreator(build, phenotype_creator, functional_annotator)


def _setup_phenotype_creator(hpo: hpotk.MinimalOntology,
                             validator: typing.Optional[hpotk.validate.ValidationRunner]) -> PhenotypeCreator:
    if validator is None:
        # This will be the default validator
        validator = hpotk.validate.ValidationRunner((
            hpotk.validate.ObsoleteTermIdsValidator(hpo),
            hpotk.validate.AnnotationPropagationValidator(hpo),
            hpotk.validate.PhenotypicAbnormalityValidator(hpo)
        ))
    else:
        validator = hpotk.util.validate_instance(validator, hpotk.validate.ValidationRunner, 'validator')
    return PhenotypeCreator(hpo, validator)


def _configure_functional_annotator(cache_dir: str,
                                    variant_fallback: str,
                                    protein_fallback: str) -> FunctionalAnnotator:
    # (1) ProteinMetadataService
    # Setup fallback
    if protein_fallback == 'UNIPROT':
        fallback1 = UniprotProteinMetadataService()
    else:
        raise ValueError(f'Unknown protein fallback annotator type {protein_fallback}')
    # Setup protein metadata cache
    prot_cache_dir = os.path.join(cache_dir, 'protein_cache')
    os.makedirs(prot_cache_dir, exist_ok=True)
    prot_cache = ProteinAnnotationCache(prot_cache_dir)
    # Assemble the final protein metadata service
    protein_metadata_service = ProtCachingFunctionalAnnotator(prot_cache, fallback1)

    # (2) FunctionalAnnotator
    # Setup fallback
    if variant_fallback == 'VEP':
        fallback = VepFunctionalAnnotator(protein_metadata_service)
    else:
        raise ValueError(f'Unknown variant fallback annotator type {variant_fallback}')

    # Setup variant cache
    var_cache_dir = os.path.join(cache_dir, 'variant_cache')
    os.makedirs(var_cache_dir, exist_ok=True)
    var_cache = VariantAnnotationCache(var_cache_dir)

    # Assemble the final functional annotator
    return VarCachingFunctionalAnnotator(var_cache, fallback)


