import logging
import os
import typing
import warnings

import hpotk
# pyright: reportGeneralTypeIssues=false
from google.protobuf.json_format import Parse
from phenopackets import Phenopacket
from tqdm import tqdm

from genophenocorr.model import Cohort
from genophenocorr.model.genome import GRCh37, GRCh38, GenomeBuild
from ._api import FunctionalAnnotator, ProteinMetadataService
from ._audit import AuditReport
from ._patient import CohortCreator
from ._phenopacket import PhenopacketPatientCreator
from ._phenotype import PhenotypeCreator
from ._protein import ProteinAnnotationCache, ProtCachingFunctionalAnnotator
from ._uniprot import UniprotProteinMetadataService
from ._variant import VarCachingFunctionalAnnotator, VariantAnnotationCache
from ._vep import VepFunctionalAnnotator
from ._vv import VVHgvsVariantCoordinateFinder

VALIDATION_POLICIES = {'none', 'lenient', 'strict'}


def configure_caching_cohort_creator(hpo: hpotk.MinimalOntology,
                                     genome_build: str = 'GRCh38.p13',
                                     validation_runner: typing.Optional[hpotk.validate.ValidationRunner] = None,
                                     cache_dir: typing.Optional[str] = None,
                                     variant_fallback: str = 'VEP',
                                     protein_fallback: str = 'UNIPROT') -> CohortCreator[Phenopacket]:
    """
    A convenience function for configuring a caching :class:`genophenocorr.preprocessing.PhenopacketPatientCreator`.

    To create the patient creator, we need hpo-toolkit's representation of HPO. Other options are optional.

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
    pc = configure_caching_patient_creator(hpo,
                                           genome_build,
                                           validation_runner,
                                           cache_dir,
                                           variant_fallback,
                                           protein_fallback)
    if cache_dir is None:
        cache_dir = os.path.join(os.getcwd(), '.genophenocorr_cache')

    build = _configure_build(genome_build)
    phenotype_creator = _setup_phenotype_creator(hpo, validation_runner)
    functional_annotator = _configure_functional_annotator(cache_dir, variant_fallback, protein_fallback)
    hgvs_annotator = VVHgvsVariantCoordinateFinder(build)
    pc = PhenopacketPatientCreator(build, phenotype_creator, functional_annotator, hgvs_annotator)

    return CohortCreator(pc)


def configure_caching_patient_creator(hpo: hpotk.MinimalOntology,
                                      genome_build: str = 'GRCh38.p13',
                                      validation_runner: typing.Optional[hpotk.validate.ValidationRunner] = None,
                                      cache_dir: typing.Optional[str] = None,
                                      variant_fallback: str = 'VEP',
                                      protein_fallback: str = 'UNIPROT') -> PhenopacketPatientCreator:
    """
    A convenience function for configuring a caching :class:`genophenocorr.preprocessing.PhenopacketPatientCreator`.

    To create the patient creator, we need hpo-toolkit's representation of HPO. Other options are optional.

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
    warnings.warn('`configure_caching_patient_creator` was deprecated. '
                  'Use `configure_caching_cohort_creator` instead', DeprecationWarning, stacklevel=2)
    if cache_dir is None:
        cache_dir = os.path.join(os.getcwd(), '.genophenocorr_cache')

    build = _configure_build(genome_build)
    phenotype_creator = _setup_phenotype_creator(hpo, validation_runner)
    functional_annotator = _configure_functional_annotator(cache_dir, variant_fallback, protein_fallback)
    hgvs_annotator = VVHgvsVariantCoordinateFinder(build)
    return PhenopacketPatientCreator(build, phenotype_creator, functional_annotator, hgvs_annotator)


def configure_cohort_creator(hpo: hpotk.MinimalOntology,
                             genome_build: str = 'GRCh38.p13',
                             validation_runner: typing.Optional[hpotk.validate.ValidationRunner] = None,
                             variant_fallback: str = 'VEP',
                             protein_fallback: str = 'UNIPROT') -> CohortCreator[Phenopacket]:
    """
    A convenience function for configuring a non-caching :class:`genophenocorr.preprocessing.PhenopacketPatientCreator`.

    To create the patient creator, we need hpo-toolkit's representation of HPO. Other options are optional

    :param hpo: a HPO instance.
    :param genome_build: name of the genome build to use, choose from `{'GRCh37.p13', 'GRCh38.p13'}`.
    :param validation_runner: an instance of the validation runner.
     if the data should be cached in `.cache` folder in the current working directory.
     In any case, the directory will be created if it does not exist (including non-existing parents).
    :param variant_fallback: the fallback variant annotator to use if we cannot find the annotation locally.
     Choose from ``{'VEP'}`` (just one fallback implementation is available at the moment).
    :param protein_fallback: the fallback protein metadata annotator to use if we cannot find the annotation locally.
     Choose from ``{'UNIPROT'}`` (just one fallback implementation is available at the moment).
    """
    build = _configure_build(genome_build)

    phenotype_creator = _setup_phenotype_creator(hpo, validation_runner)
    protein_metadata_service = _configure_fallback_protein_service(protein_fallback)
    functional_annotator = _configure_fallback_functional(protein_metadata_service, variant_fallback)
    hgvs_annotator = VVHgvsVariantCoordinateFinder(build)
    pc = PhenopacketPatientCreator(build, phenotype_creator, functional_annotator, hgvs_annotator)

    return CohortCreator(pc)


def configure_patient_creator(hpo: hpotk.MinimalOntology,
                              genome_build: str = 'GRCh38.p13',
                              validation_runner: typing.Optional[hpotk.validate.ValidationRunner] = None,
                              variant_fallback: str = 'VEP',
                              protein_fallback: str = 'UNIPROT') -> PhenopacketPatientCreator:
    """
    A convenience function for configuring a non-caching :class:`genophenocorr.preprocessing.PhenopacketPatientCreator`.

    To create the patient creator, we need hpo-toolkit's representation of HPO. Other options are optional

    :param hpo: a HPO instance.
    :param genome_build: name of the genome build to use, choose from `{'GRCh37.p13', 'GRCh38.p13'}`.
    :param validation_runner: an instance of the validation runner.
     if the data should be cached in `.cache` folder in the current working directory.
     In any case, the directory will be created if it does not exist (including non-existing parents).
    :param variant_fallback: the fallback variant annotator to use if we cannot find the annotation locally.
     Choose from ``{'VEP'}`` (just one fallback implementation is available at the moment).
    :param protein_fallback: the fallback protein metadata annotator to use if we cannot find the annotation locally.
     Choose from ``{'UNIPROT'}`` (just one fallback implementation is available at the moment).
    """
    warnings.warn('`configure_patient_creator` was deprecated. '
                  'Use `configure_cohort_creator` instead', DeprecationWarning, stacklevel=2)

    build = _configure_build(genome_build)

    phenotype_creator = _setup_phenotype_creator(hpo, validation_runner)
    protein_metadata_service = _configure_fallback_protein_service(protein_fallback)
    functional_annotator = _configure_fallback_functional(protein_metadata_service, variant_fallback)
    hgvs_annotator = VVHgvsVariantCoordinateFinder(build)
    return PhenopacketPatientCreator(build, phenotype_creator, functional_annotator, hgvs_annotator)


def _configure_build(genome_build: str) -> GenomeBuild:
    if genome_build == 'GRCh38.p13':
        return GRCh38
    elif genome_build == 'GRCh37.p13':
        return GRCh37
    else:
        raise ValueError(f'Unknown build {genome_build}. Choose from [\'GRCh37.p13\', \'GRCh38.p13\']')


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
    protein_fallback = _configure_fallback_protein_service(protein_fallback)
    # Setup protein metadata cache
    prot_cache_dir = os.path.join(cache_dir, 'protein_cache')
    os.makedirs(prot_cache_dir, exist_ok=True)
    prot_cache = ProteinAnnotationCache(prot_cache_dir)
    # Assemble the final protein metadata service
    protein_metadata_service = ProtCachingFunctionalAnnotator(prot_cache, protein_fallback)

    # (2) FunctionalAnnotator
    # Setup fallback
    fallback = _configure_fallback_functional(protein_metadata_service, variant_fallback)

    # Setup variant cache
    var_cache_dir = os.path.join(cache_dir, 'variant_cache')
    os.makedirs(var_cache_dir, exist_ok=True)
    var_cache = VariantAnnotationCache(var_cache_dir)

    # Assemble the final functional annotator
    return VarCachingFunctionalAnnotator(var_cache, fallback)


def _configure_fallback_protein_service(protein_fallback: str) -> ProteinMetadataService:
    if protein_fallback == 'UNIPROT':
        fallback1 = UniprotProteinMetadataService()
    else:
        raise ValueError(f'Unknown protein fallback annotator type {protein_fallback}')
    return fallback1


def _configure_fallback_functional(protein_metadata_service: ProteinMetadataService,
                                   variant_fallback: str) -> FunctionalAnnotator:
    if variant_fallback == 'VEP':
        fallback = VepFunctionalAnnotator(protein_metadata_service)
    else:
        raise ValueError(f'Unknown variant fallback annotator type {variant_fallback}')
    return fallback


def load_phenopacket_folder(pp_directory: str,
                            cohort_creator: CohortCreator[Phenopacket],
                            validation_policy: str = 'none') -> Cohort:
    """
    Creates a Patient object for each phenopacket formatted JSON file in the given directory `pp_directory`.

    :param pp_directory: path to a folder with phenopacket JSON files. An error is raised if the path does not point to
      a directory with at least one phenopacket.
    :param cohort_creator: cohort creator for turning a sequence of phenopacket into a :class:`genophenocorr.model.Cohort`.
    :param validation_policy: a `str` with the validation policy. The value must be one of `{'none', 'lenient', 'strict'}`
    :return: a cohort made of the phenopackets
    """
    # Check inputs before moving a finger.
    hpotk.util.validate_instance(cohort_creator, CohortCreator, 'cohort_creator')
    if not os.path.isdir(pp_directory):
        raise ValueError(f'`{pp_directory}` does not point to a directory')
    if validation_policy.lower() not in VALIDATION_POLICIES:
        raise ValueError(f'{validation_policy} must be one of {VALIDATION_POLICIES}')

    # Load phenopackets
    pps = _load_phenopacket_dir(pp_directory)
    if len(pps) == 0:
        raise ValueError(f"No phenopackets could be parsed from {pp_directory}")

    # Turn phenopackets into a cohort using the cohort creator.
    # Keep track of the progress by wrapping the list of phenopackets
    # with TQDM 😎
    cohort_iter = tqdm(pps, desc='Patients Created')
    cohort_report = cohort_creator.process(cohort_iter)

    logger = logging.getLogger('genophenocorr.preprocessing')
    _summarize_validation(validation_policy, cohort_report, logger)
    if validation_policy == 'none':
        # No validation
        return cohort_report.outcome
    elif validation_policy == 'lenient':
        if cohort_report.has_errors():
            raise ValueError(f'Cannot load cohort due to errors')
    elif validation_policy == 'strict':
        if cohort_report.has_warnings_or_errors():
            raise ValueError(f'Cannot load cohort due to warnings or errors')
    else:
        # Bug, please report to the developers.
        raise ValueError(f'Unexpected policy {validation_policy}')

    # create cohort from patients
    return cohort_report.outcome


def _summarize_validation(policy: str,
                          audit_report: AuditReport,
                          logger):
    logger.info(f'Validated under {policy} policy')
    errors = tuple(audit_report.errors())
    warnings = tuple(audit_report.warnings())
    logger.info(f'Found {len(errors)} error(s) and {len(warnings)} warnings')
    if len(errors) > 0:
        logger.info('Errors:')
        for e in errors:
            logger.info(f' - {e.message}. {e.solution}.')
    if len(warnings) > 0:
        logger.info('Warnings:')
        for w in warnings:
            logger.info(f' - {w.message}. {w.solution}.')


def _load_phenopacket_dir(pp_dir: str) -> typing.Sequence[Phenopacket]:
    patients = []
    for patient_file in os.listdir(pp_dir):
        if patient_file.endswith('.json'):
            phenopacket_path = os.path.join(pp_dir, patient_file)
            pp = load_phenopacket(phenopacket_path)
            patients.append(pp)
    return patients


def load_phenopacket(phenopacket_path: str) -> Phenopacket:
    """
    Load phenopacket JSON file.

    :param phenopacket_path: a `str` pointing to phenopacket JSON file.
    """
    with open(phenopacket_path) as f:
        return Parse(f.read(), Phenopacket())
