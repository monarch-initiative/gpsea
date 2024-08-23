import os
import sys
import typing
import warnings

import hpotk

# pyright: reportGeneralTypeIssues=false
from google.protobuf.json_format import Parse
from phenopackets import Phenopacket
from tqdm import tqdm

from gpsea.config import get_cache_dir_path
from gpsea.model import Cohort
from gpsea.model.genome import GRCh37, GRCh38, GenomeBuild
from ._api import FunctionalAnnotator, PreprocessingValidationResult
from ._generic import DefaultImpreciseSvFunctionalAnnotator
from ._patient import CohortCreator
from ._phenopacket import PhenopacketPatientCreator
from ._phenotype import PhenotypeCreator
from ._protein import (
    ProteinMetadataService,
    ProtCachingMetadataService,
    ProteinAnnotationCache,
)
from ._uniprot import UniprotProteinMetadataService
from ._variant import VarCachingFunctionalAnnotator, VariantAnnotationCache
from ._vep import VepFunctionalAnnotator
from ._vv import VVHgvsVariantCoordinateFinder, VVMultiCoordinateService

VALIDATION_POLICIES = {"none", "lenient", "strict"}


def configure_caching_cohort_creator(
    hpo: hpotk.MinimalOntology,
    genome_build: str = "GRCh38.p13",
    validation_runner: typing.Optional[hpotk.validate.ValidationRunner] = None,
    cache_dir: typing.Optional[str] = None,
    variant_fallback: str = "VEP",
    timeout: float = 30.0,
) -> CohortCreator[Phenopacket]:
    """
    A convenience function for configuring a caching :class:`~gpsea.preprocessing.PhenopacketPatientCreator`.

    To create the patient creator, we need hpo-toolkit's representation of HPO. Other options are optional.

    :param hpo: a HPO instance.
    :param genome_build: name of the genome build to use, choose from `{'GRCh37.p13', 'GRCh38.p13'}`.
    :param validation_runner: an instance of the validation runner.
    :param cache_dir: path to the folder where we will cache the results fetched from the remote APIs or `None`
        if the cache location should be determined as described in :func:`~gpsea.config.get_cache_dir_path`.
        In any case, the directory will be created if it does not exist (including non-existing parents).
    :param variant_fallback: the fallback variant annotator to use if we cannot find the annotation locally.
     Choose from ``{'VEP'}`` (just one fallback implementation is available at the moment).
    :param timeout: timeout in seconds for the REST APIs
    """
    cache_dir = _configure_cache_dir(cache_dir)

    build = _configure_build(genome_build)
    phenotype_creator = _setup_phenotype_creator(hpo, validation_runner)
    functional_annotator = _configure_functional_annotator(
        cache_dir, variant_fallback, timeout
    )
    imprecise_sv_functional_annotator = _configure_imprecise_sv_annotator(
        build, cache_dir, timeout
    )
    hgvs_annotator = VVHgvsVariantCoordinateFinder(build)
    pc = PhenopacketPatientCreator(
        build=build,
        phenotype_creator=phenotype_creator,
        functional_annotator=functional_annotator,
        imprecise_sv_functional_annotator=imprecise_sv_functional_annotator,
        hgvs_coordinate_finder=hgvs_annotator,
    )

    return CohortCreator(pc)


def configure_caching_patient_creator(
    hpo: hpotk.MinimalOntology,
    genome_build: str = "GRCh38.p13",
    validation_runner: typing.Optional[hpotk.validate.ValidationRunner] = None,
    cache_dir: typing.Optional[str] = None,
    variant_fallback: str = "VEP",
    timeout: float = 30.0,
) -> PhenopacketPatientCreator:
    """
    A convenience function for configuring a caching :class:`~gpsea.preprocessing.PhenopacketPatientCreator`.

    To create the patient creator, we need hpo-toolkit's representation of HPO. Other options are optional.

    :param hpo: a HPO instance.
    :param genome_build: name of the genome build to use, choose from `{'GRCh37.p13', 'GRCh38.p13'}`.
    :param validation_runner: an instance of the validation runner.
    :param cache_dir: path to the folder where we will cache the results fetched from the remote APIs or `None`
        if the cache location should be determined as described in :func:`~gpsea.config.get_cache_dir_path`.
        In any case, the directory will be created if it does not exist (including non-existing parents).
    :param variant_fallback: the fallback variant annotator to use if we cannot find the annotation locally.
     Choose from ``{'VEP'}`` (just one fallback implementation is available at the moment).
    :param timeout: timeout in seconds for the REST APIs
    """
    warnings.warn(
        "`configure_caching_patient_creator` was deprecated. "
        "Use `configure_caching_cohort_creator` instead",
        DeprecationWarning,
        stacklevel=2,
    )
    cache_dir = _configure_cache_dir(cache_dir)

    build = _configure_build(genome_build)
    phenotype_creator = _setup_phenotype_creator(hpo, validation_runner)
    functional_annotator = _configure_functional_annotator(
        cache_dir, variant_fallback, timeout
    )
    imprecise_sv_functional_annotator = _configure_imprecise_sv_annotator(
        build, timeout
    )
    hgvs_annotator = VVHgvsVariantCoordinateFinder(build)
    return PhenopacketPatientCreator(
        build=build,
        phenotype_creator=phenotype_creator,
        functional_annotator=functional_annotator,
        imprecise_sv_functional_annotator=imprecise_sv_functional_annotator,
        hgvs_coordinate_finder=hgvs_annotator,
    )


def configure_cohort_creator(
    hpo: hpotk.MinimalOntology,
    genome_build: str = "GRCh38.p13",
    validation_runner: typing.Optional[hpotk.validate.ValidationRunner] = None,
    variant_fallback: str = "VEP",
    timeout: float = 30.0,
) -> CohortCreator[Phenopacket]:
    """
    A convenience function for configuring a non-caching :class:`~gpsea.preprocessing.PhenopacketPatientCreator`.

    To create the patient creator, we need hpo-toolkit's representation of HPO. Other options are optional

    :param hpo: a HPO instance.
    :param genome_build: name of the genome build to use, choose from `{'GRCh37.p13', 'GRCh38.p13'}`.
    :param validation_runner: an instance of the validation runner.
     if the data should be cached in `.cache` folder in the current working directory.
     In any case, the directory will be created if it does not exist (including non-existing parents).
    :param variant_fallback: the fallback variant annotator to use if we cannot find the annotation locally.
     Choose from ``{'VEP'}`` (just one fallback implementation is available at the moment).
    :param timeout: timeout in seconds for the VEP API
    """
    build = _configure_build(genome_build)

    phenotype_creator = _setup_phenotype_creator(hpo, validation_runner)
    functional_annotator = _configure_fallback_functional(variant_fallback, timeout)
    imprecise_sv_functional_annotator = _configure_imprecise_sv_annotator(
        build, timeout
    )
    hgvs_annotator = VVHgvsVariantCoordinateFinder(build)
    pc = PhenopacketPatientCreator(
        build=build,
        phenotype_creator=phenotype_creator,
        functional_annotator=functional_annotator,
        imprecise_sv_functional_annotator=imprecise_sv_functional_annotator,
        hgvs_coordinate_finder=hgvs_annotator,
    )

    return CohortCreator(pc)


def configure_patient_creator(
    hpo: hpotk.MinimalOntology,
    genome_build: str = "GRCh38.p13",
    validation_runner: typing.Optional[hpotk.validate.ValidationRunner] = None,
    variant_fallback: str = "VEP",
    validation: str = "lenient",
    timeout: float = 30.0,
) -> PhenopacketPatientCreator:  # Rename to something more understandable by user
    """
    A convenience function for configuring a non-caching :class:`~gpsea.preprocessing.PhenopacketPatientCreator`.

    To create the patient creator, we need hpo-toolkit's representation of HPO. Other options are optional

    :param hpo: a HPO instance.
    :param genome_build: name of the genome build to use, choose from `{'GRCh37.p13', 'GRCh38.p13'}`.
    :param validation_runner: an instance of the validation runner.
     if the data should be cached in `.gpsea_cache` folder in the current working directory.
     In any case, the directory will be created if it does not exist (including non-existing parents).
    :param variant_fallback: the fallback variant annotator to use if we cannot find the annotation locally.
     Choose from ``{'VEP'}`` (just one fallback implementation is available at the moment).
    :param timeout: timeout in seconds for the REST APIs
    """
    warnings.warn(
        "`configure_patient_creator` was deprecated. "
        "Use `configure_cohort_creator` instead",
        DeprecationWarning,
        stacklevel=2,
    )

    build = _configure_build(genome_build)

    phenotype_creator = _setup_phenotype_creator(hpo, validation_runner)
    functional_annotator = _configure_fallback_functional(variant_fallback, timeout)
    imprecise_sv_functional_annotator = _configure_imprecise_sv_annotator(
        build, timeout
    )
    hgvs_annotator = VVHgvsVariantCoordinateFinder(build)
    return PhenopacketPatientCreator(
        build=build,
        phenotype_creator=phenotype_creator,
        functional_annotator=functional_annotator,
        imprecise_sv_functional_annotator=imprecise_sv_functional_annotator,
        hgvs_coordinate_finder=hgvs_annotator,
    )


def configure_protein_metadata_service(
    cache_dir: typing.Optional[str] = None,
    timeout: float = 30.0,
) -> ProteinMetadataService:
    """
    Configure default protein metadata service.

    The service will cache the responses in `cache_dir` and reach out to UNIPROT API for cache misses.

    :param cache_dir: path to the folder where we will cache the results fetched from the remote APIs or `None`
     if the data should be cached in `.gpsea_cache` folder in the current working directory.
     In any case, the directory will be created if it does not exist (including any non-existing parents).
    :param timeout: timeout in seconds for the REST APIs.
    """
    cache_dir = _configure_cache_dir(cache_dir)

    protein_cache_dir = os.path.join(cache_dir, "protein_cache")
    os.makedirs(protein_cache_dir, exist_ok=True)

    cache = ProteinAnnotationCache(protein_cache_dir)
    fallback = UniprotProteinMetadataService(timeout=timeout)
    return ProtCachingMetadataService(cache=cache, fallback=fallback)


def _configure_cache_dir(
    cache_dir: typing.Optional[str] = None,
) -> str:
    cache_dir = get_cache_dir_path(cache_dir)
    os.makedirs(cache_dir, exist_ok=True)

    return str(cache_dir)


def _configure_build(genome_build: str) -> GenomeBuild:
    if genome_build == "GRCh38.p13":
        return GRCh38
    elif genome_build == "GRCh37.p13":
        return GRCh37
    else:
        raise ValueError(
            f"Unknown build {genome_build}. Choose from ['GRCh37.p13', 'GRCh38.p13']"
        )


def _setup_phenotype_creator(
    hpo: hpotk.MinimalOntology,
    validator: typing.Optional[hpotk.validate.ValidationRunner],
) -> PhenotypeCreator:
    if validator is None:
        # This will be the default validator
        validator = hpotk.validate.ValidationRunner(
            (
                hpotk.validate.ObsoleteTermIdsValidator(hpo),
                hpotk.validate.AnnotationPropagationValidator(hpo),
                hpotk.validate.PhenotypicAbnormalityValidator(hpo),
            )
        )
    else:
        validator = hpotk.util.validate_instance(
            validator, hpotk.validate.ValidationRunner, "validator"
        )
    return PhenotypeCreator(hpo, validator)


def _configure_functional_annotator(
    cache_dir: str,
    variant_fallback: str,
    timeout: float,
) -> FunctionalAnnotator:

    # (2) FunctionalAnnotator
    # Setup fallback
    fallback = _configure_fallback_functional(variant_fallback, timeout)

    # Setup variant cache
    var_cache_dir = os.path.join(cache_dir, "variant_cache")
    os.makedirs(var_cache_dir, exist_ok=True)
    var_cache = VariantAnnotationCache(var_cache_dir)

    # Assemble the final functional annotator
    return VarCachingFunctionalAnnotator(var_cache, fallback)


def _configure_fallback_functional(
    variant_fallback: str,
    timeout: float,
) -> FunctionalAnnotator:
    if variant_fallback == "VEP":
        fallback = VepFunctionalAnnotator(timeout=timeout)
    else:
        raise ValueError(f"Unknown variant fallback annotator type {variant_fallback}")
    return fallback


def _configure_imprecise_sv_annotator(
    genome_build: GenomeBuild,
    cache_dir: str,
    timeout: float,
):
    # Setup cache for SVs
    sv_cache_dir = os.path.join(cache_dir, "sv_cache")
    # TODO: implement the cache.
    # os.makedirs(sv_cache_dir, exist_ok=True)
    # var_cache = VariantAnnotationCache(sv_cache_dir)

    return DefaultImpreciseSvFunctionalAnnotator(
        gene_coordinate_service=VVMultiCoordinateService(
            genome_build=genome_build,
            timeout=timeout,
        )
    )


def load_phenopacket_folder(
    pp_directory: str,
    cohort_creator: CohortCreator[Phenopacket],
    validation_policy: typing.Literal["none", "lenient", "strict"] = "none",
) -> typing.Tuple[Cohort, PreprocessingValidationResult]:
    """
    Load phenopacket JSON files from a directory, validate the patient data, and assemble the patients into a cohort.

    :param pp_directory: path to a folder with phenopacket JSON files. An error is raised if the path does not point to
      a directory with at least one phenopacket.
    :param cohort_creator: cohort creator for turning a sequence of phenopacket
      into a :class:`~gpsea.model.Cohort`.
    :param validation_policy: a `str` with the validation policy.
      The value must be one of `{'none', 'lenient', 'strict'}`
    :return: a tuple with the cohort and the preprocessing validation result.
    """
    # Load phenopackets
    phenopackets = _load_phenopacket_dir(pp_directory)

    # Map to patients
    return load_phenopackets(
        phenopackets=phenopackets,
        cohort_creator=cohort_creator,
        validation_policy=validation_policy,
    )


def load_phenopackets(
    phenopackets: typing.Iterator[Phenopacket],
    cohort_creator: CohortCreator[Phenopacket],
    validation_policy: typing.Literal["none", "lenient", "strict"] = "none",
) -> typing.Tuple[Cohort, PreprocessingValidationResult]:
    """
    Map phenopacket JSON file into patient, validate the patient data, and assemble the patients into a cohort.

    :param pp_directory: path to a folder with phenopacket JSON files. An error is raised if the path does not point to
      a directory with at least one phenopacket.
    :param cohort_creator: cohort creator for turning a sequence of phenopacket
      into a :class:`~gpsea.model.Cohort`.
    :param validation_policy: a `str` with the validation policy.
      The value must be one of `{'none', 'lenient', 'strict'}`
    :return: a tuple with the cohort and the preprocessing validation result.
    """
    # Check inputs before doing anything
    hpotk.util.validate_instance(cohort_creator, CohortCreator, "cohort_creator")
    if validation_policy.lower() not in VALIDATION_POLICIES:
        raise ValueError(f"{validation_policy} must be one of {VALIDATION_POLICIES}")

    # Turn phenopackets into a cohort using the cohort creator.
    # Keep track of the progress by wrapping the list of phenopackets
    # with TQDM 😎
    cohort_iter = tqdm(phenopackets, desc="Patients Created", file=sys.stdout)
    notepad = cohort_creator.prepare_notepad("Phenopackets")
    cohort = cohort_creator.process(cohort_iter, notepad)

    validation_result = PreprocessingValidationResult(
        policy=validation_policy,
        notepad=notepad,
    )

    return cohort, validation_result


def _load_phenopacket_dir(
    pp_dir: str,
) -> typing.Iterator[Phenopacket]:
    fpath_pp_abs = os.path.abspath(pp_dir)
    if not os.path.isdir(fpath_pp_abs):
        raise ValueError(f"`{fpath_pp_abs}` does not point to a directory")

    for patient_file in os.listdir(pp_dir):
        if patient_file.endswith(".json"):
            phenopacket_path = os.path.join(pp_dir, patient_file)
            yield _load_phenopacket(phenopacket_path)


def _load_phenopacket(phenopacket_path: str) -> Phenopacket:
    """
    Load phenopacket JSON file.

    :param phenopacket_path: a `str` pointing to phenopacket JSON file.
    """
    with open(phenopacket_path) as f:
        return Parse(f.read(), Phenopacket())
