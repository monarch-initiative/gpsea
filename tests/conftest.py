import json
import os
import typing

import hpotk
import pytest

from hpotk.validate import (
    AnnotationPropagationValidator,
    ObsoleteTermIdsValidator,
    PhenotypicAbnormalityValidator,
    ValidationRunner,
)

from gpsea.analysis.clf import GenotypeClassifier, biallelic_classifier, allele_count
from gpsea.analysis.clf import PhenotypeClassifier, HpoClassifier
from gpsea.analysis.predicate import variant_effect
from gpsea.io import GpseaJSONDecoder
from gpsea.model import (
    Cohort,
    Disease,
    Genotype,
    Genotypes,
    Patient,
    Phenotype,
    ProteinMetadata,
    SampleLabels,
    Sex,
    TranscriptAnnotation,
    TranscriptCoordinates,
    Variant,
    VariantCoordinates,
    VariantEffect,
    VariantInfo,
)
from gpsea.model.genome import GRCh38, GenomicRegion, Region, Strand, GenomeBuild


def pytest_addoption(parser):
    parser.addoption(
        "--runonline", action="store_true", default=False, help="run online tests"
    )


def pytest_configure(config):
    config.addinivalue_line(
        "markers", "online: mark test that require internet access to run"
    )


def pytest_collection_modifyitems(config, items):
    if config.getoption("--runonline"):
        # --runonline given in cli: do not skip online tests
        return
    skip_online = pytest.mark.skip(reason="need --runonline option to run")
    for item in items:
        if "online" in item.keywords:
            item.add_marker(skip_online)


@pytest.fixture(scope="session")
def fpath_project_dir(fpath_test_dir: str) -> str:
    """
    Path to project folder, where `pyproject.toml`, `README.md`,
    as well as `src` and `tests` folders are located.
    """
    return os.path.dirname(fpath_test_dir)


@pytest.fixture(scope="session")
def fpath_test_dir() -> str:
    """
    Path to `tests` folder.
    """
    return os.path.dirname(os.path.abspath(__file__))


@pytest.fixture(scope="session")
def fpath_docs_dir(fpath_project_dir: str) -> str:
    """
    Path to `docs` folder.
    """
    return os.path.join(fpath_project_dir, "docs")


@pytest.fixture(scope="session")
def fpath_cohort_data_dir(fpath_docs_dir: str) -> str:
    return os.path.join(fpath_docs_dir, "cohort-data")


@pytest.fixture(scope="session")
def fpath_test_data_dir(fpath_test_dir: str) -> str:
    return os.path.join(fpath_test_dir, "test_data")


@pytest.fixture(scope="session")
def fpath_phenopacket_dir(fpath_test_data_dir: str) -> str:
    return os.path.join(fpath_test_data_dir, "phenopackets")


@pytest.fixture(scope="session")
def fpath_toy_hpo(fpath_test_data_dir: str) -> str:
    return os.path.join(fpath_test_data_dir, "hp.toy.json")


@pytest.fixture(scope="session")
def toy_hpo(fpath_toy_hpo: str) -> hpotk.MinimalOntology:
    return hpotk.load_minimal_ontology(fpath_toy_hpo)


@pytest.fixture(scope="session")
def hpo(fpath_test_data_dir: str) -> hpotk.MinimalOntology:
    fpath_hpo = os.path.join(fpath_test_data_dir, "hp.v2024-04-26.json.gz")
    return hpotk.load_minimal_ontology(fpath_hpo)


@pytest.fixture(scope="session")
def validation_runner(hpo: hpotk.MinimalOntology) -> ValidationRunner:
    validators = (
        ObsoleteTermIdsValidator(hpo),
        AnnotationPropagationValidator(hpo),
        PhenotypicAbnormalityValidator(hpo),
    )
    return ValidationRunner(validators)


def make_region(contig: str, start: int, end: int) -> GenomicRegion:
    a_contig = GRCh38.contig_by_name(contig)
    assert a_contig is not None
    return GenomicRegion(a_contig, start, end, Strand.POSITIVE)


@pytest.fixture(scope="session")
def fpath_suox_cohort(
    fpath_test_data_dir: str,
) -> str:
    return os.path.join(fpath_test_data_dir, "SUOX.json")


@pytest.fixture(scope="session")
def suox_cohort(
    fpath_suox_cohort: str,
) -> Cohort:
    with open(fpath_suox_cohort) as fh:
        return json.load(fh, cls=GpseaJSONDecoder)


@pytest.fixture(scope="session")
def suox_mane_tx_id() -> str:
    return "NM_001032386.2"


@pytest.fixture(scope="session")
def suox_gt_clf(
    suox_mane_tx_id: str,
) -> GenotypeClassifier:
    return allele_count(
        counts=((0,), (1,)),
        target=variant_effect(
            effect=VariantEffect.MISSENSE_VARIANT, tx_id=suox_mane_tx_id
        ),
    )


@pytest.fixture(scope="session")
def suox_pheno_clfs(
    hpo: hpotk.MinimalOntology,
) -> typing.Sequence[PhenotypeClassifier[hpotk.TermId]]:
    """
    Get classifiers to test for presence of 5 HPO terms:
    Seizure, Ectopia lentis, Sulfocysteinuria, Neurodevelopmental delay, and Hypertonia.

    Note, these are just a *SUBSET* of all phenotypes that can be tested for in the *SUOX* cohort.
    """
    return (
        HpoClassifier(
            hpo=hpo,
            query=hpotk.TermId.from_curie("HP:0001250"),  # Seizure
        ),
        HpoClassifier(
            hpo=hpo,
            query=hpotk.TermId.from_curie("HP:0001083"),  # Ectopia lentis
        ),
        HpoClassifier(
            hpo=hpo,
            query=hpotk.TermId.from_curie("HP:0032350"),  # Sulfocysteinuria
        ),
        HpoClassifier(
            hpo=hpo,
            query=hpotk.TermId.from_curie("HP:0012758"),  # Neurodevelopmental delay
        ),
        HpoClassifier(
            hpo=hpo,
            query=hpotk.TermId.from_curie("HP:0001276"),  # Hypertonia
        ),
    )


@pytest.fixture(scope="session")
def fpath_suox_tx_coordinates(fpath_test_data_dir: str) -> str:
    suox_mane_tx_id = "NM_001032386.2"
    return os.path.join(fpath_test_data_dir, f"SUOX-{suox_mane_tx_id}.json")


@pytest.fixture(scope="session")
def suox_mane_tx_coordinates(
    fpath_suox_tx_coordinates: str,
) -> TranscriptCoordinates:
    with open(fpath_suox_tx_coordinates) as fh:
        return json.load(fh, cls=GpseaJSONDecoder)


@pytest.fixture(scope="session")
def fpath_suox_protein_metadata(fpath_test_data_dir: str) -> str:
    suox_mane_tx_protein_id = "NP_001027558.1"
    return os.path.join(fpath_test_data_dir, f"SUOX-{suox_mane_tx_protein_id}.json")


@pytest.fixture(scope="session")
def suox_protein_metadata(
    fpath_suox_protein_metadata: str,
) -> ProteinMetadata:
    with open(fpath_suox_protein_metadata) as fh:
        return json.load(fh, cls=GpseaJSONDecoder)


@pytest.fixture(scope="session")
def fpath_cyp21a2_cohort(
    fpath_test_data_dir: str,
) -> str:
    # Generated from Phenopacket Store `0.1.20`.
    return os.path.join(fpath_test_data_dir, "CYP21A2.0.1.20.json")


@pytest.fixture(scope="session")
def cyp21a2_cohort(
    fpath_cyp21a2_cohort: str,
) -> Cohort:
    with open(fpath_cyp21a2_cohort) as fh:
        return json.load(fh, cls=GpseaJSONDecoder)


@pytest.fixture(scope="session")
def cyp21a2_mane_tx_id() -> str:
    return "NM_000500.9"


@pytest.fixture(scope="session")
def cyp21a2_gt_clf(
    cyp21a2_mane_tx_id: str,
) -> GenotypeClassifier:
    return biallelic_classifier(
        a_predicate=variant_effect(
            effect=VariantEffect.MISSENSE_VARIANT,
            tx_id=cyp21a2_mane_tx_id,
        ),
        a_label="Missense",
        b_label="Other",
    )


@pytest.fixture(scope="session")
def cyp21a2_testosterone_label() -> str:
    return "LOINC:2986-8"


@pytest.fixture(scope="session")
def toy_cohort(
    test_phenotypes: typing.Mapping[str, Phenotype],
    test_diseases: typing.Mapping[str, Disease],
) -> Cohort:
    prot_id = "NP_037407.4"

    dup = Variant(
        variant_info=VariantInfo(
            variant_coordinates=VariantCoordinates(
                make_region("16", 89279849, 89279850),
                ref="G",
                alt="GC",
                change_length=1,
            )
        ),
        tx_annotations=[
            TranscriptAnnotation(
                "ANKRD11",
                "NM_013275.6",
                "NM_013275.6:c.6691dup",
                False,
                [VariantEffect.FRAMESHIFT_VARIANT],
                [9],
                prot_id,
                "NP_001243112.1:p.Ala2231GlyfsTer29",
                Region(2230, 2231),
            )
        ],
        genotypes=Genotypes.from_mapping(
            {SampleLabels("HetSingleVar"): Genotype.HETEROZYGOUS}
        ),
    )
    indel = Variant(
        variant_info=VariantInfo(
            variant_coordinates=VariantCoordinates(
                make_region("16", 89284600, 89284602),
                ref="GG",
                alt="A",
                change_length=-1,
            ),
        ),
        tx_annotations=[
            TranscriptAnnotation(
                "ANKRD11",
                "NM_013275.6",
                "NM_013275.6:c.1940_1941delinsT",
                False,
                [VariantEffect.FRAMESHIFT_VARIANT],
                [9],
                prot_id,
                "NP_001243111.1:p.Ser647LeufsTer6",
                Region(646, 647),
            )
        ],
        genotypes=Genotypes.from_mapping(
            {SampleLabels("HetDoubleVar1"): Genotype.HETEROZYGOUS}
        ),
    )
    snv_stop_gain = Variant(
        variant_info=VariantInfo(
            variant_coordinates=VariantCoordinates(
                make_region("16", 89280751, 89280752),
                ref="G",
                alt="T",
                change_length=0,
            )
        ),
        tx_annotations=[
            TranscriptAnnotation(
                "ANKRD11",
                "NM_013275.6",
                "NM_013275.6:c.5790C>A",
                False,
                [VariantEffect.STOP_GAINED],
                [9],
                prot_id,
                "NP_037407.4:p.Tyr1930Ter",
                Region(1929, 1930),
            )
        ],
        genotypes=Genotypes.from_mapping(
            {SampleLabels("HetDoubleVar1"): Genotype.HETEROZYGOUS}
        ),
    )
    snv_missense = Variant(
        variant_info=VariantInfo(
            variant_coordinates=VariantCoordinates(
                make_region("16", 89275127, 89275128),
                ref="G",
                alt="A",
                change_length=0,
            )
        ),
        tx_annotations=[
            TranscriptAnnotation(
                "ANKRD11",
                "NM_013275.6",
                "NM_013275.6:c.7534C>T",
                False,
                [VariantEffect.MISSENSE_VARIANT],
                [10],
                prot_id,
                "NP_037407.4:p.Arg2512Trp",
                Region(2511, 2512),
            )
        ],
        genotypes=Genotypes.from_mapping(
            {SampleLabels("HetDoubleVar2"): Genotype.HETEROZYGOUS}
        ),
    )
    del_frameshift = Variant(
        variant_info=VariantInfo(
            variant_coordinates=VariantCoordinates(
                make_region("16", 89279707, 89279725),
                ref="AGTGTTCGGGGCGGGGCC",
                alt="A",
                change_length=-17,
            )
        ),
        tx_annotations=[
            TranscriptAnnotation(
                "ANKRD11",
                "NM_013275.6",
                "NM_013275.6:c.6817_6833del",
                False,
                [VariantEffect.FRAMESHIFT_VARIANT],
                [9],
                prot_id,
                "NP_037407.4:p.Gly2273CysfsTer17",
                Region(2272, 2278),
            )
        ],
        genotypes=Genotypes.from_mapping(
            {SampleLabels("HetDoubleVar2"): Genotype.HETEROZYGOUS}
        ),
    )
    del_small = Variant(
        variant_info=VariantInfo(
            variant_coordinates=VariantCoordinates(
                make_region("16", 89279457, 89279459),
                ref="TG",
                alt="T",
                change_length=-1,
            )
        ),
        tx_annotations=[
            TranscriptAnnotation(
                "ANKRD11",
                "NM_013275.6",
                "NM_013275.6:c.7083del",
                False,
                [VariantEffect.FRAMESHIFT_VARIANT],
                [9],
                prot_id,
                "NP_037407.4:p.Thr2362ProfsTer39",
                Region(2360, 2362),
            )
        ],
        genotypes=Genotypes.from_mapping(
            {SampleLabels("HomoVar"): Genotype.HOMOZYGOUS_ALTERNATE}
        ),
    )
    del_large = Variant(
        variant_info=VariantInfo(
            variant_coordinates=VariantCoordinates(
                make_region("16", 89_190_070, 89_439_815),
                ref="N",
                alt="<DEL>",
                change_length=-249_745,
            )
        ),
        tx_annotations=[
            TranscriptAnnotation(
                "ANKRD11",
                "NM_013275.6",
                None,
                False,
                [
                    VariantEffect.STOP_LOST,
                    VariantEffect.FEATURE_TRUNCATION,
                    VariantEffect.CODING_SEQUENCE_VARIANT,
                    VariantEffect.FIVE_PRIME_UTR_VARIANT,
                    VariantEffect.THREE_PRIME_UTR_VARIANT,
                    VariantEffect.INTRON_VARIANT,
                ],
                [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13],
                prot_id,
                None,
                None,
            )
        ],
        genotypes=Genotypes.from_mapping(
            {SampleLabels("LargeCNV"): Genotype.HETEROZYGOUS}
        ),
    )

    patients = (
        Patient.from_raw_parts(
            SampleLabels("HetSingleVar"),
            sex=Sex.UNKNOWN_SEX,
            age=None,
            vital_status=None,
            phenotypes=(
                test_phenotypes["arachnodactyly_T"],
                test_phenotypes["spasticity_F"],
                test_phenotypes["focal_clonic_seizure_T"],
            ),
            measurements=(),
            variants=(dup,),
            diseases=(test_diseases["KBG_T"],),
        ),
        Patient.from_raw_parts(
            SampleLabels("HetDoubleVar1"),
            sex=Sex.UNKNOWN_SEX,
            age=None,
            vital_status=None,
            phenotypes=(
                test_phenotypes["arachnodactyly_T"],
                test_phenotypes["seizure_T"],
                test_phenotypes["spasticity_T"],
            ),
            measurements=(),
            variants=(indel, snv_stop_gain),
            diseases=(test_diseases["KBG_T"],),
        ),
        Patient.from_raw_parts(
            SampleLabels("HetDoubleVar2"),
            sex=Sex.UNKNOWN_SEX,
            age=None,
            vital_status=None,
            phenotypes=(
                test_phenotypes["arachnodactyly_F"],
                test_phenotypes["spasticity_T"],
                test_phenotypes["seizure_T"],
            ),
            measurements=(),
            variants=(snv_missense, del_frameshift),
            diseases=(test_diseases["KBG_T"],),
        ),
        Patient.from_raw_parts(
            SampleLabels("HomoVar"),
            sex=Sex.UNKNOWN_SEX,
            age=None,
            vital_status=None,
            phenotypes=(
                test_phenotypes["arachnodactyly_T"],
                test_phenotypes["spasticity_T"],
                test_phenotypes["seizure_T"],
            ),
            measurements=(),
            variants=(del_small,),
            diseases=(),
        ),
        Patient.from_raw_parts(
            SampleLabels("LargeCNV"),
            sex=Sex.UNKNOWN_SEX,
            age=None,
            vital_status=None,
            phenotypes=(
                test_phenotypes["arachnodactyly_T"],
                test_phenotypes["spasticity_T"],
                test_phenotypes["seizure_F"],
            ),
            measurements=(),
            variants=(del_large,),
            diseases=(),
        ),
    )

    return Cohort.from_patients(patients)


@pytest.fixture(scope="session")
def test_phenotypes() -> typing.Mapping[str, Phenotype]:
    return {
        "arachnodactyly_T": Phenotype.from_raw_parts(
            hpotk.TermId.from_curie("HP:0001166"), True
        ),
        "seizure_T": Phenotype.from_raw_parts(
            hpotk.TermId.from_curie("HP:0001250"), True
        ),
        "focal_clonic_seizure_T": Phenotype.from_raw_parts(
            hpotk.TermId.from_curie("HP:0002266"),
            True,
        ),
        "spasticity_T": Phenotype.from_raw_parts(
            hpotk.TermId.from_curie("HP:0001257"), True
        ),
        "arachnodactyly_F": Phenotype.from_raw_parts(
            hpotk.TermId.from_curie("HP:0001166"), False
        ),
        "seizure_F": Phenotype.from_raw_parts(
            hpotk.TermId.from_curie("HP:0001250"), False
        ),
        "spasticity_F": Phenotype.from_raw_parts(
            hpotk.TermId.from_curie("HP:0001257"), False
        ),
        "focal_clonic_seizure_F": Phenotype.from_raw_parts(
            hpotk.TermId.from_curie("HP:0002266"),
            False,
        ),
    }


@pytest.fixture(scope="session")
def test_diseases() -> typing.Mapping[str, Disease]:
    return {
        "KBG_T": Disease.from_raw_parts(
            hpotk.TermId.from_curie("OMIM:148050"), "KBG syndrome", True
        ),
        "KBG_F": Disease.from_raw_parts(
            hpotk.TermId.from_curie("OMIM:148050"), "KBG syndrome", False
        ),
    }


@pytest.fixture(scope="session")
def genome_build() -> GenomeBuild:
    return GRCh38
