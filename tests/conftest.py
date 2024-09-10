import json
import os
import typing

import hpotk
import numpy as np
import pandas as pd
import pytest

from gpsea.analysis.mtc_filter import PhenotypeMtcResult
from gpsea.analysis.pcats import HpoTermAnalysisResult
from gpsea.analysis.predicate.genotype import GenotypePolyPredicate, VariantPredicates, boolean_predicate
from gpsea.analysis.predicate.phenotype import PhenotypePolyPredicate, HpoPredicate
from gpsea.io import GpseaJSONDecoder
from gpsea.model import *
from gpsea.model.genome import GRCh38, GenomicRegion, Region, Strand, GenomeBuild
from ._protein_test_service import ProteinTestMetadataService


def pytest_addoption(parser):
    parser.addoption(
        "--runonline", action="store_true", default=False, help="run online tests"
    )


def pytest_configure(config):
    config.addinivalue_line("markers", "online: mark test that require internet access to run")


def pytest_collection_modifyitems(config, items):
    if config.getoption("--runonline"):
        # --runonline given in cli: do not skip online tests
        return
    skip_online = pytest.mark.skip(reason="need --runonline option to run")
    for item in items:
        if "online" in item.keywords:
            item.add_marker(skip_online)


@pytest.fixture(scope='session')
def fpath_project_dir(fpath_test_dir: str) -> str:
    """
    Path to project folder, where `pyproject.toml`, `README.md`,
    as well as `src` and `tests` folders are located.
    """
    return os.path.dirname(fpath_test_dir)


@pytest.fixture(scope='session')
def fpath_test_dir() -> str:
    """
    Path to `tests` folder
    """
    return os.path.dirname(os.path.abspath(__file__))

@pytest.fixture(scope='session')
def fpath_test_data_dir(fpath_test_dir: str) -> str:
    return os.path.join(fpath_test_dir, 'test_data')


@pytest.fixture(scope='session')
def fpath_toy_hpo(fpath_test_data_dir: str) -> str:
    return os.path.join(fpath_test_data_dir, 'hp.toy.json')


@pytest.fixture(scope='session')
def toy_hpo(fpath_toy_hpo: str) -> hpotk.MinimalOntology:
    return hpotk.load_minimal_ontology(fpath_toy_hpo)


@pytest.fixture(scope='session')
def hpo(fpath_test_data_dir: str) -> hpotk.MinimalOntology:
    fpath_hpo = os.path.join(fpath_test_data_dir, 'hp.v2024-04-26.json.gz')
    return hpotk.load_minimal_ontology(fpath_hpo)


@pytest.fixture(scope='session')
def validation_runner(hpo: hpotk.MinimalOntology) -> hpotk.validate.ValidationRunner:
    validators = (
        hpotk.validate.ObsoleteTermIdsValidator(hpo),
        hpotk.validate.AnnotationPropagationValidator(hpo),
        hpotk.validate.PhenotypicAbnormalityValidator(hpo)
    )
    return hpotk.validate.ValidationRunner(validators)


def make_region(contig: str, start: int, end: int) -> GenomicRegion:
    a_contig = GRCh38.contig_by_name(contig)
    assert a_contig is not None
    return GenomicRegion(a_contig, start, end, Strand.POSITIVE)


@pytest.fixture(scope='session')
def protein_test_service() -> ProteinTestMetadataService:
    return ProteinTestMetadataService.create()


@pytest.fixture(scope='session')
def fpath_suox_cohort(
        fpath_test_data_dir: str,
) -> str:
    return os.path.join(fpath_test_data_dir, 'SUOX.json')


@pytest.fixture(scope='session')
def suox_cohort(
        fpath_suox_cohort: str,
) -> Cohort:
    with open(fpath_suox_cohort) as fh:
        return json.load(fh, cls=GpseaJSONDecoder)


@pytest.fixture(scope='session')
def suox_mane_tx_id() -> str:
    return 'NM_001032386.2'


@pytest.fixture(scope='session')
def suox_gt_predicate(
    suox_mane_tx_id: str,
) -> GenotypePolyPredicate:
    # To bin the patients to a group with >1 MISSENSE variant or 0 MISSENSE variants.
    
    return boolean_predicate(
        variant_predicate=VariantPredicates.variant_effect(
            effect=VariantEffect.MISSENSE_VARIANT,
            tx_id=suox_mane_tx_id
        )
    )


@pytest.fixture(scope='session')
def suox_pheno_predicates(
        hpo: hpotk.MinimalOntology,
) -> typing.Sequence[PhenotypePolyPredicate[hpotk.TermId]]:
    """
    Get predicates for test for presence of 5 HPO terms:
    Seizure, Ectopia lentis, Sulfocysteinuria, Neurodevelopmental delay, and Hypertonia.

    Note, these are just a *SUBSET* of all phenotypes that can be tested for in the *SUOX* cohort.
    """
    return (
        HpoPredicate(
            hpo=hpo,
            query=hpotk.TermId.from_curie('HP:0001250'),  # Seizure
        ),
        HpoPredicate(
            hpo=hpo,
            query=hpotk.TermId.from_curie('HP:0001083'),  # Ectopia lentis
        ),
        HpoPredicate(
            hpo=hpo,
            query=hpotk.TermId.from_curie('HP:0032350'),  # Sulfocysteinuria
        ),
        HpoPredicate(
            hpo=hpo,
            query=hpotk.TermId.from_curie('HP:0012758'),  # Neurodevelopmental delay
        ),
        HpoPredicate(
            hpo=hpo,
            query=hpotk.TermId.from_curie('HP:0001276'),  # Hypertonia
        ),
    )


@pytest.fixture
def hpo_result(
    suox_pheno_predicates: typing.Sequence[PhenotypePolyPredicate[hpotk.TermId]],
    suox_gt_predicate: GenotypePolyPredicate,
) -> HpoTermAnalysisResult:
    return HpoTermAnalysisResult(
        pheno_predicates=suox_pheno_predicates,
        n_usable=(40, 20, 30, 10, 100),
        all_counts=tuple(
            make_count_df(count, suox_gt_predicate, ph_pred)
            for count, ph_pred in zip(
                [
                    (10, 5, 15, 10),
                    (5, 2, 8, 5),
                    (10, 5, 5, 10),
                    (2, 3, 2, 3),
                    (10, 25, 35, 30),
                ],
                suox_pheno_predicates,
            )
        ),
        pvals=(0.01, 0.04, 0.3, 0.2, 0.7),
        corrected_pvals=(0.04, 0.1, 0.7, 0.5, 1.0),
        gt_predicate=suox_gt_predicate,
        mtc_filter_name="MTC filter name",
        mtc_filter_results=(
            PhenotypeMtcResult.ok(),
            PhenotypeMtcResult.ok(),
            PhenotypeMtcResult.ok(),
            PhenotypeMtcResult.ok(),
            PhenotypeMtcResult.ok(),
        ),
        mtc_name="MTC name",
    )


def make_count_df(
    counts: typing.Tuple[int, int, int, int],
    gt_predicate: GenotypePolyPredicate,
    ph_predicate: PhenotypePolyPredicate[hpotk.TermId],
) -> pd.DataFrame:
    rows = tuple(ph_predicate.get_categories())
    cols = tuple(gt_predicate.get_categories())
    data = np.array(counts).reshape((len(rows), len(cols)))
    return pd.DataFrame(
        data=data,
        index=pd.Index(rows),
        columns=pd.Index(cols),
    )


@pytest.fixture(scope='session')
def fpath_suox_tx_coordinates(fpath_test_data_dir: str) -> str:
    suox_mane_tx_id = 'NM_001032386.2'
    return os.path.join(fpath_test_data_dir, f'SUOX-{suox_mane_tx_id}.json')


@pytest.fixture(scope='session')
def suox_mane_tx_coordinates(
        fpath_suox_tx_coordinates: str,
) -> TranscriptCoordinates:
    with open(fpath_suox_tx_coordinates) as fh:
        return json.load(fh, cls=GpseaJSONDecoder)


@pytest.fixture(scope='session')
def fpath_suox_protein_metadata(fpath_test_data_dir: str) -> str:
    suox_mane_tx_protein_id = 'NP_001027558.1'
    return os.path.join(fpath_test_data_dir, f'SUOX-{suox_mane_tx_protein_id}.json')


@pytest.fixture(scope='session')
def suox_protein_metadata(
        fpath_suox_protein_metadata: str,
) -> ProteinMetadata:
    with open(fpath_suox_protein_metadata) as fh:
        return json.load(fh, cls=GpseaJSONDecoder)


@pytest.fixture(scope='session')
def toy_cohort(
        test_phenotypes: typing.Mapping[str, Phenotype],
        test_diseases: typing.Mapping[str, Disease],
) -> Cohort:
    prot_id = 'NP_037407.4'

    dup = Variant(
        variant_info=VariantInfo(
            variant_coordinates=VariantCoordinates(make_region("16", 89279849, 89279850), ref='G', alt='GC', change_length=1)
        ),
        tx_annotations=[
            TranscriptAnnotation(
                'ANKRD11', 'NM_013275.6', 'NM_013275.6:c.6691dup', False, [VariantEffect.FRAMESHIFT_VARIANT], [9],
                prot_id, "NP_001243112.1:p.Ala2231GlyfsTer29", Region(2230, 2231))
        ],
        genotypes=Genotypes.from_mapping({SampleLabels('HetSingleVar'): Genotype.HETEROZYGOUS}),
    )
    indel = Variant(
        variant_info=VariantInfo(
            variant_coordinates=VariantCoordinates(make_region("16", 89284600, 89284602), ref='GG', alt='A', change_length=-1),
        ),
        tx_annotations=[
            TranscriptAnnotation(
                'ANKRD11', 'NM_013275.6', 'NM_013275.6:c.1940_1941delinsT', False, [VariantEffect.FRAMESHIFT_VARIANT], [9], 
                prot_id, "NP_001243111.1:p.Ser647LeufsTer6", Region(646, 647))
        ],
        genotypes=Genotypes.from_mapping({SampleLabels('HetDoubleVar1'): Genotype.HETEROZYGOUS}),
    )
    snv_stop_gain = Variant(
        variant_info=VariantInfo(
            variant_coordinates=VariantCoordinates(make_region("16", 89280751, 89280752), ref='G', alt='T', change_length=0)
        ),
        tx_annotations=[
            TranscriptAnnotation(
                'ANKRD11', 'NM_013275.6', 'NM_013275.6:c.5790C>A', False, [VariantEffect.STOP_GAINED], [9], 
                prot_id, "NP_037407.4:p.Tyr1930Ter", Region(1929, 1930))
        ],
        genotypes=Genotypes.from_mapping({SampleLabels('HetDoubleVar1'): Genotype.HETEROZYGOUS}),
    )
    snv_missense = Variant(
        variant_info=VariantInfo(
            variant_coordinates=VariantCoordinates(make_region("16", 89275127, 89275128), ref='G', alt='A', change_length=0)
        ),
        tx_annotations=[
            TranscriptAnnotation(
                'ANKRD11', 'NM_013275.6', 'NM_013275.6:c.7534C>T', False, [VariantEffect.MISSENSE_VARIANT], [10], 
                prot_id, "NP_037407.4:p.Arg2512Trp", Region(2511, 2512))
        ],
        genotypes=Genotypes.from_mapping({SampleLabels('HetDoubleVar2'): Genotype.HETEROZYGOUS}),
    )
    del_frameshift = Variant(
        variant_info=VariantInfo(
            variant_coordinates=VariantCoordinates(make_region("16", 89279707, 89279725), ref='AGTGTTCGGGGCGGGGCC', alt='A', change_length=-17)
        ),
        tx_annotations=[
            TranscriptAnnotation(
                'ANKRD11', 'NM_013275.6', 'NM_013275.6:c.6817_6833del', False, [VariantEffect.FRAMESHIFT_VARIANT], [9], 
                prot_id, "NP_037407.4:p.Gly2273CysfsTer17", Region(2272, 2278),
            )
        ],
        genotypes=Genotypes.from_mapping({SampleLabels('HetDoubleVar2'): Genotype.HETEROZYGOUS}),
    )
    del_small = Variant(
        variant_info=VariantInfo(
            variant_coordinates=VariantCoordinates(make_region("16", 89279457, 89279459), ref='TG', alt='T', change_length=-1)
        ),
        tx_annotations=[
            TranscriptAnnotation(
                'ANKRD11', 'NM_013275.6', 'NM_013275.6:c.7083del', False, [VariantEffect.FRAMESHIFT_VARIANT], [9],
                prot_id, "NP_037407.4:p.Thr2362ProfsTer39", Region(2360, 2362))
        ],
        genotypes=Genotypes.from_mapping({SampleLabels('HomoVar'): Genotype.HOMOZYGOUS_ALTERNATE}),
    )
    del_large = Variant(
        variant_info=VariantInfo(
            variant_coordinates=VariantCoordinates(make_region("16", 89_190_070, 89_439_815), ref='N', alt='<DEL>', change_length=-249_745)
        ),
        tx_annotations=[
            TranscriptAnnotation(
                'ANKRD11', 'NM_013275.6', None, False, 
                [VariantEffect.STOP_LOST, VariantEffect.FEATURE_TRUNCATION, VariantEffect.CODING_SEQUENCE_VARIANT, VariantEffect.FIVE_PRIME_UTR_VARIANT, VariantEffect.THREE_PRIME_UTR_VARIANT, VariantEffect.INTRON_VARIANT], 
                [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13],
                prot_id, None, None
            )
        ],
        genotypes=Genotypes.from_mapping({SampleLabels('LargeCNV'): Genotype.HETEROZYGOUS}))

    patients = (
        Patient.from_raw_parts(
            SampleLabels('HetSingleVar'),
            sex=Sex.UNKNOWN_SEX,
            phenotypes=(
                test_phenotypes['arachnodactyly_T'],
                test_phenotypes['spasticity_F'],
                test_phenotypes['focal_clonic_seizure_T']
            ),
            variants=(dup,),
            diseases=(test_diseases['KBG_T'],)
        ),
        Patient.from_raw_parts(
            SampleLabels('HetDoubleVar1'),
            sex=Sex.UNKNOWN_SEX,
            phenotypes=(
                test_phenotypes['arachnodactyly_T'], test_phenotypes['seizure_T'], test_phenotypes['spasticity_T'],
            ),
            variants=(indel, snv_stop_gain),
            diseases=(test_diseases['KBG_T'],)
        ),
        Patient.from_raw_parts(
            SampleLabels('HetDoubleVar2'),
            sex=Sex.UNKNOWN_SEX,        
            phenotypes=(
                test_phenotypes['arachnodactyly_F'], test_phenotypes['spasticity_T'], test_phenotypes['seizure_T'],
            ),
            variants=(snv_missense, del_frameshift),
            diseases=(test_diseases['KBG_T'],)
        ),
        Patient.from_raw_parts(
            SampleLabels('HomoVar'),
            sex=Sex.UNKNOWN_SEX,
            phenotypes=(
                test_phenotypes['arachnodactyly_T'], test_phenotypes['spasticity_T'], test_phenotypes['seizure_T'],
            ),
            variants=(del_small,),
            diseases=()
        ),
        Patient.from_raw_parts(
            SampleLabels('LargeCNV'),
            sex=Sex.UNKNOWN_SEX,
            phenotypes=(
                test_phenotypes['arachnodactyly_T'], test_phenotypes['spasticity_T'], test_phenotypes['seizure_F'],
            ),
            variants=(del_large,),
            diseases=()
        ),
    )

    return Cohort.from_patients(patients)


@pytest.fixture(scope='session')
def test_phenotypes() -> typing.Mapping[str, Phenotype]:
    return {
        'arachnodactyly_T': Phenotype(hpotk.TermId.from_curie('HP:0001166'), True),
        'seizure_T': Phenotype(hpotk.TermId.from_curie('HP:0001250'), True),
        'focal_clonic_seizure_T': Phenotype(
            hpotk.TermId.from_curie('HP:0002266'), True,
        ),
        'spasticity_T': Phenotype(hpotk.TermId.from_curie('HP:0001257'), True),

        'arachnodactyly_F': Phenotype(hpotk.TermId.from_curie('HP:0001166'), False),
        'seizure_F': Phenotype(hpotk.TermId.from_curie('HP:0001250'), False),
        'spasticity_F': Phenotype(hpotk.TermId.from_curie('HP:0001257'), False),
        'focal_clonic_seizure_F': Phenotype(
            hpotk.TermId.from_curie('HP:0002266'), False,
        ),
    }


@pytest.fixture(scope='session')
def test_diseases() -> typing.Mapping[str, Disease]:
    return {
        'KBG_T': Disease(hpotk.TermId.from_curie("OMIM:148050"), "KBG syndrome", True),
        'KBG_F': Disease(hpotk.TermId.from_curie("OMIM:148050"), "KBG syndrome", False),
    }


@pytest.fixture(scope='session')
def genome_build() -> GenomeBuild:
    return GRCh38
