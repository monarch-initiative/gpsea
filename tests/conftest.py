import json
import os
import typing

import hpotk
import pytest

from ._protein_test_service import ProteinTestMetadataService
from genophenocorr.model import *
from genophenocorr.model.genome import GRCh38, GenomicRegion, Region, Strand, GenomeBuild
from genophenocorr.io import GenophenocorrJSONEncoder, GenophenocorrJSONDecoder



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

def fpath_test_data() -> str:
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), 'test_data')


@pytest.fixture(scope='session')
def fpath_toy_hpo(fpath_test_data: str) -> str:
    return os.path.join(fpath_test_data, 'hp.toy.json')


@pytest.fixture(scope='session')
def fpath_test_zn462_human_uniprot(fpath_test_data: str) -> str:
    return os.path.join(fpath_test_data, "uniprot", "ZN462_HUMAN.json")


@pytest.fixture(scope='session')
def toy_hpo(fpath_toy_hpo: str) -> hpotk.MinimalOntology:
    return hpotk.load_minimal_ontology(fpath_toy_hpo)


@pytest.fixture(scope='session')
def genome_build_hg38() -> GenomeBuild:
    return GRCh38



@pytest.fixture(scope='session')
def hpo(fpath_test_data_dir: str) -> hpotk.MinimalOntology:
    fpath_hpo = os.path.join(fpath_test_data_dir, 'hp.v2024-04-26.json.gz')
    return hpotk.load_minimal_ontology(fpath_hpo)


@pytest.fixture(scope='session')
def toy_validation_runner(hpo: hpotk.MinimalOntology) -> hpotk.validate.ValidationRunner:
    validators = (
        hpotk.validate.ObsoleteTermIdsValidator(hpo),
        hpotk.validate.AnnotationPropagationValidator(hpo),
        hpotk.validate.PhenotypicAbnormalityValidator(hpo)
    )
    return hpotk.validate.ValidationRunner(validators)

def make_region(contig: str, start: int, end: int) -> GenomicRegion:
    return GenomicRegion(GRCh38.contig_by_name(contig), start, end, Strand.POSITIVE)


@pytest.fixture(scope='session')
def protein_test_service() -> ProteinTestMetadataService:
    return ProteinTestMetadataService()


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
        return json.load(fh, cls=GenophenocorrJSONDecoder)


@pytest.mark.skip('Run manually to regenerate `suox_cohort`')
def test_regenerate_cohort(
        fpath_suox_cohort: str,
        hpo: hpotk.MinimalOntology,
):
    """
    The test for regenerating the `SUOX.json` file based on a cohort of phenopackets.
    The test needs path to a folder with phenopacket JSON files (empty `str` below).
    """
    fpath_suox_pp_dir = '/path/to/SUOX/phenopackets'

    from genophenocorr.preprocessing import configure_caching_cohort_creator, load_phenopacket_folder

    cohort_creator = configure_caching_cohort_creator(hpo)
    cohort = load_phenopacket_folder(fpath_suox_pp_dir, cohort_creator, validation_policy='strict')
    with open(fpath_suox_cohort, 'w') as fh:
        json.dump(cohort, fh, cls=GenophenocorrJSONEncoder, indent=2)


@pytest.fixture(scope='session')
def toy_cohort(
        test_phenotypes: typing.Mapping[str, Phenotype],
        test_diseases: typing.Mapping[str, Disease],
) -> Cohort:
    prot_id = 'NP_037407.4'

    dup = Variant(VariantCoordinates(make_region("16", 89279849, 89279850), ref='G', alt='GC', change_length=1),
                  [
                      TranscriptAnnotation('ANKRD11', 'NM_013275.6', 'NM_013275.6:c.6691dup', False, [VariantEffect.FRAMESHIFT_VARIANT], [9],
                                           prot_id, Region(2230, 2231))
                  ],
                  Genotypes.from_mapping({SampleLabels('HetSingleVar'): Genotype.HETEROZYGOUS}))
    indel = Variant(VariantCoordinates(make_region("16", 89284600, 89284602), ref='GG', alt='A', change_length=-1),
                    [
                        TranscriptAnnotation('ANKRD11', 'NM_013275.6', 'NM_013275.6:c.1940_1941delinsT', False, [VariantEffect.FRAMESHIFT_VARIANT],
                                             [9], prot_id, Region(646, 647))
                    ],
                    Genotypes.from_mapping({SampleLabels('HetDoubleVar1'): Genotype.HETEROZYGOUS}))
    snv_stop_gain = Variant(VariantCoordinates(make_region("16", 89280751, 89280752), ref='G', alt='T', change_length=0),
                            [
                                TranscriptAnnotation('ANKRD11', 'NM_013275.6', 'NM_013275.6:c.5790C>A', False, [VariantEffect.STOP_GAINED], [9], prot_id,
                             Region(1929, 1930))],
                            Genotypes.from_mapping({SampleLabels('HetDoubleVar1'): Genotype.HETEROZYGOUS}))
    snv_missense = Variant(VariantCoordinates(make_region("16", 89275127, 89275128), ref='G', alt='A', change_length=0),
                           [
                               TranscriptAnnotation('ANKRD11', 'NM_013275.6', 'NM_013275.6:c.7534C>T', False, [VariantEffect.MISSENSE_VARIANT], [10],
                             prot_id, Region(2511, 2512))
                           ],
                           Genotypes.from_mapping({SampleLabels('HetDoubleVar2'): Genotype.HETEROZYGOUS}))
    del_frameshift = Variant(VariantCoordinates(make_region("16", 89279707, 89279725), ref='AGTGTTCGGGGCGGGGCC', alt='A', change_length=-17),
                             [
                                 TranscriptAnnotation('ANKRD11', 'NM_013275.6', 'NM_013275.6:c.6817_6833del', False, [VariantEffect.FRAMESHIFT_VARIANT],
                              [9], prot_id, Region(2272, 2278))
                             ],
                             Genotypes.from_mapping({SampleLabels('HetDoubleVar2'): Genotype.HETEROZYGOUS}))
    del_small = Variant(VariantCoordinates(make_region("16", 89279457, 89279459), ref='TG', alt='T', change_length=-1),
                        [
                            TranscriptAnnotation('ANKRD11', 'NM_013275.6', 'NM_013275.6:c.7083del', False, [VariantEffect.FRAMESHIFT_VARIANT], [9],
                             prot_id, Region(2360, 2362))
                        ],
                        Genotypes.from_mapping({SampleLabels('HomoVar'): Genotype.HOMOZYGOUS_ALTERNATE}))
    del_large = Variant(VariantCoordinates(make_region("16", 89_190_070, 89_439_815), ref='N', alt='<DEL>', change_length=-249_745),
                        [
                            TranscriptAnnotation('ANKRD11', 'NM_013275.6', None, False,
                                 [VariantEffect.STOP_LOST, VariantEffect.FEATURE_TRUNCATION, VariantEffect.CODING_SEQUENCE_VARIANT, VariantEffect.FIVE_PRIME_UTR_VARIANT,
                                  VariantEffect.THREE_PRIME_UTR_VARIANT, VariantEffect.INTRON_VARIANT], [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13],
                                 prot_id, None)
                        ],
                        Genotypes.from_mapping({SampleLabels('LargeCNV'): Genotype.HETEROZYGOUS}))

    patients = (
        Patient(SampleLabels('HetSingleVar'),
                phenotypes=(
                    test_phenotypes['arachnodactyly_T'],
                    test_phenotypes['spasticity_F'],
                    test_phenotypes['focal_clonic_seizure_T']),
                variants=(dup,),
                diseases=(test_diseases['KBG_T'],)
                ),
        Patient(SampleLabels('HetDoubleVar1'),
                phenotypes=(
                    test_phenotypes['arachnodactyly_T'], test_phenotypes['seizure_T'], test_phenotypes['spasticity_T'],
                ),
                variants=(indel, snv_stop_gain),
                diseases=(test_diseases['KBG_T'],)
                ),
        Patient(SampleLabels('HetDoubleVar2'),
                phenotypes=(
                    test_phenotypes['arachnodactyly_F'], test_phenotypes['spasticity_T'], test_phenotypes['seizure_T'],
                ),
                variants=(snv_missense, del_frameshift),
                diseases=(test_diseases['KBG_T'],)
                ),
        Patient(SampleLabels('HomoVar'),
                phenotypes=(
                    test_phenotypes['arachnodactyly_T'], test_phenotypes['spasticity_T'], test_phenotypes['seizure_T'],
                ),
                variants=(del_small,),
                diseases=()
                ),
        Patient(SampleLabels('LargeCNV'),
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
        'arachnodactyly_T': Phenotype(hpotk.TermId.from_curie('HP:0001166'), "Arachnodactyly", True),
        'seizure_T': Phenotype(hpotk.TermId.from_curie('HP:0001250'), "Seizure", True),
        'focal_clonic_seizure_T': Phenotype(
            hpotk.TermId.from_curie('HP:0002266'), "Focal clonic seizure", True,
        ),
        'spasticity_T': Phenotype(hpotk.TermId.from_curie('HP:0001257'), "Spasticity", True),

        'arachnodactyly_F': Phenotype(hpotk.TermId.from_curie('HP:0001166'), "Arachnodactyly", False),
        'seizure_F': Phenotype(hpotk.TermId.from_curie('HP:0001250'), "Seizure", False),
        'spasticity_F': Phenotype(hpotk.TermId.from_curie('HP:0001257'), "Spasticity", False),
        'focal_clonic_seizure_F': Phenotype(
            hpotk.TermId.from_curie('HP:0002266'), "Focal clonic seizure", False,
        ),
    }


@pytest.fixture(scope='session')
def test_diseases() -> typing.Mapping[str, Disease]:
    return {
        'KBG_T': Disease(hpotk.TermId.from_curie("OMIM:148050"), "KBG syndrome", True),
        'KBG_F': Disease(hpotk.TermId.from_curie("OMIM:148050"), "KBG syndrome", False),
    }
