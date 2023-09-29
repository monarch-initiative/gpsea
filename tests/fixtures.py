import os

import hpotk
import pytest

from genophenocorr.model import *
from genophenocorr.model.genome import GRCh38, GenomicRegion, Region, Strand


def make_region(contig: str, start: int, end: int) -> GenomicRegion:
    return GenomicRegion(GRCh38.contig_by_name(contig), start, end, Strand.POSITIVE)


@pytest.fixture
def toy_cohort() -> Cohort:

    prot = ProteinMetadata(protein_id='NP_037407.4', label='Ankyrin repeat domain-containing protein 11', 
            protein_features=(ProteinFeature.create(feature_type=FeatureType.REPEAT, info=FeatureInfo('ANK 1', Region(start=167, end=196))),
            ProteinFeature.create(feature_type=FeatureType.REPEAT, info=FeatureInfo('ANK 2', Region(start=200, end=229))),
            ProteinFeature.create(feature_type=FeatureType.REPEAT, info=FeatureInfo('ANK 3', Region(start=233, end=262))),
            ProteinFeature.create(feature_type=FeatureType.REPEAT, info=FeatureInfo('ANK 4', Region(start=266, end=292))),
            ProteinFeature.create(feature_type=FeatureType.REGION, info=FeatureInfo('Disordered', Region(start=1, end=90))),
            ProteinFeature.create(feature_type=FeatureType.REGION, info=FeatureInfo('Disordered', Region(start=128, end=169))),
            ProteinFeature.create(feature_type=FeatureType.REGION, info=FeatureInfo('Disordered', Region(start=289, end=380))),
            ProteinFeature.create(feature_type=FeatureType.REGION, info=FeatureInfo('Disordered', Region(start=398, end=647))),
            ProteinFeature.create(feature_type=FeatureType.REGION, info=FeatureInfo('Disordered', Region(start=723, end=783))),
            ProteinFeature.create(feature_type=FeatureType.REGION, info=FeatureInfo('Disordered', Region(start=881, end=1043))),
            ProteinFeature.create(feature_type=FeatureType.REGION, info=FeatureInfo('Disordered', Region(start=1059, end=1393))),
            ProteinFeature.create(feature_type=FeatureType.REGION, info=FeatureInfo('Disordered', Region(start=1424, end=1710))),
            ProteinFeature.create(feature_type=FeatureType.REGION, info=FeatureInfo('Disordered', Region(start=1814, end=1836))),
            ProteinFeature.create(feature_type=FeatureType.REGION, info=FeatureInfo('Disordered', Region(start=1988, end=2019))),
            ProteinFeature.create(feature_type=FeatureType.REGION, info=FeatureInfo('Disordered', Region(start=2131, end=2406))),
            ProteinFeature.create(feature_type=FeatureType.REGION, info=FeatureInfo('Important for protein degradation', Region(start=2369, end=2663)))))

    HetSingleVar = [Variant('16_89279851_-/C', 'insertion', 
                VariantCoordinates(make_region("16", 89279849, 89279851),
                                   ref='G', alt='GC', change_length=1),
                [TranscriptAnnotation('ANKRD11', 'NM_013275.6', 'NM_013275.6:c.6691dup', False, ['frameshift_variant'],
                [9], [prot], 2231, 2231)])]
    HetDoubleVar1 = [Variant('16_89284601_GG/A', 'indel', 
                VariantCoordinates(make_region("16", 89284600, 89284602),
                                   ref='GG', alt='A', change_length=-1),
                [TranscriptAnnotation('ANKRD11', 'NM_013275.6', 'NM_013275.6:c.1940_1941delinsT', False, ['frameshift_variant'],
                [9], [prot], 647, 647)]),
                Variant('16_89280752_G/T', 'SNV', 
                VariantCoordinates(make_region("16", 89280751, 89280752),
                                   ref='G', alt='T', change_length=0),
                [TranscriptAnnotation('ANKRD11', 'NM_013275.6', 'NM_013275.6:c.5790C>A', False, ['stop_gained'],
                [9], [prot], 1930, 1930)])]
    HetDoubleVar2 = [Variant('16_89275128_G/A', 'SNV', 
                VariantCoordinates(make_region("16", 89275127, 89275128),
                                   ref='G', alt='A', change_length=0),
                [TranscriptAnnotation('ANKRD11', 'NM_013275.6', 'NM_013275.6:c.7534C>T', False, ['missense_variant'],
                [10], [prot], 2512, 2512)]),
                Variant('16_89279708_AGTGTTCGGGGCGGGGCC/A', 'indel', 
                VariantCoordinates(make_region("16", 89279707, 89279725),
                                   ref='AGTGTTCGGGGCGGGGCC', alt='A', change_length=-17),
                [TranscriptAnnotation('ANKRD11', 'NM_013275.6', 'NM_013275.6:c.6817_6833del', False, ['frameshift_variant'],
                [9], [prot], 2273, 2278)])]
    HomoVar = [Variant('16_89279458_TG/T', 'indel', 
                VariantCoordinates(make_region("16", 89279457, 89279459),
                                   ref='TG', alt='T', change_length=-1),
                [TranscriptAnnotation('ANKRD11', 'NM_013275.6', 'NM_013275.6:c.7083del', False, ['frameshift_variant'],
                [9], [prot], 2361, 2362)])]
    LargeCNV = [Variant('16_89190071_deletion', 'deletion', 
                VariantCoordinates(make_region("16", 89190070, 89439815),
                                   ref='N', alt='<DEL>', change_length=4),
                [TranscriptAnnotation('ANKRD11', 'NM_013275.6', None, False, ['stop_lost', 'feature_truncation', 'coding_sequence_variant', '5_prime_UTR_variant', '3_prime_UTR_variant', 'intron_variant'],
                [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13], [prot], None, None)])]

    phenos = get_test_phenotypes()

    patients = (
        Patient('HetSingleVar',
                phenotypes=(phenos['arachnodactyly_T'], phenos['spasticity_F'], phenos['focal_clonic_seizure_T']),
                variants=HetSingleVar,
                proteins=[prot]
                ),
        Patient('HetDoubleVar1',
                phenotypes=(phenos['arachnodactyly_T'], phenos['seizure_T'], phenos['spasticity_T']),
                variants=HetDoubleVar1,
                proteins=[prot]
                ),
        Patient('HetDoubleVar2',
                phenotypes=(phenos['arachnodactyly_F'], phenos['spasticity_T'], phenos['seizure_T']),
                variants=HetDoubleVar2,
                proteins=[prot]
                ),
        Patient('HomoVar',
                phenotypes=(phenos['arachnodactyly_T'], phenos['spasticity_T'], phenos['seizure_T']),
                variants=HomoVar,
                proteins=[prot]
                ),
        Patient('LargeCNV',
                phenotypes=(phenos['arachnodactyly_T'], phenos['spasticity_T'], phenos['seizure_F']),
                variants=LargeCNV,
                proteins=[prot]
                ),
    )

    return Cohort.from_patients(patients)


def get_test_phenotypes():
    phenotypes = {}

    phenotypes['arachnodactyly_T'] = Phenotype(hpotk.TermId.from_curie('HP:0001166'), "Arachnodactyly", True)
    phenotypes['seizure_T'] = Phenotype(hpotk.TermId.from_curie('HP:0001250'), "Seizure", True)
    phenotypes['focal_clonic_seizure_T'] = Phenotype(hpotk.TermId.from_curie('HP:0002266'), "Focal clonic seizure", True)
    phenotypes['spasticity_T'] = Phenotype(hpotk.TermId.from_curie('HP:0001257'), "Spasticity", True)
    phenotypes['arachnodactyly_F'] = Phenotype(hpotk.TermId.from_curie('HP:0001166'), "Arachnodactyly", False)
    phenotypes['seizure_F'] = Phenotype(hpotk.TermId.from_curie('HP:0001250'), "Seizure", False)
    phenotypes['spasticity_F'] = Phenotype(hpotk.TermId.from_curie('HP:0001257'), "Spasticity", False)
    phenotypes['focal_clonic_seizure_F'] = Phenotype(hpotk.TermId.from_curie('HP:0002266'), "Focal clonic seizure", False)

    return phenotypes


@pytest.fixture
def toy_hpo() -> hpotk.Ontology:
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'testingDefaults', 'hp.toy.json')
    return hpotk.ontology.load.obographs.load_ontology(path)
