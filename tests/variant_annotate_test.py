import pytest
import json
from genophenocorr import variant
from phenopackets import GenomicInterpretation
from google.protobuf.json_format import Parse


@pytest.fixture
def VariantCoordFinder():
    return variant.PhenopacketVariantCoordinateFinder()

@pytest.fixture
def VariantAnnotatorFixure():
    return variant.VepFunctionalAnnotator()

@pytest.fixture
def CachingAnnotator():
    return variant.CachingFunctionalAnnotator('testSamples', variant.VariantAnnotationCache(), variant.VepFunctionalAnnotator())
    

@pytest.mark.parametrize('patient, start, hgvsc, effects',
    [('testSamples/deletion_test.json', 89284129, 'NM_013275.6:c.2408_2412del', None),
    ('testSamples/insertion_test.json', 89280829, 'NM_013275.6:c.5712_5713insT', None),
    ('testSamples/missense_test.json', 89279135, 'NM_013275.6:c.7407C>G', None),
    ('testSamples/duplication_test.json', 89279850, 'NM_013275.6:c.6691dup', None),
    ('testSamples/delinsert_test.json', 89284601, 'NM_013275.6:c.1940_1941delinsT', None),
    ('testSamples/CVDup_test.json', 89284523, None, [
                    "coding_sequence_variant",
                    "5_prime_UTR_variant",
                    "intron_variant",
                    "feature_elongation"
                ]),
    ('testSamples/CVDel_test.json', 89217281, None, ['transcript_ablation'])
    ])

def test_non_CNV(patient, start, hgvsc, effects, VariantCoordFinder, VariantAnnotatorFixure):
    with open(patient) as f:
        data = f.read()
        data = json.loads(data)
    phenopack = Parse(json.dumps(data), GenomicInterpretation())
    var_coords = VariantCoordFinder.find_coordinates(phenopack)
    assert var_coords.start == start
    var_anno = VariantAnnotatorFixure.annotate(var_coords)
    if hgvsc is not None:
        trans = [tx_anno.hgvsc_id for tx_anno in var_anno.tx_annotations if tx_anno.transcript_id == 'NM_013275.6'][0]
        assert trans == hgvsc
    elif effects is not None:
        trans = [tx_anno.variant_effects for tx_anno in var_anno.tx_annotations if tx_anno.transcript_id == 'NM_013275.6'][0]
        assert trans == effects

def test_caching(VariantCoordFinder, VariantAnnotatorFixure, CachingAnnotator):
    test_patient = 'testSamples/deletion_test.json'
    with open(test_patient) as f:
        data = f.read()
    data = json.loads(data)
    phenopack = Parse(json.dumps(data), GenomicInterpretation())
    var_coords = VariantCoordFinder.find_coordinates(phenopack)
    var_anno_uncached = VariantAnnotatorFixure.annotate(var_coords)
    var_anno_cached = CachingAnnotator.annotate(var_coords)
    assert var_anno_cached.__eq__(var_anno_uncached)
    var_anno_recached = CachingAnnotator.annotate(var_coords)
    assert var_anno_uncached.__eq__(var_anno_recached)
