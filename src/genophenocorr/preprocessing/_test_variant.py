import json
import os

import pytest

from pkg_resources import resource_filename

from google.protobuf.json_format import Parse

# pyright: reportGeneralTypeIssues=false
from phenopackets import Phenopacket, GenomicInterpretation

from genophenocorr.model import VariantCoordinates
from genophenocorr.model.genome import GRCh38, GenomicRegion, Strand

from ._phenopacket import PhenopacketVariantCoordinateFinder
from ._protein import ProteinAnnotationCache, ProtCachingFunctionalAnnotator
from ._uniprot import UniprotProteinMetadataService
from ._variant import verify_start_end_coordinates, VepFunctionalAnnotator, VariantAnnotationCache, VarCachingFunctionalAnnotator


@pytest.mark.parametrize('contig_name, start, end, ref, alt, chlen, expected',
                         (['1', 100, 101, 'G', 'C', 0, '1:101-101/C'],  # SNP

                          ['1', 100, 102, 'GG', 'G', -1, '1:101-102/G'],  # DEL
                          ['16', 89_284_128, 89_284_134, 'CTTTTT', 'C', 5, '16:89284129-89284134/C'],  # DEL
                          ['16', 89_284_087, 89_284_089, 'AC', 'A', -1, '16:89284088-89284089/A'],  # DEL

                          #['1', 100, 101, 'G', 'GC', 1, '1:102-101/C'],  # INS
                          #['16', 89_283_999, 89_284_000, 'A', 'AT', 1, '16:89284001-89284000/T'],  # custom INS

                          ['1', 100, 120, 'G', 'C', 0, '1:101-120/C'],  # MNV
                          ['16', 89_283_999, 89_284_002, 'AAT', 'AGCG', 1, '16:89284000-89284002/AGCG'],  # custom MNV

                          # symbolic DEL
                          ['1', 100, 120, 'N', '<DEL>', -20, '1:101-120/DEL'],
                          ['9', 133_359_999, 133_360_011, 'N', '<DEL>', -12, '9:133360000-133360011/DEL'],

                          # symbolic DUP
                          ['1', 100, 103, 'N', '<DUP>', 3, '1:101-103/DUP'],
                          ['9', 133_359_999, 133_360_002, 'N', '<DUP>', 3, '9:133360000-133360002/DUP'],

                          # symbolic INS
                          # We can describe an insertion of 20 bases in two ways, depending on if
                          # the variant coordinates are trimmed. In the untrimmed case,
                          # the INS replaces the REF allele, in the trimmed case it does not.
                          # Therefore, the change length of the untrimmed is 19 - one base (`N`)
                          # is gone and 20 bases are inserted. The change length of the trimmed is 20 - zero bases
                          # are gone and 20 bases are inserted.
                          # In practical terms, the `cdna_start` and `cdna_end` fields of the VEP response will change
                          # depending on the trimming status.

                          # TODO:
                          # ['1', 100, 101, 'N', '<INS>', 19, '1:102-101/INS'],
                          ['1', 101, 101, '', '<INS>', 20, '1:102-101/INS'],

                          # ['9', 133_359_999, 133_360_000, 'N', '<INS>', 29, '9:133360001-133360000/INS'])
                          ['9', 133_360_000, 133_360_000, '', '<INS>', 30, '9:133360001-133360000/INS'])
                         )
def test_verify_start_end_coordinates(contig_name, start, end, ref, alt, chlen, expected):
    contig = GRCh38.contig_by_name(contig_name)
    region = GenomicRegion(contig, start, end, Strand.POSITIVE)
    vc = VariantCoordinates(region, ref, alt, chlen)
    out = verify_start_end_coordinates(vc)
    assert out == expected


@pytest.fixture
def pp_vc_finder() -> PhenopacketVariantCoordinateFinder:
    return PhenopacketVariantCoordinateFinder(GRCh38)


@pytest.mark.parametrize("pp_path, expected",
                         [('test_data/deletion_test.json', '16_89284129_89284134_CTTTTT_C'),
                          ('test_data/insertion_test.json', '16_89280829_89280830_C_CA'),
                          ('test_data/missense_test.json', '16_89279135_89279135_G_C'),
                          ('test_data/duplication_test.json', '16_89279850_89279851_G_GC'),
                          ('test_data/delinsert_test.json', '16_89284601_89284602_GG_A'),
                          ('test_data/CVDup_test.json', '16_89284524_89373231_DUP'),
                          ('test_data/CVDel_test.json', '16_89217282_89506042_DEL')
                          ])
def test_find_coordinates(pp_path, expected, pp_vc_finder):
    fname = resource_filename(__name__, pp_path)
    gi = read_genomic_interpretation_json(fname)

    vc, gt = pp_vc_finder.find_coordinates(gi)

    assert expected == vc.variant_key


def read_genomic_interpretation_json(fpath: str) -> GenomicInterpretation:
    with open(fpath) as fh:
        return Parse(fh.read(), GenomicInterpretation())

@pytest.fixture
def variant_annotator(tmp_path):
    pm = UniprotProteinMetadataService()
    pac = ProteinAnnotationCache(tmp_path)
    pfa = ProtCachingFunctionalAnnotator(pac, pm)
    return VepFunctionalAnnotator(pfa)

@pytest.fixture
def caching_annotator(variant_annotator, tmp_path):
    vac = VariantAnnotationCache(tmp_path)
    return VarCachingFunctionalAnnotator(vac, variant_annotator)

def test_caching_full_circle(caching_annotator, pp_vc_finder, variant_annotator):
    fname = resource_filename(__name__, 'test_data/missense_test.json')
    gi = read_genomic_interpretation_json(fname)
    var_coords, gt = pp_vc_finder.find_coordinates(gi)
    var_anno_results = variant_annotator.annotate(var_coords)
    cache_anno_results = caching_annotator.annotate(var_coords)
    assert var_anno_results == cache_anno_results
    assert cache_anno_results == caching_annotator.annotate(var_coords)

@pytest.fixture
def oldfile_cache_annotator(variant_annotator):
    data = resource_filename(__name__, 'test_data/annotations')
    vac = VariantAnnotationCache(data)
    return VarCachingFunctionalAnnotator(vac, variant_annotator)


@pytest.mark.skip("Skipped until new Pickle files are generated")
def test_cache_from_older_file(oldfile_cache_annotator, pp_vc_finder, variant_annotator):
    fname = resource_filename(__name__, 'test_data/missense_test.json')
    gi = read_genomic_interpretation_json(fname)
    var_coords, gt = pp_vc_finder.find_coordinates(gi)
    var_anno_results = variant_annotator.annotate(var_coords)
    cached_file_results = oldfile_cache_annotator.annotate(var_coords)
    assert var_anno_results == cached_file_results


class TestVepFunctionalAnnotator:

    TEST_DATA_DIR = resource_filename(__name__, os.path.join('test_data', 'vep_response'))

    def test__process_item_missense(self, variant_annotator: VepFunctionalAnnotator):
        response = self._load_response_json('missense.json')
        first = response[0]

        annotations = list(filter(lambda i: i is not None, map(lambda item: variant_annotator._process_item(item), first['transcript_consequences'])))

        ann_by_tx = {ann.transcript_id: ann for ann in annotations}

        assert {'NM_013275.6', 'NM_001256183.2', 'NM_001256182.2'} == set(ann_by_tx.keys())

        preferred = ann_by_tx['NM_013275.6']
        assert preferred.transcript_id == 'NM_013275.6'
        assert preferred.is_preferred == True
        assert preferred.hgvsc_id == 'NM_013275.6:c.7407C>G'
        assert preferred.variant_effects == ('stop_gained',)

    def test__process_item_deletion(self, variant_annotator: VepFunctionalAnnotator):
        response = self._load_response_json('deletion.json')
        first = response[0]

        annotations = [variant_annotator._process_item(item) for item in first['transcript_consequences']]

        # TODO - finish
        for a in annotations:
            if a is not None:
                print(a)
        print('Done')

    def _load_response_json(self, test_name: str):
        response_fpath = os.path.join(self.TEST_DATA_DIR, test_name)
        with open(response_fpath) as fh:
            return json.load(fh)
