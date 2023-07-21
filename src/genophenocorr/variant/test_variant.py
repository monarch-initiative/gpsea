
import pytest
from google.protobuf.json_format import Parse
# pyright: reportGeneralTypeIssues=false
from phenopackets import Phenopacket, GenomicInterpretation

from ._annotators import VariantCoordinates, verify_start_end_coordinates, PhenopacketVariantCoordinateFinder, VepFunctionalAnnotator, VariantAnnotationCache, VarCachingFunctionalAnnotator
from genophenocorr.protein import UniprotProteinMetadataService, ProteinAnnotationCache, ProtCachingFunctionalAnnotator
from pkg_resources import resource_filename


@pytest.mark.parametrize('contig, start, end, ref, alt, chlen, expected',
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
def test_verify_start_end_coordinates(contig, start, end, ref, alt, chlen, expected):
    vc = VariantCoordinates(contig, start, end, ref, alt, chlen, 'Whatever')
    out = verify_start_end_coordinates(vc)
    assert out == expected


@pytest.fixture
def pp_vc_finder() -> PhenopacketVariantCoordinateFinder:
    return PhenopacketVariantCoordinateFinder()


@pytest.mark.parametrize("pp_path, expected",
                         [('test_data/deletion_test.json', '16_89284128_89284134_CTTTTT_C_heterozygous'),
                          ('test_data/insertion_test.json', '16_89280828_89280830_C_CA_heterozygous'),
                          ('test_data/missense_test.json', '16_89279134_89279135_G_C_heterozygous'),
                          ('test_data/duplication_test.json', '16_89279849_89279851_G_GC_heterozygous'),
                          ('test_data/delinsert_test.json', '16_89284600_89284602_GG_A_heterozygous'),
                          ('test_data/CVDup_test.json', '16_89284523_89373231_N_<DUP>_heterozygous'),
                          ('test_data/CVDel_test.json', '16_89217281_89506042_N_<DEL>_heterozygous')
                          ])
def test_find_coordinates(pp_path, expected, pp_vc_finder):
    fname = resource_filename(__name__, pp_path)
    gi = read_genomic_interpretation_json(fname)

    assert pp_vc_finder.find_coordinates(gi).as_string() == expected


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
    var_coords = pp_vc_finder.find_coordinates(gi)
    var_anno_results = variant_annotator.annotate(var_coords)
    cache_anno_results = caching_annotator.annotate(var_coords)
    assert var_anno_results == cache_anno_results
    assert cache_anno_results == caching_annotator.annotate(var_coords)

@pytest.fixture
def oldfile_cache_annotator(variant_annotator):
    data = resource_filename(__name__, 'test_data/annotations')
    vac = VariantAnnotationCache(data)
    return VarCachingFunctionalAnnotator(vac, variant_annotator)

def test_cache_from_older_file(oldfile_cache_annotator, pp_vc_finder, variant_annotator):
    fname = resource_filename(__name__, 'test_data/missense_test.json')
    gi = read_genomic_interpretation_json(fname)
    var_coords = pp_vc_finder.find_coordinates(gi)
    var_anno_results = variant_annotator.annotate(var_coords)
    cached_file_results = oldfile_cache_annotator.annotate(var_coords)
    assert var_anno_results == cached_file_results

