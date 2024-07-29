import json
import os
import typing

import pytest

from genophenocorr.model import VariantCoordinates, VariantEffect
from genophenocorr.model.genome import GenomicRegion, Strand, GenomeBuild

from genophenocorr.preprocessing import VepFunctionalAnnotator


LMNA_MANE_TX_ID = 'NM_170707.4'
ANKRD11_MANE_TX_ID = 'NM_013275.6'


@pytest.mark.parametrize(
    'contig_name, start, end, ref, alt, chlen, expected',
    (
            ['1', 100, 101, 'G', 'C', 0, '1:101-101/C'],  # SNP

            ['1', 100, 102, 'GG', 'G', -1, '1:101-102/G'],  # DEL
            ['16', 89_284_128, 89_284_134, 'CTTTTT', 'C', 5, '16:89284129-89284134/C'],  # DEL
            ['16', 89_284_087, 89_284_089, 'AC', 'A', -1, '16:89284088-89284089/A'],  # DEL

            # ['1', 100, 101, 'G', 'GC', 1, '1:102-101/C'],  # INS
            # ['16', 89_283_999, 89_284_000, 'A', 'AT', 1, '16:89284001-89284000/T'],  # custom INS

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
            ['9', 133_360_000, 133_360_000, '', '<INS>', 30, '9:133360001-133360000/INS']),
)
def test_verify_start_end_coordinates(
        contig_name, start, end, ref, alt, chlen, expected,
        genome_build: GenomeBuild,
):
    contig = genome_build.contig_by_name(contig_name)
    region = GenomicRegion(contig, start, end, Strand.POSITIVE)
    vc = VariantCoordinates(region, ref, alt, chlen)
    out = VepFunctionalAnnotator.format_coordinates_for_vep_query(vc)
    assert out == expected


class TestVepFunctionalAnnotator:
    """
    The test suite mainly works with VEP responses stored in a test folder.
    The test that regenerates the files is skipped and designed to be run manually.
    """

    @pytest.fixture(scope='class')
    def fpath_vep_response_dir(
        self,
        fpath_preprocessing_data_dir: str,
    ) -> str:
        return os.path.join(fpath_preprocessing_data_dir, 'vep_response')

    @pytest.fixture
    def variant_annotator(
            self,
    ) -> VepFunctionalAnnotator:
        return VepFunctionalAnnotator()

    def test_process_response_missense(
        self,
        variant_annotator: VepFunctionalAnnotator,
        fpath_vep_response_dir: str,
    ):
        # ('16', 89_279_134, 89_279_135, 'G', 'C', 0, ),  # `16_89279135_89279135_G_C`
        variant_key = '16_89279135_89279135_G_C'
        response_fpath = os.path.join(fpath_vep_response_dir, f'{variant_key}.json')
        response = load_json_response(response_fpath)

        annotations = variant_annotator.process_response(variant_key, response)

        # we should have 3 non-computational transcripts
        anns = [a for a in annotations if a is not None]
        assert len(anns) == 3

        ann_by_tx = {ann.transcript_id: ann for ann in anns}

        assert {'NM_013275.6', 'NM_001256183.2', 'NM_001256182.2'} == set(ann_by_tx.keys())

        preferred = ann_by_tx[ANKRD11_MANE_TX_ID]
        assert preferred.transcript_id == 'NM_013275.6'
        assert preferred.is_preferred is True
        assert preferred.hgvs_cdna == 'NM_013275.6:c.7407C>G'
        assert preferred.variant_effects == (VariantEffect.STOP_GAINED,)
        assert preferred.overlapping_exons == (9,)

    def test_process_response_deletion(
            self,
            variant_annotator: VepFunctionalAnnotator,
            fpath_vep_response_dir: str,
    ):
        # ('16', 89_284_128, 89_284_134, 'CTTTTT', 'C', -5, ),  # `16_89284129_89284134_CTTTTT_C`
        variant_key = '16_89284129_89284134_CTTTTT_C'
        response_fpath = os.path.join(fpath_vep_response_dir, f'{variant_key}.json')
        response = load_json_response(response_fpath)

        annotations = variant_annotator.process_response(variant_key, response)

        # we should have 3 non-computational transcripts
        anns = [a for a in annotations if a is not None]
        assert len(anns) == 3

        ann_by_tx = {ann.transcript_id: ann for ann in anns}

        assert {'NM_013275.6', 'NM_001256183.2', 'NM_001256182.2'} == set(ann_by_tx.keys())

        preferred = ann_by_tx[ANKRD11_MANE_TX_ID]
        assert preferred.transcript_id == 'NM_013275.6'
        assert preferred.is_preferred is True
        assert preferred.hgvs_cdna == 'NM_013275.6:c.2408_2412del'
        assert preferred.variant_effects == (VariantEffect.FRAMESHIFT_VARIANT,)
        assert preferred.overlapping_exons == (9,)

    @pytest.mark.parametrize(
        'contig_name, start, end, ref, alt, chlen, mane_tx_id,'
        'exp_protein_id, exp_protein_start, exp_protein_end',
        [
            (  # `1_156114921_156114920_G_A`
                    '1', 156_114_920, 156_114_920, 'G', 'A', 0, LMNA_MANE_TX_ID,
                    'NP_733821.1', 0, 1,
            ),
            (  # `1_156114215_156114218_CGCC_C`
                    '1', 156_115_214, 156_115_218, 'CGCC', 'C', -3, LMNA_MANE_TX_ID,
                    'NP_733821.1', 99, 100,
            ),
            (  # `16_89284129_89284130_CT_C`
                    '16', 89_284_128, 89_284_130, 'CT', 'C', -1, ANKRD11_MANE_TX_ID,
                    'NP_037407.4', 803, 804,
            ),
        ]
    )
    def test_parse_response(
            self,
            fpath_vep_response_dir: str,
            variant_annotator: VepFunctionalAnnotator,
            genome_build: GenomeBuild,
            contig_name: str, start: int, end: int, ref: str, alt: str, chlen: int, mane_tx_id: str,
            exp_protein_id: str, exp_protein_start: int, exp_protein_end: int,
    ):
        contig = genome_build.contig_by_name(contig_name)
        region = GenomicRegion(contig, start, end, Strand.POSITIVE)
        vc = VariantCoordinates(region, ref, alt, chlen)
        fpath_response = os.path.join(fpath_vep_response_dir, f'{vc.variant_key}.json')
        response = load_json_response(fpath_response)

        annotations = variant_annotator.process_response(vc.variant_key, response)
        ann_by_tx = {a.transcript_id: a for a in annotations}
        ann = ann_by_tx[mane_tx_id]

        assert ann.protein_id == exp_protein_id
        assert ann.protein_effect_location.start == exp_protein_start
        assert ann.protein_effect_location.end == exp_protein_end

    @pytest.mark.skip('Run manually to regenerate the VEP REST API responses')
    @pytest.mark.parametrize(
        'contig_name, start, end, ref, alt, chlen',
        [
            ('1', 156_114_920, 156_114_920, 'G', 'A', 0,),  # `1_156114921_156114920_G_A`
            ('1', 156_115_214, 156_115_218, 'CGCC', 'C', -3,),  # `1_156114215_156114218_CGCC_C`
            ('16', 89_284_128, 89_284_130, 'CT', 'C', -1, ),  # `16_89284129_89284130_CT_C`
            ('16', 89_284_128, 89_284_134, 'CTTTTT', 'C', -5, ),  # `16_89284129_89284134_CTTTTT_C`
            ('16', 89_279_134, 89_279_135, 'G', 'C', 0, ),  # `16_89279135_89279135_G_C`
            ('X',  31_180_436, 31_180_437, 'C', 'T', 0, ),  # `X_31180437_31180437_C_T`
        ]
    )
    def test_fetch_response(
            self,
            fpath_vep_response_dir: str,
            variant_annotator: VepFunctionalAnnotator,
            genome_build: GenomeBuild,
            contig_name: str, start: int, end: int, ref: str, alt: str, chlen: int,
    ):
        contig = genome_build.contig_by_name(contig_name)
        region = GenomicRegion(contig, start, end, Strand.POSITIVE)
        vc = VariantCoordinates(region, ref, alt, chlen)

        response = variant_annotator.fetch_response(vc)
        fpath_response = os.path.join(fpath_vep_response_dir, f'{vc.variant_key}.json')
        with open(fpath_response, 'w') as fh:
            json.dump(response, fh, indent=2)


def load_json_response(fpath_json: str) -> typing.Mapping[str, typing.Any]:
    with open(fpath_json) as fh:
        return json.load(fh)
