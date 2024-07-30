import json
import os
import typing

import pytest

from genophenocorr.io import GenophenocorrJSONEncoder
from genophenocorr.model.genome import GenomeBuild, Strand
from genophenocorr.preprocessing import VVMultiCoordinateService

from genophenocorr.model.genome import transpose_coordinate


@pytest.fixture(scope='module')
def fpath_vv_response_dir(
    fpath_preprocessing_data_dir: str,
) -> str:
    return os.path.join(fpath_preprocessing_data_dir, 'vv_response')


@pytest.fixture(scope='module')
def vv_multi_coordinate_service(
        genome_build: GenomeBuild,
) -> VVMultiCoordinateService:
    return VVMultiCoordinateService(genome_build)


class TestVVTranscriptCoordinateServiceProcessResponse:
    """
    Test processing of a JSON file with a cached response of the variant validator for a set of transcripts.
    """

    @pytest.mark.online
    @pytest.mark.parametrize(
        'tx_id, contig, start, end, strand',
        [
            # ANKRD11
            ('NM_013275.6', '16', 847_784, 1_070_716, Strand.NEGATIVE),
            # MAPK8IP3
            ('NM_001318852.2', '16', 1_706_194, 1_770_351, Strand.POSITIVE),
            # SUOX
            ('NM_001032386.2', '12', 55_997_275, 56_005_525, Strand.POSITIVE),
        ]
    )
    def test_fetch__end_to_end(
            self,
            vv_multi_coordinate_service: VVMultiCoordinateService,
            tx_id: str,
            contig: str, start: int, end: int, strand: Strand,
    ):
        tx_coordinates = vv_multi_coordinate_service.fetch(tx_id)

        assert tx_coordinates.identifier == tx_id
        assert tx_coordinates.region.contig.name == contig
        assert tx_coordinates.region.start == start
        assert tx_coordinates.region.end == end
        assert tx_coordinates.region.strand == strand

    @pytest.mark.skip('Run to regenerate the response JSON files if the response format changed')
    @pytest.mark.parametrize(
        'tx_id',
        [
            'NM_013275.6',  # ANKRD11
            'NM_001318852.2',  # MAPK8IP3
            'NM_001032386.2',  # SUOX
        ]
    )
    def test_get_response(
            self,
            tx_id: str,
            vv_multi_coordinate_service: VVMultiCoordinateService,
            fpath_vv_response_dir: str,
    ):
        response = vv_multi_coordinate_service.get_response(tx_id)
        fpath_vv_response_json = os.path.join(fpath_vv_response_dir, f'{tx_id}.json')
        with open(fpath_vv_response_json, 'w') as fh:
            json.dump(
                response, fh,
                indent=2,
            )

    @pytest.mark.parametrize(
        'tx_id, contig, start, end, strand',
        [
            # ANKRD11
            ('NM_013275.6', '16', 847_784, 1_070_716, Strand.NEGATIVE),
            # MAPK8IP3
            ('NM_001318852.2', '16', 1_706_194, 1_770_351, Strand.POSITIVE),
            # SUOX
            ('NM_001032386.2', '12', 55_997_275, 56_005_525, Strand.POSITIVE),
        ]
    )
    def test_parse_response(
            self,
            fpath_vv_response_dir: str,
            vv_multi_coordinate_service: VVMultiCoordinateService,
            tx_id: str,
            contig: str, start: int, end: int, strand: Strand,
    ):
        fpath_vv_response_json = os.path.join(fpath_vv_response_dir, f'{tx_id}.json')
        with open(fpath_vv_response_json) as fh:
            payload = json.load(fh)

        tx_coordinates = vv_multi_coordinate_service.parse_response(tx_id, payload)

        assert tx_coordinates.identifier == tx_id
        assert tx_coordinates.region.contig.name == contig
        assert tx_coordinates.region.start == start
        assert tx_coordinates.region.end == end
        assert tx_coordinates.region.strand == strand


@pytest.mark.online
class TestVVTranscriptCoordinateServiceGeneralProperties:
    """
    Test class to understand and document the structure of the JSON returned by VariantValidator
    queries for transcript coordinates.

    We run these tests to check that the REST API format did not change.
    """

    @pytest.fixture(scope='class')
    def response(
            self,
            vv_multi_coordinate_service: VVMultiCoordinateService,
    ) -> typing.Mapping[str, typing.Any]:
        # We test general properties of the response for a query for a transcript of *SUOX* gene
        tx_id = 'NM_001032386.2'
        return vv_multi_coordinate_service.get_response(tx_id)

    def test_keys(
            self,
            response: typing.Mapping[str, typing.Any],
    ):
        expected_keys = {"current_name", "current_symbol", "hgnc", "previous_symbol", "requested_symbol", "transcripts"}
        assert len(response.keys()) == len(expected_keys)
        assert all(key in expected_keys for key in response)

    def test_current_name(
            self,
            response: typing.Mapping[str, typing.Any],
    ):
        assert response.get("current_name") == "sulfite oxidase"

    def test_current_symbol(
            self,
            response: typing.Mapping[str, typing.Any],
    ):
        assert response.get("current_symbol") == "SUOX"

    def test_hgnc(
            self,
            response: typing.Mapping[str, typing.Any],
    ):
        expected = "HGNC:11460"
        assert response.get("hgnc") == expected

    def test_previous_symbol(
            self,
            response: typing.Mapping[str, typing.Any],
    ):
        assert response.get("previous_symbol") == ""

    def test_requested_symbol(
            self,
            response: typing.Mapping[str, typing.Any],
    ):
        # This is the transcript ID we requested from variantvalidator
        assert response.get("requested_symbol") == "NM_001032386.2"

    def test_transcripts(
            self,
            response: typing.Mapping[str, typing.Any],
    ):
        transcripts = response.get("transcripts")

        assert isinstance(transcripts, list)
        assert len(transcripts) == 6

    def test_transcript_zero(
            self,
            response: typing.Mapping[str, typing.Any],
    ):
        transcripts = response.get("transcripts")
        t0 = transcripts[0]
        assert isinstance(t0, dict)
        expected_keys = {'annotations', 'coding_end', 'coding_start', 'description', 'genomic_spans', 'length',
                         'reference', 'translation'}
        assert len(t0.keys()) == len(expected_keys)
        assert all(key in expected_keys for key in t0.keys())

        assert t0.get("coding_start") == 75
        assert t0.get("coding_end") == 1712
        assert t0.get("length") == 2210
        assert t0.get("reference") == "NM_001032387.2"
        assert t0.get("translation") == "NP_001027559.1"

        genomic_spans = t0.get("genomic_spans")
        assert isinstance(genomic_spans, dict)
        # Homo sapiens chromosome 12, GRCh37.p13
        assert "NC_000012.11" in genomic_spans
        # Homo sapiens chromosome 12, GRCh38.p14
        assert "NC_000012.12" in genomic_spans

        hg38span = genomic_spans.get("NC_000012.12")
        assert isinstance(hg38span.get("exon_structure"), list)
        assert hg38span.get("end_position") == 56005525
        assert len(hg38span.get("exon_structure")) == 4
        assert hg38span.get("total_exons") == 4




class TestVVTranscriptCoordinateServiceOffline:
    """
    Test processing of a cached JSON response to check if we parse the response well.
    """

    def test_ptpn11(
            self, 
            vv_multi_coordinate_service: VVMultiCoordinateService,
            fpath_vv_response_dir: str,
        ):
        tx_id = 'NM_002834.5'

        response_fpath = os.path.join(fpath_vv_response_dir, 'txid-NM_002834.5-PTPN11.json')
        response = load_response_json(response_fpath)

        tc = vv_multi_coordinate_service.parse_response(tx_id, response)

        assert tc.identifier == tx_id

        tx_region = tc.region
        assert tx_region.contig.name == '12'
        assert tx_region.start == 112_418_946
        assert tx_region.end == 112_509_918
        assert tx_region.strand == Strand.POSITIVE

        exons = tc.exons
        assert len(exons) == 16
        first = exons[0]
        assert first.start == tx_region.start
        assert first.end == 112_419_125
        last = exons[-1]
        assert last.start == 112_505_824
        assert last.end == tx_region.end
        assert all(exon.strand == tx_region.strand for exon in exons)

        assert tc.cds_start == 112_419_111
        assert tc.cds_end == 112_504_764

    def test_hbb(
        self, 
        vv_multi_coordinate_service: VVMultiCoordinateService,
        fpath_vv_response_dir: str,
    ):
        tx_id = 'NM_000518.4'

        response_fpath = os.path.join(fpath_vv_response_dir, 'txid-NM_000518.4-HBB.json')
        response = load_response_json(response_fpath)

        tc = vv_multi_coordinate_service.parse_response(tx_id, response)

        assert tc.identifier == tx_id

        tx_region = tc.region
        assert tx_region.contig.name == '11'
        assert tx_region.start_on_strand(Strand.POSITIVE) == 5_225_465
        assert tx_region.end_on_strand(Strand.POSITIVE) == 5_227_071
        assert tx_region.strand == Strand.NEGATIVE

        exons = tc.exons
        assert len(exons) == 3
        first = exons[0]
        assert first.start == tx_region.start
        assert first.start_on_strand(Strand.POSITIVE) == 5_226_929
        assert first.end_on_strand(Strand.POSITIVE) == 5_227_071


        last = exons[-1]
        assert last.end == tx_region.end
        assert last.start_on_strand(Strand.POSITIVE) == 5_225_465
        assert last.end_on_strand(Strand.POSITIVE) == 52_25_726
        assert all(exon.strand == tx_region.strand for exon in exons)

        cds_start = tc.cds_start
        assert cds_start is not None
        assert transpose_coordinate(tc.region.contig, cds_start) == 5_227_021
        cds_end = tc.cds_end
        assert cds_end is not None
        assert transpose_coordinate(tc.region.contig, cds_end) == 5_225_597
        assert tc.cds_start == 129_859_601
        assert tc.cds_end == 129_861_025

    @pytest.mark.skip('Online tests are disabled by default')
    def test_fetch_ptpn11(
        self, 
        vv_multi_coordinate_service: VVMultiCoordinateService,
    ):
        tx_id = 'NM_002834.5'

        tc = vv_multi_coordinate_service.fetch(tx_id)

        assert tc.identifier == tx_id

        tx_region = tc.region
        assert tx_region.contig.name == '12'
        assert tx_region.start == 112_418_946
        assert tx_region.end == 112_509_918
        assert tx_region.strand == Strand.POSITIVE

        exons = tc.exons
        assert len(exons) == 16
        first = exons[0]
        assert first.start == tx_region.start
        assert first.end == 112_419_125
        last = exons[-1]
        assert last.start == 112_505_824
        assert last.end == tx_region.end
        assert all(exon.strand == tx_region.strand for exon in exons)

        assert tc.cds_start == 112_419_111
        assert tc.cds_end == 112_504_764

    @pytest.mark.skip('Run to manually regenerate SUOX MANE transcript JSON dump')
    def test_fetch_suox(
        self,
        vv_multi_coordinate_service: VVMultiCoordinateService,
        fpath_suox_tx_coordinates: str,
    ):
        tx_id = 'NM_001032386.2'
        response = vv_multi_coordinate_service.fetch(tx_id)
        with open(fpath_suox_tx_coordinates, 'w') as fh:
            json.dump(response, fh, cls=GenophenocorrJSONEncoder, indent=2)


class TestVVMultiCoordinateService_as_GeneCoordinateService:

    def test_parse_hbb(
        self,
        vv_multi_coordinate_service: VVMultiCoordinateService,
        fpath_vv_response_dir: str,
    ):
        response_fpath = os.path.join(fpath_vv_response_dir, 'gene-HBB.json')
        response = load_response_json(response_fpath)

        txs = vv_multi_coordinate_service.parse_multiple(response=response)
        
        assert len(txs) == 2
        
        tx_ids = [txc.identifier for txc in txs]
        assert tx_ids == ['NM_000518.5', 'NM_000518.4']
        
        preferred = [txc.is_preferred for txc in txs]
        assert preferred == [True, False]

    def test_parse_ptpn11(
        self,
        vv_multi_coordinate_service: VVMultiCoordinateService,
        fpath_vv_response_dir: str,
    ):
        response_fpath = os.path.join(fpath_vv_response_dir, 'gene-PTPN11.json')
        response = load_response_json(response_fpath)

        txs = vv_multi_coordinate_service.parse_multiple(response=response)
        
        assert len(txs) == 9
        
        tx_ids = [txc.identifier for txc in txs]
        assert tx_ids == [
            'NM_001374625.1', 
            'NM_001330437.2', 
            'NM_080601.3', 
            'NM_002834.5', 
            'NM_002834.3', 
            'NM_001330437.1', 
            'NM_080601.2', 
            'NM_080601.1', 
            'NM_002834.4',
        ]
        
        preferred = [txc.is_preferred for txc in txs]
        assert preferred == [False, False, False, True, False, False, False, False, False]


def load_response_json(path: str):
    with open(path) as fh:
        return json.load(fh)
