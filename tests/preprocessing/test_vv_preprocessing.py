import json
import os
import typing

import pytest

from genophenocorr.model.genome import GenomeBuild, Strand
from genophenocorr.preprocessing import VVTranscriptCoordinateService


@pytest.fixture(scope='module')
def fpath_vv_response_dir(fpath_test_data: str) -> str:
    return os.path.join(fpath_test_data, 'vv')


@pytest.fixture(scope='module')
def coordinate_service(
        genome_build_hg38: GenomeBuild,
) -> VVTranscriptCoordinateService:
    return VVTranscriptCoordinateService(genome_build_hg38)


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
            coordinate_service: VVTranscriptCoordinateService,
            tx_id: str,
            contig: str, start: int, end: int, strand: Strand,
    ):
        tx_coordinates = coordinate_service.fetch(tx_id)

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
            coordinate_service: VVTranscriptCoordinateService,
            fpath_vv_response_dir: str,
    ):
        response = coordinate_service.get_response(tx_id)
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
            coordinate_service: VVTranscriptCoordinateService,
            tx_id: str,
            contig: str, start: int, end: int, strand: Strand,
    ):
        fpath_vv_response_json = os.path.join(fpath_vv_response_dir, f'{tx_id}.json')
        with open(fpath_vv_response_json) as fh:
            payload = json.load(fh)

        tx_coordinates = coordinate_service.parse_response(tx_id, payload)

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
    def coordinate_service(
            self,
            genome_build_hg38: GenomeBuild,
    ) -> VVTranscriptCoordinateService:
        return VVTranscriptCoordinateService(genome_build_hg38)

    @pytest.fixture(scope='class')
    def response(
            self,
            coordinate_service: VVTranscriptCoordinateService,
    ) -> typing.Mapping[str, typing.Any]:
        # We test general properties of the response for a query for a transcript of *SUOX* gene
        tx_id = 'NM_001032386.2'
        return coordinate_service.get_response(tx_id)

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
