import json
import os

import pytest


from genophenocorr.io import GenophenocorrJSONEncoder
from genophenocorr.model.genome import GRCh38, Strand, transpose_coordinate

from ._vv import VVHgvsVariantCoordinateFinder, VVTranscriptCoordinateService


@pytest.fixture
def coordinate_finder() -> VVHgvsVariantCoordinateFinder:
    return VVHgvsVariantCoordinateFinder(genome_build=GRCh38)


class TestVVHgvsVariantCoordinateFinder:
    TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), 'test_data', 'vv_response')

    def test_load_hgvs_MC4R_missense(self, coordinate_finder: VVHgvsVariantCoordinateFinder):
        response_fpath = os.path.join(self.TEST_DATA_DIR, 'hgvs_MC4R_missense.json')
        response = load_response_json(response_fpath)

        vc = coordinate_finder._extract_variant_coordinates(response)

        assert vc is not None
        assert vc.chrom == '18'
        assert vc.start == 60_372_096
        assert vc.end == 60_372_097
        assert vc.ref == 'T'
        assert vc.alt == 'C'
        assert vc.change_length == 0


    def test_load_hgvs_SURF2_del(self, coordinate_finder: VVHgvsVariantCoordinateFinder):
        response_fpath = os.path.join(self.TEST_DATA_DIR, 'hgvs_SURF2_del.json')
        response = load_response_json(response_fpath)

        vc = coordinate_finder._extract_variant_coordinates(response)

        assert vc is not None
        assert vc.chrom == '9'
        assert vc.start == 133_357_809
        assert vc.end == 133_357_813
        assert vc.ref == 'TAAA'
        assert vc.alt == 'T'
        assert vc.change_length == -3

    def test_load_hgvs_SURF2_dup(self, coordinate_finder: VVHgvsVariantCoordinateFinder):
        response_fpath = os.path.join(self.TEST_DATA_DIR, 'hgvs_SURF2_dup.json')
        response = load_response_json(response_fpath)

        vc = coordinate_finder._extract_variant_coordinates(response)

        assert vc is not None
        assert vc.chrom == '9'
        assert vc.start == 133_357_809
        assert vc.end == 133_357_810
        assert vc.ref == 'T'
        assert vc.alt == 'TAAA'
        assert vc.change_length == 3

    def test_load_hgvs_MC4R_dup(self, coordinate_finder: VVHgvsVariantCoordinateFinder):
        response_fpath = os.path.join(self.TEST_DATA_DIR, 'hgvs_MC4R_dup.json')
        response = load_response_json(response_fpath)

        vc = coordinate_finder._extract_variant_coordinates(response)

        assert vc is not None
        assert vc.chrom == '18'
        assert vc.start == 60_372_345
        assert vc.end == 60_372_346
        assert vc.ref == 'C'
        assert vc.alt == 'CCAT'
        assert vc.change_length == 3

    def test_load_hgvs_SURF2_ins(self, coordinate_finder: VVHgvsVariantCoordinateFinder):
        response_fpath = os.path.join(self.TEST_DATA_DIR, 'hgvs_SURF2_ins.json')
        response = load_response_json(response_fpath)

        vc = coordinate_finder._extract_variant_coordinates(response)

        assert vc is not None
        assert vc.chrom == '9'
        assert vc.start == 133_357_809
        assert vc.end == 133_357_810
        assert vc.ref == 'T'
        assert vc.alt == 'TAGC'
        assert vc.change_length == 3

    @pytest.mark.skip('Online tests are disabled by default')
    def test_find_coordinates(self, coordinate_finder: VVHgvsVariantCoordinateFinder):
        vc = coordinate_finder.find_coordinates('NM_013275.6:c.7407C>G')
        print(vc)


class TestVVTranscriptCoordinateService:
    TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), 'test_data', 'vv_response')

    @pytest.fixture
    def tx_coordinate_service(self) -> VVTranscriptCoordinateService:
        return VVTranscriptCoordinateService(genome_build=GRCh38)

    def test_ptpn11(self, tx_coordinate_service: VVTranscriptCoordinateService):
        tx_id = 'NM_002834.5'

        response_fpath = os.path.join(self.TEST_DATA_DIR, 'txid-NM_002834.5-PTPN11.json')
        response = load_response_json(response_fpath)

        tc = tx_coordinate_service.parse_response(tx_id, response)

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

    def test_hbb(self, tx_coordinate_service: VVTranscriptCoordinateService):
        tx_id = 'NM_000518.4'

        response_fpath = os.path.join(self.TEST_DATA_DIR, 'txid-NM_000518.4-HBB.json')
        response = load_response_json(response_fpath)

        tc = tx_coordinate_service.parse_response(tx_id, response)

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
    def test_fetch_ptpn11(self, tx_coordinate_service: VVTranscriptCoordinateService):
        tx_id = 'NM_002834.5'

        tc = tx_coordinate_service.fetch(tx_id)

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
            tx_coordinate_service: VVTranscriptCoordinateService,
            fpath_suox_tx_coordinates: str,
    ):
        tx_id = 'NM_001032386.2'
        response = tx_coordinate_service.fetch(tx_id)
        with open(fpath_suox_tx_coordinates, 'w') as fh:
            json.dump(response, fh, cls=GenophenocorrJSONEncoder, indent=2)


def load_response_json(path: str):
    with open(path) as fh:
        return json.load(fh)
