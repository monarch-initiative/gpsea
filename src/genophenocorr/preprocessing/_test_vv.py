import json
import os

from pkg_resources import resource_filename

import pytest


from genophenocorr.model.genome import GRCh38

from ._vv import VVHgvsVariantCoordinateFinder


@pytest.fixture
def coordinate_finder() -> VVHgvsVariantCoordinateFinder:
    return VVHgvsVariantCoordinateFinder(genome_build=GRCh38)


class TestVVHgvsVariantCoordinateFinder:
    TEST_DATA_DIR = resource_filename(__name__, os.path.join('test_data', 'vv_response'))

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


def load_response_json(path: str):
    with open(path) as fh:
        return json.load(fh)
