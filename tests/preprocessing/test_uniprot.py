import json
import os
import typing

import pytest

from gpsea.io import GpseaJSONEncoder
from gpsea.preprocessing import UniprotProteinMetadataService


@pytest.fixture
def fpath_uniprot_response_dir(
    fpath_preprocessing_data_dir: str,
) -> str:
    return os.path.join(fpath_preprocessing_data_dir, "uniprot_response")


def read_json_payload(path: str) -> typing.Mapping[str, typing.Any]:
    with open(path) as fh:
        return json.load(fh)


class TestUniprotProteinMetadataService:

    @pytest.fixture
    def uniprot_metadata_service(self) -> UniprotProteinMetadataService:
        return UniprotProteinMetadataService()

    @pytest.fixture
    def zn462_human_uniprot_json(
        self,
        fpath_uniprot_response_dir: str,
    ) -> typing.Mapping[str, typing.Any]:
        fpath_zn462_human = os.path.join(fpath_uniprot_response_dir, "ZN462_HUMAN.json")
        return read_json_payload(fpath_zn462_human)

    def test_zn462(
        self,
        zn462_human_uniprot_json: typing.Mapping[str, typing.Any],
    ):
        protein_id = "NP_037407.4"
        protein_metadata = UniprotProteinMetadataService.parse_uniprot_json(
            payload=zn462_human_uniprot_json,
            protein_id=protein_id,
        )

        assert protein_metadata.protein_id == protein_id
        assert protein_metadata.label == "Ankyrin repeat domain-containing protein 11"
        assert len(protein_metadata.protein_features) == 16
        assert protein_metadata.protein_length == 2663
    
    def test_itpr1(
        self,
        fpath_uniprot_response_dir: str,
    ):
        response_json_path = os.path.join(fpath_uniprot_response_dir, 'ITPR1_HUMAN.json')
        response_json = read_json_payload(response_json_path)

        protein_id = "NP_001365381.1"
        protein_metadata = UniprotProteinMetadataService.parse_uniprot_json(
            payload=response_json,
            protein_id=protein_id,
        )

        assert protein_metadata.protein_id == protein_id
        assert protein_metadata.label == "Inositol 1,4,5-trisphosphate-gated calcium channel ITPR1"
        assert len(protein_metadata.protein_features) == 13
        assert protein_metadata.protein_length == 2758

    @pytest.mark.skip("Just for debugging")
    def test_fetch_annotate_online(
        self,
        uniprot_metadata_service: UniprotProteinMetadataService,
    ):
        query = "NP_852259.1"
        protein_meta = uniprot_metadata_service.annotate(protein_id=query)
        print(protein_meta)

    @pytest.mark.skip("Run manually to regenerate SUOX `NP_001027558.1` metadata")
    def test_fetch_suox_protein_metadata(
        self,
        uniprot_metadata_service: UniprotProteinMetadataService,
        fpath_suox_protein_metadata: str,
    ):
        suox_mane_tx_protein_id = "NP_001027558.1"
        metadata = uniprot_metadata_service.annotate(suox_mane_tx_protein_id)

        with open(fpath_suox_protein_metadata, "w") as fh:
            json.dump(
                metadata,
                fh,
                cls=GpseaJSONEncoder,
                indent=2,
            )
