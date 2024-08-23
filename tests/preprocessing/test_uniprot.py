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
        with open(fpath_zn462_human) as f:
            return json.load(f)

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
