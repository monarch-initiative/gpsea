import pytest
import os

from gpsea.model import ProteinMetadata, FeatureType


# Homo sapiens inositol 1,4,5-trisphosphate receptor type 1 (ITPR1), transcript variant 4, mRNA
ITPR1_protein_id = "NP_001365381.1"
ITPR1_protein_len = 2758


#P17010
ZFX_protein_len = 805
ZFX_protein_id = "NP_001171555.1"

class TestUniprotJsonToMetadata:
    """
    Test function that ingests UniProt JSON and transforms it to a ProteinMetadata object
    """

    @pytest.fixture
    def q8izt6_json_file_path(
        self,
        fpath_preprocessing_data_dir: str,
    ) -> str:
        return os.path.join(fpath_preprocessing_data_dir, "uniprot_response", "Q8IZT6_manual_download.json")

    @pytest.fixture
    def q8izt6_protein_metadata(
        self,
        q8izt6_json_file_path: str,
    ) -> ProteinMetadata:
        """
        :returns: ProteinMetadata created from a downloaded UniProt JSON file
        """
        return ProteinMetadata.from_uniprot_json(
            protein_id=ITPR1_protein_id,
            label=ITPR1_protein_id,
            uniprot_json=q8izt6_json_file_path,
            protein_length=ITPR1_protein_len,
        )

    @pytest.fixture
    def P17010_json_file_path(
            self,
            fpath_preprocessing_data_dir: str,
    ) -> str:
        return os.path.join(fpath_preprocessing_data_dir, "uniprot_response", "P17010_manual_download.json")

    @pytest.fixture
    def P17010_protein_metadata(
            self,
            P17010_json_file_path#: str,
    ) -> ProteinMetadata:
        """
        :returns: ProteinMetadata created from a downloaded UniProt JSON file
        """
        return ProteinMetadata.from_uniprot_json(
            protein_id=ZFX_protein_id,
            label=ZFX_protein_id,
            uniprot_json=P17010_json_file_path,
            protein_length=ZFX_protein_len,
        )

    def test_general_info(
        self,
        q8izt6_protein_metadata: ProteinMetadata,
    ):
        """
        grep -o location  Q8IZT6_manual_download.json | wc -l
        47
        Thus there are 47 Features in this JSON file.
        """
        assert ITPR1_protein_len == q8izt6_protein_metadata.protein_length
        assert ITPR1_protein_id == q8izt6_protein_metadata.protein_id
        assert len(q8izt6_protein_metadata.protein_features) == 47

    def test_first_feature(
        self,
        q8izt6_protein_metadata: ProteinMetadata,
    ):
        """
        By inspection, the first feature is this:
        {"type":"Domain",
            "location":
                {"start":{"value":920,"modifier":"EXACT"},
                "end":{"value":1056,"modifier":"EXACT"}},
            "description":"Calponin-homology (CH) 1",
            "evidences":[{"evidenceCode":"ECO:0000255","source":"PROSITE-ProRule","id":"PRU00044"}]},
        """
        feature_0 = q8izt6_protein_metadata.protein_features[0]
        assert feature_0.feature_type == FeatureType.DOMAIN
        # Note that we have converted to zero based internally, so this corresponds to 920 in the input
        assert feature_0.info.start == 919
        assert feature_0.info.end == 1056
        assert feature_0.info.name == "Calponin-homology (CH) 1"


    def test_ZFX(self,
                 P17010_protein_metadata):
        """
        :[{"type":"Zinc finger",
        "location":{"start":{"value":425,"modifier":"EXACT"},
        "end":{"value":447,"modifier":"EXACT"}},"
        """
        assert P17010_protein_metadata is not None
        feature_0 = P17010_protein_metadata.protein_features[0]
        assert feature_0.feature_type == FeatureType.ZINC_FINGER
        assert feature_0.info.start == 424 ## zero based open-closed
        assert feature_0.info.end == 447

