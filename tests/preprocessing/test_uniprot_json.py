import pytest
import pandas as pd
import os

from gpsea.model import ProteinMetadata, FeatureType


# Homo sapiens inositol 1,4,5-trisphosphate receptor type 1 (ITPR1), transcript variant 4, mRNA
ITPR1_protein_id = "NP_001365381.1"
ITPR1_protein_len = 2758



class TestUniprotJsonToMetadata:
    """
    Test function that ingests UniProt JSON and transforms it to a ProteinMetadata object
    """

    @pytest.fixture
    def Q8IZT6_json_file_path(
        fpath_preprocessing_data_dir: str,
    ) -> str:
        return os.path.join(fpath_preprocessing_data_dir, "Q8IZT6_manual_download.json")

    @pytest.fixture
    def q8izt6_protein_metadata(self,
                                Q8IZT6_json_file_path:str) -> ProteinMetadata:
        """
        :returns: ProteinMetadata created from a downloaded UniProt JSON file
        """
        return ProteinMetadata.from_uniprot_json(
            protein_id=ITPR1_protein_id,
            label=ITPR1_protein_id,
            uniprot_json=Q8IZT6_json_file_path,
            protein_length=ITPR1_protein_len,
        )
    
    def test_general_info(self, 
                          q8izt6_protein_metadata: ProteinMetadata):
        assert ITPR1_protein_len == q8izt6_protein_metadata.protein_length
        assert ITPR1_protein_id == q8izt6_protein_metadata.protein_id
        assert len(q8izt6_protein_metadata.protein_features) == 2
