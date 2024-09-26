import pytest
import pandas as pd

from gpsea.model import ProteinMetadata, FeatureType


# Homo sapiens inositol 1,4,5-trisphosphate receptor type 1 (ITPR1), transcript variant 4, mRNA
ITPR1_protein_id = "NP_001365381.1"
ITPR1_protein_len = 2758


class TestProteinMetadata:

    @pytest.fixture
    def itpr1_protein_metadata(self) -> ProteinMetadata:
        region_list = [
            {
                "region": "Suppresor domain",
                "category": "domain",
                "start": 1,
                "end": 223,
            },
            {
                "region": "IP3 binding",
                "category": "region",
                "start": 224,
                "end": 578,
            },
        ]
        df = pd.DataFrame(region_list)
        return ProteinMetadata.from_feature_frame(
            protein_id=ITPR1_protein_id,
            label=ITPR1_protein_id,
            features=df,
            protein_length=ITPR1_protein_len,
        )

    def test_general_info(self, itpr1_protein_metadata: ProteinMetadata):
        assert ITPR1_protein_len == itpr1_protein_metadata.protein_length
        assert ITPR1_protein_id == itpr1_protein_metadata.protein_id
        assert len(itpr1_protein_metadata.protein_features) == 2

    def test_suppressor_domain(self, itpr1_protein_metadata: ProteinMetadata):
        domain = itpr1_protein_metadata.protein_features[0]

        assert domain.info.name == "Suppresor domain"
        assert domain.info.start == 0  # 0-based!
        assert domain.info.end == 223
        assert domain.feature_type == FeatureType.DOMAIN

    def test_IP3_domain(self, itpr1_protein_metadata: ProteinMetadata):
        region = itpr1_protein_metadata.protein_features[1]
        
        assert region.info.name == "IP3 binding"
        assert region.info.start == 223  # 0-based!
        assert region.info.end == 578
        assert region.feature_type == FeatureType.REGION

    def test_malformed_protein_metadata(self):
        """
        This test tries to create a ProteinFeatures that do not include FeatureType information and thus rightly fails
        """
        with pytest.raises(Exception) as e:
            regions = [
                {
                    "region": "Suppresor domain",
                    "start": 1,
                    "end": 223,
                },
                {
                    "region": "IP3 binding",
                    "start": 224,
                    "end": 578,
                },
            ]
            df = pd.DataFrame(regions)
            ProteinMetadata.from_feature_frame(
                protein_id=ITPR1_protein_id,
                label=ITPR1_protein_id,
                features=df,
                protein_length=ITPR1_protein_len,
            )
        
        assert e.value.args[0] == "The column(s) {category} are missing from the `features` DataFrame: ('region', 'start', 'end')"
        
