import json
import os
import typing

import pytest
import pandas as pd

from gpsea.model import ProteinMetadata, FeatureType


ITPR1_protein_id = 'NP_001365381.1' # Homo sapiens inositol 1,4,5-trisphosphate receptor type 1 (ITPR1), transcript variant 4, mRNA
ITPR1_protein_len = 2758

@pytest.fixture
def itpr1_protein_metadata() -> ProteinMetadata:
    region_list = list()
    region_list.append({"region":"Suppresor domain", "start": 1, "end":223, "category": "domain"})
    region_list.append({"region":"IP3 binding", "start": 224, "end":578, "category": "region"})
    df = pd.DataFrame(region_list)
    return ProteinMetadata.from_feature_frame(protein_id=ITPR1_protein_id, label=ITPR1_protein_id, features=df, length=ITPR1_protein_len)




class TestProteinMetadata:
    
    def test_general_info(self, itpr1_protein_metadata: ProteinMetadata):
        assert ITPR1_protein_len ==  itpr1_protein_metadata.protein_length
        assert ITPR1_protein_id == itpr1_protein_metadata.protein_id
        features = itpr1_protein_metadata.protein_features
        assert 2 == len(features)
       

    def test_suppressor_domain(self, itpr1_protein_metadata: ProteinMetadata):
        features = itpr1_protein_metadata.protein_features
        assert 2 == len(features)
        domains = [ft for ft in features if ft.feature_type == FeatureType.DOMAIN]
        assert 1 == len(domains)
        domain = domains[0]
        assert "Suppresor domain" == domain.info.name
        assert 1 == domain.info.start
        assert 223 == domain.info.end
        assert FeatureType.DOMAIN == domain.feature_type

    def test_IP3_domain(self, itpr1_protein_metadata: ProteinMetadata):
        features = itpr1_protein_metadata.protein_features
        assert 2 == len(features)
        regions = [ft for ft in features if ft.feature_type == FeatureType.REGION]
        assert 1 == len(regions)
        region = regions[0]
        assert "IP3 binding" == region.info.name
        assert 224 == region.info.start
        assert 578 == region.info.end
        assert FeatureType.REGION == region.feature_type


    def test_malformed_protein_metadata(self):
        """
        This test tries to create a ProteinFeatures that do not include FeatureType information and thus rightly fails
        """
        with pytest.raises(Exception) as e_info:
            region_list = list()
            region_list.append({"region":"Suppresor domain", "start": 1, "end":223,})
            region_list.append({"region":"IP3 binding", "start": 224, "end":578, })
            df = pd.DataFrame(region_list)
            pmd = ProteinMetadata.from_feature_frame(protein_id=ITPR1_protein_id, label=ITPR1_protein_id, features=df, length=ITPR1_protein_len)



