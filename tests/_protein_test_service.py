import typing
from typing import Sequence

from gpsea.model import ProteinMetadata, ProteinFeature, FeatureInfo, FeatureType
from gpsea.model.genome import Region
from gpsea.preprocessing import ProteinMetadataService


def create_protein_test_data() -> Sequence[ProteinMetadata]:
    prots = []
    prot_feat_1 = ProteinFeature.create(FeatureInfo('domain', Region(1, 75)), FeatureType.DOMAIN)
    prot_feat_2 = ProteinFeature.create(FeatureInfo('region', Region(50, 100)), FeatureType.REGION)
    prot_doc = ProteinMetadata('NP_09876.5', 'FakeProtein', [prot_feat_1, prot_feat_2])
    prot_test = ProteinMetadata(
        protein_id='NP_037407.4', label='Ankyrin repeat domain-containing protein 11',
        protein_features=(
            ProteinFeature.create(
                feature_type=FeatureType.REPEAT,
                info=FeatureInfo('ANK 1', Region(start=167, end=196)),
            ),
            ProteinFeature.create(
                feature_type=FeatureType.REPEAT,
                info=FeatureInfo('ANK 2', Region(start=200, end=229)),
            ),
            ProteinFeature.create(
                feature_type=FeatureType.REPEAT,
                info=FeatureInfo('ANK 3', Region(start=233, end=262)),
                ),
            ProteinFeature.create(
                feature_type=FeatureType.REPEAT,
                info=FeatureInfo('ANK 4', Region(start=266, end=292)),
                ),
            ProteinFeature.create(
                feature_type=FeatureType.REGION,
                info=FeatureInfo('Disordered', Region(start=1, end=90)),
                ),
            ProteinFeature.create(
                feature_type=FeatureType.REGION,
                info=FeatureInfo('Disordered', Region(start=128, end=169)),
            ),
            ProteinFeature.create(
                feature_type=FeatureType.REGION,
                info=FeatureInfo('Disordered', Region(start=289, end=380)),
            ),
            ProteinFeature.create(
                feature_type=FeatureType.REGION,
                info=FeatureInfo('Disordered', Region(start=398, end=647)),
            ),
            ProteinFeature.create(
                feature_type=FeatureType.REGION,
                info=FeatureInfo('Disordered', Region(start=723, end=783)),
            ),
            ProteinFeature.create(
                feature_type=FeatureType.REGION,
                info=FeatureInfo('Disordered', Region(start=881, end=1043)),
            ),
            ProteinFeature.create(
                feature_type=FeatureType.REGION,
                info=FeatureInfo('Disordered', Region(start=1059, end=1393)),
            ),
            ProteinFeature.create(
                feature_type=FeatureType.REGION,
                info=FeatureInfo('Disordered', Region(start=1424, end=1710)),
            ),
            ProteinFeature.create(
                feature_type=FeatureType.REGION,
                info=FeatureInfo('Disordered', Region(start=1814, end=1836)),
            ),
            ProteinFeature.create(
                feature_type=FeatureType.REGION,
                info=FeatureInfo('Disordered', Region(start=1988, end=2019)),
            ),
            ProteinFeature.create(
                feature_type=FeatureType.REGION,
                info=FeatureInfo('Disordered', Region(start=2131, end=2406)),
            ),
            ProteinFeature.create(
                feature_type=FeatureType.REGION,
                info=FeatureInfo('Important for protein degradation', Region(start=2369, end=2663)),
            ))
    )
    prots.append(prot_test)
    prots.append(prot_doc)
    return prots


class ProteinTestMetadataService(ProteinMetadataService):

    @staticmethod
    def create():
        return ProteinTestMetadataService(create_protein_test_data())

    def __init__(
            self,
            known_protein_metadata: typing.Iterable[ProteinMetadata],
    ):
        self._known = tuple(known_protein_metadata)

    def annotate(self, protein_id: str) -> ProteinMetadata:
        for metadatum in self._known:
            if metadatum.protein_id == protein_id:
                return metadatum
        raise ValueError(f'Could not find metadata for {protein_id}')
