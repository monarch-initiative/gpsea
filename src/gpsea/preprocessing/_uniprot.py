import logging
import requests
import typing

from gpsea.model import FeatureInfo, FeatureType, ProteinFeature, ProteinMetadata
from gpsea.model.genome import Region

from ._api import ProteinMetadataService


class UniprotProteinMetadataService(ProteinMetadataService):
    """A class that creates ProteinMetadata objects from data found with the Uniprot REST API.
    More info on the Uniprot REST API here - https://www.uniprot.org/help/programmatic_access

    Methods:
        annotate(protein_id:str): Gets metadata and creates ProteinMetadata for given protein ID
    """

    def __init__(
        self,
        timeout: float = 30.,
    ):
        self._logger = logging.getLogger(__name__)
        self._url = 'https://rest.uniprot.org/uniprotkb/search?query=(%s)AND(reviewed:true)&fields=accession,id,' \
                    'gene_names,gene_primary,protein_name,ft_domain,ft_motif,ft_region,ft_repeat,xref_refseq,length'
        self._timeout = timeout

    @staticmethod
    def parse_uniprot_json(
            payload: typing.Mapping[str, typing.Any],
            protein_id: str,
    ) -> ProteinMetadata:
        """
        Try to extract `ProteinMetadata` corresponding to `protein_id` from the Uniprot JSON `payload`.

        Args:
            payload: a JSON object corresponding to Uniprot response
            protein_id: a `str` with the accession the protein of interest
        Returns:
            a complete instance of `ProteinMetadata`
        Raises:
            `ValueError` if unable to parse a complete instance
        """
        results = payload['results']
        if len(results) == 0:
            raise ValueError(f"No proteins found for ID {protein_id}. Please verify refseq ID.")

        for protein in results:
            if any(uni['id'] == protein_id for uni in protein['uniProtKBCrossReferences']):
                # `protein` has a cross-reference to the `protein_id` of interest
                try:
                    protein_name = protein['proteinDescription']['recommendedName']['fullName']['value']
                except KeyError:
                    protein_name = protein['proteinDescription']['submissionNames'][0]['fullName']['value']

                all_features_list = []
                try:
                    for feature in protein['features']:
                        feature_start = int(feature['location']['start']['value'])
                        feature_end = int(feature['location']['end']['value'])
                        feature_name = feature['description']
                        feature_info = FeatureInfo(
                            feature_name,
                            Region(start=feature_start, end=feature_end),
                        )
                        feature_type = FeatureType[feature['type'].upper()]
                        protein_feature = ProteinFeature.create(feature_info, feature_type)
                        all_features_list.append(protein_feature)
                except KeyError:
                    raise ValueError(f"No features for {protein_id}")

                try:
                    protein_length = protein["sequence"]["length"]
                except KeyError as e:
                    raise ValueError(e)

                return ProteinMetadata(protein_id, protein_name, all_features_list, protein_length)

        raise ValueError(f'Could not find an entry for {protein_id} in Uniprot response')

    def annotate(self, protein_id: str) -> ProteinMetadata:
        """
        Get metadata for given protein ID.
        This class specifically only works with a RefSeq database ID (e.g. NP_037407.4).

        Args:
            protein_id (string): A protein ID
        Returns:
            Sequence[ProteinMetadata]: A sequence of ProteinMetadata objects, or an empty sequence if no data was found.
        """
        if not isinstance(protein_id, str):
            raise ValueError(f'Protein ID must be a str but it was {type(protein_id)}')
        if protein_id.startswith(" ") or protein_id.endswith(" "):
            raise ValueError(f"Please remove whitespace from protein id: \"{protein_id}\" and try again!")
        if not protein_id.startswith("NP_"):
            raise ValueError(f"only works with a RefSeq database ID (e.g. NP_037407.4), but we got {protein_id}")
        api_url = self._url % protein_id
        response = requests.get(api_url, timeout=self._timeout).json()
        return UniprotProteinMetadataService.parse_uniprot_json(response, protein_id)
