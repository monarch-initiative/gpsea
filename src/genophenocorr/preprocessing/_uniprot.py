import logging
import requests
import typing

from genophenocorr.model import FeatureInfo, FeatureType, ProteinFeature, ProteinMetadata

from ._api import ProteinMetadataService


class UniprotProteinMetadataService(ProteinMetadataService):
    """A class that creates ProteinMetadata objects from data found with the Uniprot REST API.
    More info on the Uniprot REST API here - https://www.uniprot.org/help/programmatic_access

    Methods:
        annotate(protein_id:str): Gets metadata and creates ProteinMetadata for given protein ID
    """

    def __init__(self):
        """Constructs all necessary attributes for a UniprotProteinMetadataService object
        """
        self._logger = logging.getLogger(__name__)
        self._url = 'https://rest.uniprot.org/uniprotkb/search?query=(%s)AND(reviewed:true)&fields=accession,id,' \
                    'gene_names,gene_primary,protein_name,ft_domain,ft_motif,ft_region,ft_repeat,xref_refseq'

    def annotate(self, protein_id: str) -> typing.Sequence[ProteinMetadata]:
        """Get metadata for given protein ID. This class specifically only works with a RefSeq database ID (e.g. NP_037407.4)

        Args:
            protein_id (string): A protein ID
        Returns:
            Sequence[ProteinMetadata]: A sequence of ProteinMetadata objects, or an empty sequence if no data was found.
        """
        if not isinstance(protein_id, str):
            # TODO - log a warning?
            raise ValueError(f'Protein ID must be a str but it was {type(protein_id)}')
        api_url = self._url % protein_id
        r = requests.get(api_url).json()
        results = r['results']
        if len(results) == 0:
            self._logger.warning(f"No proteins found for ID {protein_id}. Please verify refseq ID.")
            return []
        protein_list = []
        for protein in results:
            verify = False
            for uni in protein['uniProtKBCrossReferences']:
                if uni['id'] == protein_id:
                    verify = True
            if verify:
                try:
                    protein_name = protein['proteinDescription']['recommendedName']['fullName']['value']
                except KeyError:
                    protein_name = protein['proteinDescription']['submissionNames'][0]['fullName']['value']
                all_features_list = []
                try:
                    for feature in protein['features']:
                        feat_type = feature['type']
                        feat_name = feature['description']
                        feat_start = int(feature['location']['start']['value'])
                        feat_end = int(feature['location']['end']['value'])
                        feat = ProteinFeature.create(FeatureInfo(feat_name, feat_start, feat_end), FeatureType[feat_type.upper()])
                        all_features_list.append(feat)
                except KeyError:
                    self._logger.warning(f"No features for {protein_id}")
                protein_list.append(ProteinMetadata(protein_id, protein_name, all_features_list))
            else:
                self._logger.warning(f"ID {protein_id} did not match")
        self._logger.warning(f'Protein ID {protein_id} got {len(protein_list)} results')

        # TODO - DD would like to discuss an example when there are >1 items in this list.
        return protein_list
