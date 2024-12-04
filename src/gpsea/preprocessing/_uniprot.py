import logging
import requests
import typing

from gpsea.model import FeatureInfo, ProteinFeature, ProteinMetadata
from gpsea.model.genome import Region

from ._api import ProteinMetadataService


class UniprotProteinMetadataService(ProteinMetadataService):
    """A class that creates ProteinMetadata objects from data found with the Uniprot REST API.
    More info on the Uniprot REST API are
    in the `Programmatic access <https://www.uniprot.org/help/programmatic_access>`_ section.
    """

    def __init__(
        self,
        timeout: float = 30.,
    ):
        self._logger = logging.getLogger(__name__)
        self._headers = {"Content-type": "application/json"}
        
        # Query is parameterized during lookup
        query = "(%s)AND(reviewed:true)&(organism_id:9606)"
        
        # An opinionated list of fields of interest.
        # Reference: https://rest.uniprot.org/configure/uniprotkb/result-fields
        fields = ",".join(
            (
                "accession",
                "id",
                "gene_names",
                "gene_primary",
                "protein_name",
                "xref_refseq",
                "length",
                # Family & Domains
                "ft_coiled",
                "ft_compbias",
                "ft_domain",
                "ft_motif",
                "ft_region",
                "ft_repeat",
                "ft_zn_fing",
                # Sequences
                "ft_non_std",
                # Function
                "ft_act_site",
                "ft_binding",
                "cc_cofactor",
                "ft_dna_bind",
            )
        )
        self._url = f"https://rest.uniprot.org/uniprotkb/search?query={query}&fields={fields}"
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
        elif len(results) == 1:
            # In case we get only one result, let's use it!
            return UniprotProteinMetadataService._extract_metadata(
                protein_id=protein_id,
                data=results[0],
            )
        else:
            for protein in results:
                if any(uni['id'] == protein_id for uni in protein['uniProtKBCrossReferences']):
                    return UniprotProteinMetadataService._extract_metadata(
                        protein_id=protein_id,
                        data=protein,
                    )

        raise ValueError(f'Could not find an entry for {protein_id} in Uniprot response')

    @staticmethod
    def _extract_metadata(
        protein_id: str,
        data: typing.Mapping[str, typing.Any],
    ) -> ProteinMetadata:
        # `protein` has a cross-reference to the `protein_id` of interest
        try:
            protein_name = data['proteinDescription']['recommendedName']['fullName']['value']
        except KeyError:
            protein_name = data['proteinDescription']['submissionNames'][0]['fullName']['value']

        all_features_list = []
        try:
            for feature in data["features"]:
                feature_start = int(feature["location"]["start"]["value"])
                feature_end = int(feature["location"]["end"]["value"])
                feature_name = feature["description"]
                feature_info = FeatureInfo(
                    feature_name,
                    Region(start=feature_start, end=feature_end),
                )
                feature_type = feature["type"]
                protein_feature = ProteinFeature.create(feature_info, feature_type)
                all_features_list.append(protein_feature)
        except KeyError:
            raise ValueError(f"No features for {protein_id}")

        try:
            protein_length = data["sequence"]["length"]
        except KeyError as e:
            raise ValueError(e)

        return ProteinMetadata(
            protein_id,
            protein_name,
            all_features_list,
            protein_length,
        )

    def _fetch_uniprot_response(
        self,
        protein_id: str,
    ) -> typing.Mapping[str, typing.Any]:
        api_url = self._url % protein_id
        return requests.get(
            api_url,
            headers=self._headers,
            timeout=self._timeout,
        ).json()

    def annotate(self, protein_id: str) -> ProteinMetadata:
        """
        Get metadata for given protein ID.
        This class specifically only works with a RefSeq database ID (e.g. NP_037407.4).

        Args:
            protein_id (str): A protein ID
        Returns:
            ProteinMetadata: A :class:`~gpsea.model.ProteinMetadata` corresponding to the input `protein_id`.
        Raises:
            ValueError: in case of issues with `protein_id`, I/O issues, or parsing the REST response.
        """
        if not isinstance(protein_id, str):
            raise ValueError(f"Protein ID must be a str but it was {type(protein_id)}")
        if protein_id.startswith(" ") or protein_id.endswith(" "):
            raise ValueError(
                f"Please remove whitespace from protein id: `{protein_id}` and try again!"
            )
        if not protein_id.startswith("NP_"):
            raise ValueError(
                f"only works with a RefSeq database ID (e.g. NP_037407.4), but we got {protein_id}"
            )

        response = self._fetch_uniprot_response(protein_id)
        return UniprotProteinMetadataService.parse_uniprot_json(response, protein_id)
