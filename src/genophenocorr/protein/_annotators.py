import abc
import logging
import os
import pickle
import typing
import requests

from ._protein_data import FeatureInfo, ProteinMetadata, SimpleProteinFeature, FeatureType


class ProteinMetadataService(metaclass=abc.ABCMeta):

    @abc.abstractmethod
    def annotate(self, protein_id: str) -> typing.Sequence[ProteinMetadata]:
        """
        Get metadata for given protein ID (e.g. NP_037407.4).

        The function returns an empty sequence if there is no info for the protein ID.
        """
        pass


class UniprotProteinMetadataService(ProteinMetadataService):

    def __init__(self):
        self._logger = logging.getLogger(__name__)
        self._url = 'https://rest.uniprot.org/uniprotkb/search?query=%s&fields=accession,id,' \
                    'gene_names,gene_primary,protein_name,ft_domain,ft_motif,ft_region,ft_repeat,xref_refseq'

    def annotate(self, protein_id: str) -> typing.Sequence[ProteinMetadata]:
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
                    feat = SimpleProteinFeature(FeatureInfo(feat_name, feat_start, feat_end), FeatureType[feat_type.upper()])
                    all_features_list.append(feat)
            except KeyError:
                self._logger.warning(f"No features for {protein_id}")
            protein_list.append(ProteinMetadata(protein_id, protein_name, all_features_list))

        # TODO: DD would like to discuss an example when there are >1 items in this list.
        return protein_list

class ProteinAnnotationCache:

    def get_annotations(self, protein_id: str, directory) -> typing.Optional[typing.Sequence[ProteinMetadata]]:
        if os.path.exists(f'{directory}/{protein_id}.pickle'):
            with open(f'{directory}/{protein_id}.pickle', 'rb') as f:
                protein = pickle.load(f)
        else:
            protein = None
        return protein

    def store_annotations(self, protein_id: str, annotation: typing.Sequence[ProteinMetadata], directory):
        with open(f'{directory}/{protein_id}.pickle', 'wb') as f:
            pickle.dump(annotation, f)
        return None


class CachingFunctionalAnnotator(ProteinMetadataService):

    def __init__(self, output_directory: str, cache: ProteinAnnotationCache, fallback: ProteinMetadataService):
        if not os.path.exists(output_directory):
            raise ValueError(f"Invalid path: {output_directory}")
        directory = os.path.join(output_directory, 'annotations')
        os.makedirs(directory, exist_ok=True)
        self._output = directory
        if not isinstance(cache, ProteinAnnotationCache):
            raise ValueError(f"cache must be type VariantAnnotationCache but was type {type(cache)}")
        self._cache = cache
        if not isinstance(fallback, ProteinMetadataService):
            raise ValueError(f"fallback must be type FunctionalAnnotator but was type {type(fallback)}")
        self._fallback = fallback

    def annotate(self, protein_id: str) -> typing.Sequence[ProteinMetadata]:
        if not isinstance(protein_id, str):
            raise ValueError(f"protein_id must be type VariantCoordinates but was type {type(protein_id)}")
        annotations = self._cache.get_annotations(protein_id, self._output)
        if annotations is not None:
            # we had cached annotations
            return annotations
        else:
            ann = self._fallback.annotate(protein_id)
            self._cache.store_annotations(protein_id, ann, self._output)
            return ann
