import abc
import logging
import os
import pickle
import typing
import requests
import hpotk

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
            protein_list.append(ProteinMetadata(protein_id, protein_name, tuple(all_features_list)))

        # TODO - DD would like to discuss an example when there are >1 items in this list.
        return tuple(protein_list)

class ProteinAnnotationCache:

    def __init__(self, datadir:str) -> None:
        if not os.path.isdir(datadir):
            raise ValueError(f'datadir {datadir} must be an existing directory')
        self._datadir = datadir

    def get_annotations(self, protein_id: str) -> typing.Optional[typing.Sequence[ProteinMetadata]]:
        fpath = self._create_file_name(protein_id)
        if os.path.isfile(fpath):
            with open(fpath, 'rb') as fh:
                return pickle.load(fh)
        else:
            return None

    def store_annotations(self, protein_id: str, annotation: typing.Sequence[ProteinMetadata]):
        fpath = self._create_file_name(protein_id)
        with open(fpath, 'wb') as f:
            pickle.dump(annotation, f)

    def _create_file_name(self, prot_id: str) -> str:
        fname = f'{prot_id}.pickle'
        return os.path.join(self._datadir, fname)

class ProtCachingFunctionalAnnotator(ProteinMetadataService):

    def __init__(self, cache: ProteinAnnotationCache, fallback: ProteinMetadataService):
        self._cache = hpotk.util.validate_instance(cache, ProteinAnnotationCache, 'cache')
        self._fallback = hpotk.util.validate_instance(fallback, ProteinMetadataService, 'fallback')

    def annotate(self, protein_id: str) -> ProteinMetadata:
        hpotk.util.validate_instance(protein_id, str, 'variant_coordinates')
        annotations = self._cache.get_annotations(protein_id)
        if annotations is not None:
            # we had cached annotations
            return annotations
        else:
            ann = self._fallback.annotate(protein_id)
            self._cache.store_annotations(protein_id, ann)
            return ann