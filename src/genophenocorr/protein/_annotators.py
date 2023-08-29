import abc
import logging
import os
import pickle
import typing
import requests
import hpotk

from ._protein_data import FeatureInfo, ProteinMetadata, SimpleProteinFeature, FeatureType


class ProteinMetadataService(metaclass=abc.ABCMeta):
    """A metaclass that can be used to establish a class that creates ProteinMetadata objects

    Methods:
        annotate(protein_id:str): Gets metadata and creates ProteinMetadata for given protein ID
    """

    @abc.abstractmethod
    def annotate(self, protein_id: str) -> typing.Sequence[ProteinMetadata]:
        """Get metadata for given protein ID

        Args:
            protein_id (string): A protein ID
        Returns:
            Sequence[ProteinMetadata]: A sequence of ProteinMetadata objects, or an empty sequence if no data was found 
        """
        pass


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
                        feat = SimpleProteinFeature(FeatureInfo(feat_name, feat_start, feat_end), FeatureType[feat_type.upper()])
                        all_features_list.append(feat)
                except KeyError:
                    self._logger.warning(f"No features for {protein_id}")
                protein_list.append(ProteinMetadata(protein_id, protein_name, all_features_list))
            else:
                self._logger.warning(f"ID {protein_id} did not match")
        self._logger.warning(f'Protein ID {protein_id} got {len(protein_list)} results')

        # TODO - DD would like to discuss an example when there are >1 items in this list.
        return protein_list

class ProteinAnnotationCache:
    """A class that stores or retrieves ProteinMetadata objects using pickle format

    Methods:
        get_annotations(protein_id:str): Searches a given data directory for a pickle file with given ID and returns ProteinMetadata
        store_annotations(protein_id:str, annotation:Sequence[ProteinMetadata]): Creates a pickle file with given ID and stores the given ProteinMetadata into that file  
    """
    def __init__(self, datadir:str) -> None:
        """Constructs all necessary attributes for a ProteinAnnotationCache object

        Args:
            datadir (string): A string that references an existing directory that does or will contain all pickle files being stored
        """
        if not os.path.isdir(datadir):
            raise ValueError(f'datadir {datadir} must be an existing directory')
        self._datadir = datadir

    def get_annotations(self, protein_id: str) -> typing.Optional[typing.Sequence[ProteinMetadata]]:
        """Searches a given data directory for a pickle file with given ID and returns ProteinMetadata from file. Returns None if no file is found.

        Args:
            protein_id (string): The protein_id associated with the desired ProteinMetadata
        """
        fpath = self._create_file_name(protein_id)
        if os.path.isfile(fpath):
            with open(fpath, 'rb') as fh:
                return pickle.load(fh)
        else:
            return None

    def store_annotations(self, protein_id: str, annotation: typing.Sequence[ProteinMetadata]):
        """Creates a pickle file with the given protein id in the file name. Loads the ProteinMetadata given into the file for storage.

        Args:
            protein_id (string): The protein_id associated with the ProteinMetadata
            annotation (Sequence[ProteinMetadata]): A sequence of ProteinMetadata objects that will be stored under the given protein id
        """
        fpath = self._create_file_name(protein_id)
        with open(fpath, 'wb') as f:
            pickle.dump(annotation, f)

    def _create_file_name(self, prot_id: str) -> str:
        """Creates a file name with full location and the protein id (e.g. "/path/to/desired/directory/NP_037407.4.pickle")

        Args:
            protein_id (string): The protein_id associated with the ProteinMetadata
        """
        fname = f'{prot_id}.pickle'
        return os.path.join(self._datadir, fname)

class ProtCachingFunctionalAnnotator(ProteinMetadataService):
    """A class that retrieves ProteinMetadata if it exists or will run the fallback Fuctional Annotator if it does not exist.

    Methods:
        annotate(protein_id:str): Gets metadata and returns ProteinMetadata for given protein ID
    """

    def __init__(self, cache: ProteinAnnotationCache, fallback: ProteinMetadataService):
        """Constructs all necessary attributes for a ProtCachingFunctionalAnnotator object

        Args:
            cache (ProteinAnnotationCache): A ProteinAnnotationCache object
            fallback (ProteinMetadataService): A ProteinMetadataService object
        """
        self._cache = hpotk.util.validate_instance(cache, ProteinAnnotationCache, 'cache')
        self._fallback = hpotk.util.validate_instance(fallback, ProteinMetadataService, 'fallback')

    def annotate(self, protein_id: str) -> ProteinMetadata:
        """Gets metadata for given protein ID

        Args:
            protein_id (string): A protein ID
        Returns:
            ProteinMetadata: A ProteinMetadata object
        """
        hpotk.util.validate_instance(protein_id, str, 'protein_id')
        annotations = self._cache.get_annotations(protein_id)
        if annotations is not None:
            # we had cached annotations
            return annotations
        else:
            ann = self._fallback.annotate(protein_id)
            self._cache.store_annotations(protein_id, ann)

            return ann