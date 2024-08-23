import os
import pickle
import typing
import hpotk


from gpsea.model import ProteinMetadata

from ._api import ProteinMetadataService


class ProteinAnnotationCache:
    """A class that stores or retrieves ProteinMetadata objects using pickle format

    Methods:
        get_annotations(protein_id:str): Searches a given data directory for a pickle file with given ID and returns ProteinMetadata
        store_annotations(protein_id:str, annotation:Sequence[ProteinMetadata]): Creates a pickle file with given ID and stores the given ProteinMetadata into that file
    """
    def __init__(self, datadir: str) -> None:
        """Constructs all necessary attributes for a ProteinAnnotationCache object

        Args:
            datadir (string): A string that references an existing directory that does or will contain all pickle files being stored
        """
        if not os.path.isdir(datadir):
            raise ValueError(f'datadir {datadir} must be an existing directory')
        self._datadir = datadir

    def get_annotations(self, protein_id: str) -> typing.Optional[ProteinMetadata]:
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

    def store_annotations(self, protein_id: str, annotation: ProteinMetadata):
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
            prot_id (string): The protein_id associated with the ProteinMetadata
        """
        fname = f'{prot_id}.pickle'
        return os.path.join(self._datadir, fname)


class ProtCachingMetadataService(ProteinMetadataService):
    """A class that retrieves ProteinMetadata if it exists or will run the fallback Fuctional Annotator if it does not exist.

    Methods:
        annotate(protein_id:str): Gets metadata and returns ProteinMetadata for given protein ID
    """

    def __init__(self, cache: ProteinAnnotationCache, fallback: ProteinMetadataService):
        """Constructs all necessary attributes for a ProtCachingMetadataService object

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
