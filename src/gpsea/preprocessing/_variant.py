import os
import pickle
import typing

import hpotk

from gpsea.model import VariantCoordinates, TranscriptAnnotation, Variant, TranscriptInfoAware, \
    TranscriptCoordinates
from ._api import FunctionalAnnotator, TranscriptCoordinateService


class VariantAnnotationCache:
    """A class that stores or retrieves Variant objects using pickle format

    Methods:
        get_annotations(variant_coordinates:VariantCoordinates): Searches a given data directory for a pickle file with variant coordinates and returns a Variant object
        store_annotations(variant_coordinates:VariantCoordinates, annotation:Variant): Creates a pickle file with variant coordinates and stores the given Variant object into that file  
    """

    def __init__(self, datadir: str):
        """Constructs all necessary attributes for a VariantAnnotationCache object

        Args:
            datadir (string): A string that references an existing directory that does or will contain all pickle files being stored
        """
        if not os.path.isdir(datadir):
            raise ValueError(f'datadir {datadir} must be an existing directory')
        self._datadir = datadir

    def get_annotations(self, variant_coordinates: VariantCoordinates) -> typing.Optional[typing.Sequence[TranscriptAnnotation]]:
        """Searches a given data directory for a pickle file with given variant coordinates and returns Variant from file. Returns None if no file is found.

        Args:
            variant_coordinates (VariantCoordinates): The variant_coordinates associated with the desired Variant
        """
        fpath = self._create_file_name(variant_coordinates)
        if os.path.isfile(fpath):
            with open(fpath, 'rb') as fh:
                return pickle.load(fh)
        else:
            return None

    def store_annotations(self, variant_coordinates: VariantCoordinates, annotations: typing.Sequence[TranscriptAnnotation]):
        """Creates a pickle file with the given variant coordinates in the file name. Loads the Variant object given into the file for storage.

        Args:
            variant_coordinates (VariantCoordinates): The variant_coordinates associated with the desired Variant
            annotations (typing.Sequence[TranscriptAnnotation]): Annotations that will be stored under the given variant coordinates
        """
        fpath = self._create_file_name(variant_coordinates)
        with open(fpath, 'wb') as f:
            pickle.dump(annotations, f)

    def _create_file_name(self, vc: VariantCoordinates) -> str:
        """Creates a file name with full location and the variant coordinates (e.g. "/path/to/desired/directory/1_2345_G_C_heterozygous.pickle")

        Args:
            vc (VariantCoordinates): The variant coordinates associated with the Variant
        """
        vk = vc.variant_key
        if len(vk) <= 50:
            fname = f'{vk}.pickle'
        else:
            # long INDELs in sequence notation
            fname = f'{vc.chrom}_{vc.start}_{vc.end}_{vc.variant_class}.pickle'
        return os.path.join(self._datadir, fname)


class VarCachingFunctionalAnnotator(FunctionalAnnotator):
    """A class that retrieves a Variant object if it exists or will run the fallback Fuctional Annotator if it does not exist.

    Methods:
        annotate(variant_coordinates:VariantCoordinates): Gets data and returns a Variant object for given variant coordinates
    """

    @staticmethod
    def with_cache_folder(
        fpath_cache_dir: str, 
        fallback: FunctionalAnnotator,
    ):
        """
        Create caching functional annotator that will store the data in `fpath_cache_dir` 
        and use `fallback` to annotate the missing variants.
        """
        cache = VariantAnnotationCache(fpath_cache_dir)
        return VarCachingFunctionalAnnotator(
            cache=cache,
            fallback=fallback,
        )

    def __init__(self, cache: VariantAnnotationCache, fallback: FunctionalAnnotator):
        """Constructs all necessary attributes for a VarCachingFunctionalAnnotator object

        Args:
            cache (VariantAnnotationCache): A VariantAnnotationCache object
            fallback (FunctionalAnnotator): A FunctionalAnnotator object
        """
        self._cache = hpotk.util.validate_instance(cache, VariantAnnotationCache, 'cache')
        self._fallback = hpotk.util.validate_instance(fallback, FunctionalAnnotator, 'fallback')

    def annotate(self, variant_coordinates: VariantCoordinates) -> typing.Sequence[TranscriptAnnotation]:
        """Gets Variant for given variant coordinates

        Args:
            variant_coordinates (VariantCoordinates): A VariantCoordinates object
        Returns:
            Variant: A Variant object
        """
        hpotk.util.validate_instance(variant_coordinates, VariantCoordinates, 'variant_coordinates')
        annotations = self._cache.get_annotations(variant_coordinates)
        if annotations is not None:
            # we had cached annotations
            return annotations
        else:
            ann = self._fallback.annotate(variant_coordinates)
            self._cache.store_annotations(variant_coordinates, ann)
            return ann

