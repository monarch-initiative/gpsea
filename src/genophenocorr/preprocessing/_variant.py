import logging
import os
import pickle
import typing

import hpotk
import requests


from genophenocorr.model import VariantCoordinates, TranscriptAnnotation, Variant, TranscriptInfoAware, TranscriptCoordinates

from ._api import FunctionalAnnotator, ProteinMetadataService, TranscriptCoordinateService


def verify_start_end_coordinates(vc: VariantCoordinates):
    """
    Converts the 0-based VariantCoordinates to ones that will be interpreted correctly by VEP
    
    Example - an insertion/duplication of G after the given G at coordinate 3:
    1 2 3 4 5
    A C G T A

    0-based: 2 3 G GG       1-based: 3 G GG         VEP: 4 3 - G 

    Args:
        vc (VariantCoordinates): A VariantCoordinates object
    Returns:
        string: The variant coordinates formatted to work with VEP
    """

    chrom = vc.chrom
    end = vc.end
    start = vc.start + 1
    alt = vc.alt
    if vc.is_structural():
        alt = vc.alt[1:-1]
        # TODO: Verify <INS> are working correctly
    else:
        if len(vc.ref) == 0 or len(vc.alt) == 0:
            raise ValueError(f'Trimmed alleles are not yet supported!')
        if len(vc.ref) == 1 and len(vc.alt) != 1:
            # INS/DUP
            start = start + 1  # we must "trim"
            end = end - 1
            alt = vc.alt[1:]
            # 100 AC AGT
            # MNV

    return f'{chrom}:{start}-{end}/{alt}'


class VepFunctionalAnnotator(FunctionalAnnotator):
    """A class that peforms functional annotation of variant coordinates using Variant Effect Predictor (VEP) REST API.

    Methods:
        annotate(variant_coordinates: VariantCoordinates): the variant to annotate.
    """

    def __init__(self, protein_annotator: ProteinMetadataService):
        """Constructs all necessary attributes for a VepFunctionalAnnotator object

        Args:
            protein_annotator (ProteinMetadataService): A ProteinMetadataService object for ProteinMetadata creation
        """
        self._logging = logging.getLogger(__name__)
        self._protein_annotator = protein_annotator
        self._url = 'https://rest.ensembl.org/vep/human/region/%s?LoF=1&canonical=1&domains=1&hgvs=1' \
                    '&mutfunc=1&numbers=1&protein=1&refseq=1&mane=1&transcript_version=1&variant_class=1'

    def annotate(self, variant_coordinates: VariantCoordinates) -> typing.Sequence[TranscriptAnnotation]:
        """Perform functional annotation using Variant Effect Predictor (VEP) REST API.

        Args:
            variant_coordinates (VariantCoordinates): A VariantCoordinates object
        Returns:
            typing.Sequence[TranscriptAnnotation]: A sequence of transcript annotations for the variant coordinates
        """
        variant = self._query_vep(variant_coordinates)
        # canon_tx = None
        annotations = []
        for trans in variant.get('transcript_consequences'):
            trans_id = trans.get('transcript_id')
            if not trans_id.startswith('NM'):
                continue
            # TODO - implement
            is_preferred = False
            # if trans.get('canonical') == 1:
            #     canon_tx = trans_id
            hgvsc_id = trans.get('hgvsc')
            consequences = trans.get('consequence_terms')
            gene_name = trans.get('gene_symbol')
            protein_id = trans.get('protein_id')
            protein = self._protein_annotator.annotate(protein_id)
            protein_effect_start = trans.get('protein_start')
            protein_effect_end = trans.get('protein_end')
            if protein_effect_start is None and protein_effect_end is not None:
                protein_effect_start = 1
            if protein_effect_end is not None:
                protein_effect_end = int(protein_effect_end)
            if protein_effect_start is not None:
                protein_effect_start = int(protein_effect_start)
            exons_effected = trans.get('exon')
            if exons_effected is not None:
                exons_effected = exons_effected.split('/')[0].split('-')
                if len(exons_effected) == 2:
                    exons_effected = range(int(exons_effected[0]), int(exons_effected[1]) + 1)
                exons_effected = [int(x) for x in exons_effected]
            annotations.append(
                TranscriptAnnotation(gene_name,
                                     trans_id,
                                     hgvsc_id,
                                     is_preferred,
                                     consequences,
                                     exons_effected,
                                     protein,
                                     protein_effect_start,
                                     protein_effect_end)
            )

        return annotations

    def _query_vep(self, variant_coordinates) -> dict:
        api_url = self._url % (verify_start_end_coordinates(variant_coordinates))
        r = requests.get(api_url, headers={'Content-Type': 'application/json'})
        if not r.ok:
            self._logging.error(f"Expected a result but got an Error for variant: {variant_coordinates.as_string()}")
            r.raise_for_status()
        results = r.json()
        if not isinstance(results, list):
            self._logging.error(results.get('error'))
            raise ConnectionError(f"Expected a result but got an Error. See log for details.")
        if len(results) > 1:
            self._logging.error([result.id for result in results])
            raise ValueError(f"Expected only one variant per request but received {len(results)} different variants.")
        return results[0]


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

    def _create_file_name(self, variant_coordinates: VariantCoordinates) -> str:
        """Creates a file name with full location and the variant coordinates (e.g. "/path/to/desired/directory/1_2345_G_C_heterozygous.pickle")

        Args:
            variant_coordinates (VariantCoordinates): The variant coordinates associated with the Variant
        """
        if len(variant_coordinates.as_string()) <= 50:
            fname = f'{variant_coordinates.as_string()}.pickle'
        else:
            fname = f'{variant_coordinates.chrom}_{variant_coordinates.start}_{variant_coordinates.end}_{variant_coordinates.genotype}.pickle'
        return os.path.join(self._datadir, fname)


class VarCachingFunctionalAnnotator(FunctionalAnnotator):
    """A class that retrieves a Variant object if it exists or will run the fallback Fuctional Annotator if it does not exist.

    Methods:
        annotate(variant_coordinates:VariantCoordinates): Gets data and returns a Variant object for given variant coordinates
    """

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


class EnsemblTranscriptCoordinateService(TranscriptCoordinateService):

    def fetch(self, tx: TranscriptInfoAware) -> TranscriptCoordinates:
        # TODO - implement
        pass
