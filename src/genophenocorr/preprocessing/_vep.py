# A module with classes that interact with Ensembl's REST API to fetch required data.
import logging
import typing

import requests

from genophenocorr.model import VariantCoordinates, TranscriptAnnotation, VariantEffect
from genophenocorr.model.genome import Region
from ._api import FunctionalAnnotator, ProteinMetadataService


def verify_start_end_coordinates(vc: VariantCoordinates):
    """
    Converts the 0-based VariantCoordinates to ones that will be interpreted
    correctly by VEP

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
            alt = vc.alt[1:]
            # 100 AC AGT
            # MNV

    return f'{chrom}:{start}-{end}/{alt}'


class VepFunctionalAnnotator(FunctionalAnnotator):
    """A `FunctionalAnnotator` that uses Variant Effect Predictor (VEP) REST API to
    do functional variant annotation.

    Args:
        protein_annotator (ProteinMetadataService): a service for getting protein data
        include_computational_txs (bool): Include computational transcripts, such as
        RefSeq `XM_`.
    """

    def __init__(self, protein_annotator: ProteinMetadataService,
                 include_computational_txs: bool = False):
        self._logger = logging.getLogger(__name__)
        self._protein_annotator = protein_annotator
        self._url = 'https://rest.ensembl.org/vep/human/region/%s?LoF=1&canonical=1' \
                    '&domains=1&hgvs=1' \
                    '&mutfunc=1&numbers=1&protein=1&refseq=1&mane=1' \
                    '&transcript_version=1&variant_class=1'
        self._include_computational_txs = include_computational_txs

    def annotate(self, variant_coordinates: VariantCoordinates) -> typing.Sequence[TranscriptAnnotation]:
        """Perform functional annotation using Variant Effect Predictor (VEP) REST API.

        Args:
            variant_coordinates (VariantCoordinates): A VariantCoordinates object
        Returns:
            typing.Sequence[TranscriptAnnotation]: A sequence of transcript
            annotations for the variant coordinates
        """
        response = self._query_vep(variant_coordinates)
        annotations = []
        if response is None:
            self._logger.error('VEP did not finish successfully.')
            return None
        if 'transcript_consequences' not in response:
            self._logger.error('The VEP response lacked the required `transcript_consequences` field. %s', response)
            return None
        for trans in response['transcript_consequences']:
            annotation = self._process_item(trans)
            if annotation is not None:
                annotations.append(annotation)

        return annotations

    def _parse_variant_effect(self, effect: str) -> typing.Optional[VariantEffect]:
        effect = effect.upper()
        if effect == "5_PRIME_UTR_VARIANT":
            effect = "FIVE_PRIME_UTR_VARIANT"
        elif effect == "3_PRIME_UTR_VARIANT":
            effect = 'THREE_PRIME_UTR_VARIANT'
        try:
            var_effect = VariantEffect[effect]
        except KeyError:
            self._logger.warning("VariantEffect %s was not found in our record of possible effects. Please report this issue to the genophenocorr GitHub.", effect)
            return None
        return var_effect

    def _process_item(self, item: typing.Dict) -> typing.Optional[TranscriptAnnotation]:
        """
        Parse one transcript annotation from the JSON response.
        """
        trans_id = item.get('transcript_id')
        if not self._include_computational_txs and not trans_id.startswith('NM_'):
            # Skipping a computational transcript
            return None
        is_preferred = True if ('canonical' in item and item['canonical'] == 1) else False
        hgvsc_id = item.get('hgvsc')
        var_effects = []
        consequences = item.get('consequence_terms')
        for con in consequences:
            var_effect = self._parse_variant_effect(con)
            if var_effect is not None:
                var_effects.append(var_effect)
        gene_name = item.get('gene_symbol')
        exons_effected = item.get('exon')
        if exons_effected is not None:
            exons_effected = exons_effected.split('/')[0].split('-')
            if len(exons_effected) == 2:
                exons_effected = range(int(exons_effected[0]),
                                       int(exons_effected[1]) + 1)
            exons_effected = (int(x) for x in exons_effected)

        protein_id = item.get('protein_id')
        protein = self._protein_annotator.annotate(protein_id)
        protein_effect_start = item.get('protein_start')
        protein_effect_end = item.get('protein_end')
        if protein_effect_start is None or protein_effect_end is None:
            # Does this ever happen? Let's log a warning for now and address the absence of a coordinate later,
            # if we see a lot of these warnings popping out.
            # Note that Lauren's version of the code had a special branch for missing start, where she set the variable
            # to `1` (1-based coordinate).
            self._logger.warning('Missing start/end coordinate for %s on protein %s', hgvsc_id, protein_id)
            protein_effect = None
        else:
            # The coordinates are in 1-based system and we need 0-based.
            protein_effect_start = int(protein_effect_start) - 1
            protein_effect_end = int(protein_effect_end)
            protein_effect = Region(protein_effect_start, protein_effect_end)


        return TranscriptAnnotation(gene_name,
                                    trans_id,
                                    hgvsc_id,
                                    is_preferred,
                                    var_effects,
                                    exons_effected,
                                    protein,
                                    protein_effect)

    def _query_vep(self, variant_coordinates: VariantCoordinates) -> dict:
        api_url = self._url % (verify_start_end_coordinates(variant_coordinates))
        r = requests.get(api_url, headers={'Content-Type': 'application/json'})
        if not r.ok:
            self._logger.error("Expected a result but got an Error for variant: %s", variant_coordinates.variant_key)
            self._logger.error(r.raise_for_status())
            return None
        results = r.json()
        if not isinstance(results, list):
            self._logger.error(results.get('error'))
            raise ConnectionError(
                f"Expected a result but got an Error. See log for details.")
        if len(results) > 1:
            self._logger.error("Expected only one variant per request but received %s different variants.", len(results))
            self._logger.error([result.id for result in results])
            raise ValueError(
                f"Expected only one variant per request but received {len(results)} "
                f"different variants.")
        return results[0]
