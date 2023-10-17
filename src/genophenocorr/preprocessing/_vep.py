# A module with classes that interact with Ensembl's REST API to fetch required data.
import logging
import re
import typing

import hpotk
import requests

from genophenocorr.model import VariantCoordinates, TranscriptAnnotation, Genotype
from genophenocorr.model.genome import GenomeBuild
from genophenocorr.model.genome import GRCh38
from genophenocorr.model.genome import GenomicRegion, Strand

from genophenocorr.preprocessing._api import FunctionalAnnotator, \
    ProteinMetadataService, \
    VariantCoordinateFinder, T


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
            end = end - 1
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
        self._logging = logging.getLogger(__name__)
        self._protein_annotator = protein_annotator
        self._url = 'https://rest.ensembl.org/vep/human/region/%s?LoF=1&canonical=1' \
                    '&domains=1&hgvs=1' \
                    '&mutfunc=1&numbers=1&protein=1&refseq=1&mane=1' \
                    '&transcript_version=1&variant_class=1'
        self._include_computational_txs = include_computational_txs

    def annotate(self, variant_coordinates: VariantCoordinates) -> typing.Sequence[
        TranscriptAnnotation]:
        """Perform functional annotation using Variant Effect Predictor (VEP) REST API.

        Args:
            variant_coordinates (VariantCoordinates): A VariantCoordinates object
        Returns:
            typing.Sequence[TranscriptAnnotation]: A sequence of transcript
            annotations for the variant coordinates
        """
        response = self._query_vep(variant_coordinates)
        annotations = []
        if 'transcript_consequences' not in response:
            raise ValueError(
                f'The VEP response lacked the required `transcript_consequences` field')

        for trans in response['transcript_consequences']:
            annotation = self._process_item(trans)
            if annotation is not None:
                annotations.append(annotation)

        return annotations

    def _process_item(self, item) -> typing.Optional[TranscriptAnnotation]:
        """
        Parse one transcript annotation from the JSON response.
        """
        trans_id = item.get('transcript_id')
        if not self._include_computational_txs and not trans_id.startswith('NM_'):
            # Skipping a computational transcript
            return None
        is_preferred = True if 'canonical' in item and item['canonical'] else False
        hgvsc_id = item.get('hgvsc')
        consequences = item.get('consequence_terms')
        gene_name = item.get('gene_symbol')
        protein_id = item.get('protein_id')
        protein = self._protein_annotator.annotate(protein_id)
        protein_effect_start = item.get('protein_start')
        protein_effect_end = item.get('protein_end')
        if protein_effect_start is None and protein_effect_end is not None:
            protein_effect_start = 1
        if protein_effect_end is not None:
            protein_effect_end = int(protein_effect_end)
        if protein_effect_start is not None:
            protein_effect_start = int(protein_effect_start)
        exons_effected = item.get('exon')
        if exons_effected is not None:
            exons_effected = exons_effected.split('/')[0].split('-')
            if len(exons_effected) == 2:
                exons_effected = range(int(exons_effected[0]),
                                       int(exons_effected[1]) + 1)
            exons_effected = [int(x) for x in exons_effected]
        return TranscriptAnnotation(gene_name,
                                    trans_id,
                                    hgvsc_id,
                                    is_preferred,
                                    consequences,
                                    exons_effected,
                                    protein,
                                    protein_effect_start,
                                    protein_effect_end)

    def _query_vep(self, variant_coordinates) -> dict:
        api_url = self._url % (verify_start_end_coordinates(variant_coordinates))
        r = requests.get(api_url, headers={'Content-Type': 'application/json'})
        if not r.ok:
            self._logging.error(
                f"Expected a result but got an Error for variant: "
                f"{variant_coordinates.as_string()}")
            r.raise_for_status()
        results = r.json()
        if not isinstance(results, list):
            self._logging.error(results.get('error'))
            raise ConnectionError(
                f"Expected a result but got an Error. See log for details.")
        if len(results) > 1:
            self._logging.error([result.id for result in results])
            raise ValueError(
                f"Expected only one variant per request but received {len(results)} "
                f"different variants.")
        return results[0]


class VepHgvsVariantCoordinateFinder(VariantCoordinateFinder[typing.Tuple[str, str]]):
    """
    `VepHgvsVariantCoordinateFinder` uses Ensembl's Variant Effect Predictor (VEP)
    REST API to build
    `VariantCoordinates` from an HGVS string.

    The finder takes a tuple with two strings:

    * HGVS `str` (e.g. `NM_005912.3:c.253A>G`)
    * genotype `str` (e.g. `heterozygous`)

    and extracts the variant coordinates from the response.

    """

    def __init__(self, genome_build: GenomeBuild):
        self._build = hpotk.util.validate_instance(genome_build, GenomeBuild,
                                                   'genome_build')
        self._url = 'https://rest.ensembl.org/vep/human/hgvs/%s?refseq=1'

    def find_coordinates(self, item: T) -> typing.Tuple[VariantCoordinates, Genotype]:
        # TODO - implement
        #  `item` is a tuple described in class docstring.
        #  - construct a URL to query the VEP endpoint, use the placeholder `self._url`.
        #    We may want to include sanity check before query execution; we may not
        #    want to query VEP
        #    with rubbish values. Perhaps a regexp to check if the input looks
        #    roughly like a HGVS string?
        #     - has RefSeq transcript identifier (NM_...) with a version
        #     - represents annotation with respect to a coding sequence: `c.123C>T`,
        #     `c.100_110del`, etc.
        #       (not a protein (p.Met123Gly))
        #  - query the endpoint to get response
        #  - extract `VariantCoordinates` from the response
        #  - process the genotype string into `Genotype` enum member
        #  - return the results in a tuple
        hgvs, genotype = item
        request_url = self._url % hgvs
        headers = {'Content-type': 'application/json'}
        pattern = re.compile(r'^NM_\d+\.\d+:c\.(?:\d+|dup|del|ins)[ACGT]*>[ACGT]*$')
        if pattern.match(hgvs):
            response = requests.get(request_url, headers=headers)
            response = response.json()
            print(response)

            variant_coordinates = self._extract_variant_coordinates(response)

            genotype = Genotype[genotype.upper()]
            return variant_coordinates, genotype
        else:
            raise ValueError(f'Invalid HGVS string: {hgvs}')

    def _extract_variant_coordinates(self, response: typing.List[typing.Dict]) \
            -> typing.Optional[VariantCoordinates]:
        # TODO - implement
        #  `response` is a dict with JSON response.
        #  An example response is at `test_data/vep_response/hgvs_missense.json`,
        #  see `test_data/vep_response/README.md`
        #  for more info.
        #  Here we need to extract the relevant pieces and form `VariantCoordinates`.
        #  Note, you can get a `Contig` from `self._build`.
        #
        #  I prepared unit tests with examples of variants which we we are most
        #  likely to encounter.
        #  We have the following variant categories: missense, deletion, insertion,
        #  and duplication.
        #  If the tests in `_test_vep/TestVepHgvsVariantCoordinateFinder` pass then
        #  we're practically done! 😎
        response = response[0]

        chrom = int(response['seq_region_name'])
        transcript_consequences = response['transcript_consequences'][0]
        strand = transcript_consequences['strand']

        variant_coordinates = VariantCoordinates(
            region=GenomicRegion(
                contig=self._build.contigs[chrom-1],
                start=response['start'],
                end=response['end'],
                strand=(
                    Strand['NEGATIVE'] if strand == -1
                    else Strand['POSITIVE'] if strand == 1
                    else None
                ),
            ),
            ref=transcript_consequences['given_ref'],
            alt=transcript_consequences['variant_allele'],
            change_length=response['end'] - response['start'] + 1,
        )
        return variant_coordinates


if __name__ == '__main__':
    vep = VepHgvsVariantCoordinateFinder(GenomeBuild('GRCh38.p13', GRCh38.contigs))
    print(vep.find_coordinates(('NM_005912.3:c.253A>G', 'heterozygous')))
