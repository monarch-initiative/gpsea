# A module with classes that interact with Variant validator's REST API to fetch required data.
import logging
import re
import typing

import hpotk
import requests
import ratelimit
import json


from gpsea.model import VariantCoordinates, TranscriptInfoAware, TranscriptCoordinates
from gpsea.model.genome import GenomeBuild, GenomicRegion, Strand, Contig, transpose_coordinate
from ._api import VariantCoordinateFinder, TranscriptCoordinateService, GeneCoordinateService

# To match strings such as `NM_1234.56`, `NM_1234`, or `XM_123456.7`
REFSEQ_TX_PT = re.compile(r'^[NX]M_\d+(\.\d+)?$')


@ratelimit.sleep_and_retry
@ratelimit.limits(calls=1, period=1.2)
def fetch_response(
    url, 
    headers, 
    timeout,
):
    # This is the only place we interact with the Variant validator REST API.
    # Per documentation at `https://rest.variantvalidator.org/`,
    # we must limit requests to up to 1 request per second (+ 100ms buffer).
    response = requests.get(url, headers=headers, timeout=timeout)
    response.raise_for_status()
    return response.json()


class VVHgvsVariantCoordinateFinder(VariantCoordinateFinder[str]):
    """
    `VVHgvsVariantCoordinateFinder` uses Variant Validator's REST API to build :class:`VariantCoordinates`
    from an HGVS string.

    The finder takes an HGVS `str` (e.g. `NM_005912.3:c.253A>G`) and extracts the variant coordinates from the response.

    :param genome_build: the genome build to use to construct :class:`VariantCoordinates`
    :param timeout: the REST API request timeout
    """

    def __init__(self, genome_build: GenomeBuild, timeout: int = 30):
        self._build = hpotk.util.validate_instance(genome_build, GenomeBuild, 'genome_build')
        self._timeout = timeout
        self._url = 'https://rest.variantvalidator.org/VariantValidator/variantvalidator/%s/%s/%s'
        self._headers = {'Content-type': 'application/json'}
        self._hgvs_pattern = re.compile(r'^(?P<tx>NM_\d+\.\d+):c.\d+(_\d+)?.*')

    def find_coordinates(self, item: str) -> typing.Optional[VariantCoordinates]:
        """
        Extracts variant coordinates from an HGVS string using Variant Validator's REST API.

        :param item: a hgvs string
        :return: variant coordinates
        """
        matcher = self._hgvs_pattern.match(item)
        if matcher:
            transcript = matcher.group('tx')
            request_url = self._url % (self._build.genome_build_id.major_assembly, item, transcript)

            try:
                response = fetch_response(request_url, self._headers, self._timeout)
                variant_coordinates = self._extract_variant_coordinates(response)
            except (requests.exceptions.RequestException, VariantValidatorDecodeException) as e:
                raise ValueError(f'Error processing {item}', e)
        else:
            # The HGVS did not pass validation by a regular expression.
            # Please submit an issue to our GitHub tracker if you believe
            # the HGVS expression is correct.
            raise ValueError(f'Invalid HGVS string: {item}')

        return variant_coordinates

    def _extract_variant_coordinates(self, response: typing.Dict) -> typing.Optional[VariantCoordinates]:
        """
        Extracts variant coordinates from a response from Variant Validator's REST API.
        :param response: Variant Validator's REST API response as a dictionary
        :return: variant coordinates
        """
        variant_identifier = list(response.keys())[0]
        response = response[variant_identifier]

        selected_assembly = response['selected_assembly']
        variant_data = response['primary_assembly_loci'][selected_assembly.lower()]['vcf']

        contig = self._build.contig_by_name(variant_data['chr'])
        if contig is None:
            raise VariantValidatorDecodeException(f'Contig {variant_data["chr"]} was not found in build '
                                                  f'{self._build.identifier}')

        pos = int(variant_data['pos'])
        ref = variant_data['ref']
        alt = variant_data['alt']

        start, end, change_length = self._extract_coordinates(pos, ref, alt)

        return VariantCoordinates(
            region=GenomicRegion(
                contig=contig,
                start=start,
                end=end,
                strand=Strand.POSITIVE,
            ),
            ref=ref,
            alt=alt,
            change_length=change_length,
        )

    @staticmethod
    def _extract_coordinates(pos: int, ref: str, alt: str) -> typing.Tuple[int, int, int]:
        """
        Extracts 0-based start, end and change_length from 1-based position, reference and alternative alleles.

        :param pos: 1-based position
        :param ref: reference allele
        :param alt: alternative allele
        :return: 0-based start, end and change_length
        """
        change_length = len(alt) - len(ref)
        start = pos - 1
        end = start + len(ref)

        return start, end, change_length


class VariantValidatorDecodeException(BaseException):
    """
    An exception thrown when parsing variant validator response fails.
    """
    pass


class VVMultiCoordinateService(TranscriptCoordinateService, GeneCoordinateService):
    """
    `VVMultiCoordinateService` uses the Variant Validator REST API to fetch transcript coordinates for
    both a *gene* ID and a specific *transcript* ID.

    :param genome_build: the genome build for constructing the transcript coordinates.
    :param timeout: a positive `float` with the REST API timeout in seconds.
    """

    def __init__(
        self,
        genome_build: GenomeBuild, 
        timeout: float = 30.,
    ):
        self._logger = logging.getLogger(__name__)
        self._genome_build = hpotk.util.validate_instance(genome_build, GenomeBuild, 'genome_build')

        self._timeout = timeout
        if self._timeout <= 0:
            raise ValueError(f'`timeout` must be a positive `float` but got {timeout}')

        self._url = "https://rest.variantvalidator.org/VariantValidator/tools/gene2transcripts/%s"
        self._headers = {'Accept': 'application/json'}

    def fetch(self, tx: typing.Union[str, TranscriptInfoAware]) -> TranscriptCoordinates:
        tx_id = self._parse_tx(tx)
        response_json = self.get_response(tx_id)
        return self.parse_response(tx_id, response_json)

    def fetch_for_gene(self, gene: str) -> typing.Sequence[TranscriptCoordinates]:
        response_json = self.get_response(gene)
        return self.parse_multiple(response_json)

    @staticmethod
    def _parse_tx(tx: typing.Union[str, TranscriptInfoAware]) -> str:
        if isinstance(tx, str):
            return tx
        elif isinstance(tx, TranscriptInfoAware):
            return tx.transcript_id
        else:
            raise ValueError(f'Expected a `str` or `TranscriptInfoAware` but got {type(tx)}: {tx}')

    def get_response(self, tx_id: str):
        api_url = self._url % tx_id
        return fetch_response(api_url, self._headers, self._timeout)

    def parse_response(self, tx_id: str, response) -> TranscriptCoordinates:
        if not isinstance(response, list):
            transcript_response = response
        else:
            if len(response) != 1:
                self._logger.warning('Response has %s!=1 items. Choosing the first', len(response))
            transcript_response = response[0]
        if 'requested_symbol' not in transcript_response or transcript_response['requested_symbol'] != tx_id:
            json_formatted_str = json.dumps(response, indent=2)
            error_string = f"Not able to find {tx_id} in the `requested_symbol` field in the response from Variant Validator API: \n{json_formatted_str}"
            raise ValueError(error_string)
        if 'transcripts' not in transcript_response:
            VVMultiCoordinateService._handle_missing_field(
                response=response, 
                field='transcripts',
            )
        tx_data = self._find_tx_data(tx_id, transcript_response['transcripts'])
        if 'genomic_spans' not in tx_data:
            raise ValueError(f'A required `genomic_spans` field is missing in the response from Variant Validator API')
        
        return VVMultiCoordinateService._parse_tx_coordinates(
            genome_build=self._genome_build,
            tx_id=tx_id,
            tx_data=tx_data,
        )

    def parse_multiple(
        self, 
        response: typing.Mapping[str, typing.Any],
    ) -> typing.Sequence[TranscriptCoordinates]:
        if 'transcripts' not in response:
            VVMultiCoordinateService._handle_missing_field(
                response=response, 
                field='transcripts',
            )
        transcripts: typing.Sequence[typing.Mapping[str, typing.Any]] = response['transcripts']
        coordinates = []
        for tx_data in transcripts:
            if 'reference' not in tx_data:
                VVMultiCoordinateService._handle_missing_field(
                    response=response, 
                    field='reference',
                )
            tx_id = tx_data['reference']
            matcher = REFSEQ_TX_PT.match(tx_id)
            if not matcher:
                self._logger.debug('Skipping processing transcript %s', tx_id)
                continue
            tx_coordinates = VVMultiCoordinateService._parse_tx_coordinates(
                genome_build=self._genome_build,
                tx_id=tx_id,
                tx_data=tx_data,
            )
            coordinates.append(tx_coordinates)

        return tuple(coordinates)

    def _find_tx_data(
        self, 
        tx_id: str, 
        txs: typing.Dict,
    ) -> typing.Dict:
        for i, tx_data in enumerate(txs):
            if 'reference' not in tx_data:
                self._logger.warning('Skipping `transcripts` element %d due to missing `reference` field', i)
                continue
            if tx_data['reference'] == tx_id:
                return tx_data

        # If we get here, we did not find the transcript data, so we raise.
        raise ValueError(f'Did not find transcript data for the requested transcript {tx_id}')

    @staticmethod
    def _handle_missing_field(response, field: str):
        json_formatted_str = json.dumps(response, indent=2)
        error_string = f"A required `{field}` field is missing in the response from Variant Validator API: \n{json_formatted_str}"
        raise ValueError(error_string)

    @staticmethod
    def _parse_tx_coordinates(
        genome_build: GenomeBuild,
        tx_id: str,
        tx_data: typing.Mapping[str, typing.Any],
    ) -> TranscriptCoordinates:
        contig, genomic_span = VVMultiCoordinateService._find_genomic_span(
            genome_build=genome_build, 
            genomic_spans=tx_data['genomic_spans'],
        )

        strand = VVMultiCoordinateService._parse_strand(genomic_span)
        start_pos = genomic_span['start_position']
        end_pos = genomic_span['end_position']
        region = VVMultiCoordinateService._parse_tx_region(contig, start_pos, end_pos, strand)

        exons = VVMultiCoordinateService._parse_exons(contig, strand, genomic_span)

        cds_start, cds_end = VVMultiCoordinateService._parse_cds(tx_data['coding_start'], tx_data['coding_end'], exons)

        is_preferred = VVMultiCoordinateService._parse_is_preferred(tx_data)

        return TranscriptCoordinates(
            identifier=tx_id, 
            region=region, 
            exons=exons, 
            cds_start=cds_start, 
            cds_end=cds_end,
            is_preferred=is_preferred,
        )

    @staticmethod
    def _find_genomic_span(
        genome_build: GenomeBuild,
        genomic_spans: typing.Dict,
    ) -> typing.Tuple[Contig, typing.Dict]:
        for contig_name, value in genomic_spans.items():
            contig = genome_build.contig_by_name(contig_name)
            if contig is not None:
                return contig, value

        contig_names = sorted(genomic_spans.keys())
        raise ValueError(f'Contigs {contig_names} were not found in genome build `{genome_build.identifier}`')

    @staticmethod
    def _parse_strand(genomic_span) -> Strand:
        orientation = genomic_span['orientation']
        if orientation == 1:
            return Strand.POSITIVE
        elif orientation == -1:
            return Strand.NEGATIVE
        else:
            raise ValueError(f'Invalid orientation value {orientation} was not {{-1, 1}}')

    @staticmethod
    def _parse_tx_region(contig: Contig, start_pos: int, end_pos: int, strand: Strand) -> GenomicRegion:
        """
        The `start_pos` and `end_pos` coordinates are on the positive strand, but the genomic region that we
        are returning must be on given `strand`. Therefore, we transpose accordingly.

        However, note that start_pos > end_pos if `strand == Strand.NEGATIVE`. This is unusual in our environment,
        where we ensure the following: `start <= end`.
        """
        if strand == Strand.POSITIVE:
            start = start_pos - 1  # Convert from 1-based to 0-based coordinate
            end = end_pos
        else:
            end = transpose_coordinate(contig, start_pos - 1)  # Convert from 1-based to 0-based coordinate
            start = transpose_coordinate(contig, end_pos)

        return GenomicRegion(contig, start, end, strand)

    @staticmethod
    def _parse_exons(contig: Contig, strand: Strand, genomic_span: typing.Dict) -> typing.Sequence[GenomicRegion]:
        # Ensure the exons are sorted in ascending order
        exons = []
        for exon in sorted(genomic_span['exon_structure'], key=lambda exon_data: exon_data['exon_number']):
            gen_start = exon['genomic_start'] - 1  # -1 to convert to 0-based coordinates.
            gen_end = exon['genomic_end']
            if strand.is_positive():
                start = gen_start
                end = gen_end
            else:
                start = transpose_coordinate(contig, gen_end)  # !
                end = transpose_coordinate(contig, gen_start)  # !
            gr = GenomicRegion(contig, start, end, strand)
            exons.append(gr)

        return exons

    @staticmethod
    def _parse_cds(coding_start: int, coding_end: int, exons: typing.Iterable[GenomicRegion]) -> typing.Tuple[int, int]:
        # Variant validator reports CDS start and end as n-th exonic bases. We need to figure out which exon overlaps
        # with `coding_start` and `coding_end` in the CDS space and then calculate the
        start = None
        end = None
        processed = 0
        for exon in exons:
            exon_len = len(exon)
            if start is None:
                if processed < coding_start <= processed + exon_len:
                    start = exon.start + coding_start - processed - 1  # `-1` to convert to 0-based coordinate!

            if end is None:
                if processed < coding_end <= processed + exon_len:
                    end = exon.start + coding_end - processed

            if start is not None and end is not None:
                return start, end

            processed += exon_len

        raise ValueError('Could not parse CDS start and end from given coordinates')

    @staticmethod
    def _parse_is_preferred(
        tx_data: typing.Mapping[str, typing.Any],
    ) -> typing.Optional[bool]:
        if 'annotations' in tx_data:
            annotations = tx_data['annotations']
            if 'mane_select' in annotations:
                assert isinstance(annotations['mane_select'], bool), '\'mane_select\' field must be `bool`'
                return annotations['mane_select']
        
        return None
