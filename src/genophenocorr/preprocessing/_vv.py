# A module with classes that interact with Variant validator's REST API to fetch required data.
import logging
import re
import typing

import hpotk
import requests
import json


from genophenocorr.model import VariantCoordinates, TranscriptInfoAware, TranscriptCoordinates
from genophenocorr.model.genome import GenomeBuild, GenomicRegion, Strand, Contig, transpose_coordinate
from ._api import VariantCoordinateFinder, TranscriptCoordinateService


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

    def find_coordinates(self, item: str) -> VariantCoordinates:
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
                response = requests.get(request_url, headers=self._headers, timeout=self._timeout)
                response.raise_for_status()
                response = response.json()
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


class VVTranscriptCoordinateService(TranscriptCoordinateService):
    """
    `VVTranscriptCoordinateService` fetches the transcript coordinates from the Variant Validator REST API.

    :param genome_build: the genome build for constructing the transcript coordinates.
    :param timeout: a positive `float` with the REST API timeout in seconds.
    """

    def __init__(self, genome_build: GenomeBuild, timeout: float = 30.):
        self._logger = logging.getLogger(__name__)
        self._genome_build = hpotk.util.validate_instance(genome_build, GenomeBuild, 'genome_build')

        self._timeout = hpotk.util.validate_instance(timeout, float, 'timeout')
        if self._timeout <= 0:
            raise ValueError(f'`timeout` must be a positive `float` but got {timeout}')

        self._url = "https://rest.variantvalidator.org/VariantValidator/tools/gene2transcripts/%s"
        self._headers = {'Accept': 'application/json'}

    def fetch(self, tx: typing.Union[str, TranscriptInfoAware]) -> TranscriptCoordinates:
        tx_id = self._parse_tx(tx)
        api_url = self._url % tx_id

        response = requests.get(api_url, headers=self._headers, timeout=self._timeout)

        if not response.ok:
            response.raise_for_status()

        return self.parse_response(tx_id, response.json())

    @staticmethod
    def _parse_tx(tx: typing.Union[str, TranscriptInfoAware]) -> str:
        if isinstance(tx, str):
            return tx
        elif isinstance(tx, TranscriptInfoAware):
            return tx.transcript_id
        else:
            raise ValueError(f'Expected a `str` or `TranscriptInfoAware` but got {type(tx)}: {tx}')

    def parse_response(self, tx_id: str, response) -> TranscriptCoordinates:
        if not isinstance(response, list):
            transcript_response = response
        else:
            if len(response) != 1:
                self._logger.warning('Response has %s!=1 items. Choosing the first', len(response))
            transcript_response = response[0]
        if 'requested_symbol' not in transcript_response or transcript_response['requested_symbol'] != tx_id:
            print("########### Not able to parse VariantValidator Response ##############")
            #print(f"     Not able to find {tx_id} in the `requested_symbol` field in the response from Variant Validator API")
            print("This is the response:")
            json_formatted_str = json.dumps(response, indent=2)
            print(json_formatted_str)
            #raise ValueError(
            #    f'Could not find {tx_id} in the `requested_symbol` field in the response from Variant Validator API'
            #)
        if 'transcripts' not in transcript_response:
            raise ValueError(f'A required `transcripts` field is missing in the response from Variant Validator API')
        tx = self._find_tx_data(tx_id, transcript_response['transcripts'])

        if 'genomic_spans' not in tx:
            raise ValueError(f'A required `genomic_spans` field is missing in the response from Variant Validator API')
        contig, genomic_span = self._find_genomic_span(tx['genomic_spans'])

        strand = self._parse_strand(genomic_span)
        start_pos = genomic_span['start_position']
        end_pos = genomic_span['end_position']
        region = self._parse_tx_region(contig, start_pos, end_pos, strand)

        exons = self._parse_exons(contig, strand, genomic_span)

        cds_start, cds_end = self._parse_cds(tx['coding_start'], tx['coding_end'], exons)

        return TranscriptCoordinates(tx_id, region, exons, cds_start, cds_end)

    def _find_tx_data(self, tx_id: str, txs: typing.Dict) -> typing.Dict:
        for i, tx_data in enumerate(txs):
            if 'reference' not in tx_data:
                self._logger.warning('Skipping `transcripts` element %d due to missing `reference` field', i)
                continue
            if tx_data['reference'] == tx_id:
                return tx_data

        # If we get here, we did not find the transcript data, so we raise.
        raise ValueError(f'Did not find transcript data for the requested transcript {tx_id}')

    def _find_genomic_span(self, genomic_spans: typing.Dict) -> typing.Tuple[Contig, typing.Dict]:
        for contig_name, value in genomic_spans.items():
            contig = self._genome_build.contig_by_name(contig_name)
            if contig is not None:
                return contig, value

        contig_names = sorted(genomic_spans.keys())
        raise ValueError(f'Contigs {contig_names} were not found in genome build `{self._genome_build.identifier}`')

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

        raise ValueError(f'Could not parse CDS start and end from given coordinates')
