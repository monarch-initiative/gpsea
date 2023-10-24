# A module with classes that interact with Variant validator's REST API to fetch required data.
import re
import typing

import hpotk
import requests

from genophenocorr.model import VariantCoordinates, Genotype
from genophenocorr.model.genome import GenomeBuild, GenomicRegion, Strand
from ._api import VariantCoordinateFinder


class VVHgvsVariantCoordinateFinder(VariantCoordinateFinder[typing.Tuple[str, str]]):
    """
    `VVHgvsVariantCoordinateFinder` uses Variant Validator's REST API to build :class:`VariantCoordinates`
    from an HGVS string.

    The finder takes a tuple with two strings:

    * HGVS `str` (e.g. `NM_005912.3:c.253A>G`)
    * genotype `str` (e.g. `heterozygous`)

    and extracts the variant coordinates from the response.

    :param genome_build: the genome build to use to construct :class:`VariantCoordinates`.
    :param timeout: the REST API request timeout (default: 10).
    """

    def __init__(self, genome_build: GenomeBuild, timeout: int = 10):
        self._build = hpotk.util.validate_instance(genome_build, GenomeBuild, 'genome_build')
        self._timeout = timeout
        self._url = 'https://rest.variantvalidator.org/VariantValidator/variantvalidator/%s/%s/%s'
        self._headers = {'Content-type': 'application/json'}
        self._hgvs_pattern = re.compile(r'^(?P<tx>NM_\d+\.\d+):c.\d+(_\d+)?.*')

    def find_coordinates(self, item: typing.Tuple[str, str]) -> typing.Tuple[VariantCoordinates, Genotype]:
        """
        Extracts variant coordinates from an HGVS string using Variant Validator's REST API.
        :param item: Tuple of hgvs string and genotype string
        :return: variant coordinates and genotype
        """
        hgvs, genotype = item
        matcher = self._hgvs_pattern.match(hgvs)
        if matcher:
            transcript = matcher.group('tx')
            request_url = self._url % (self._build.genome_build_id.major_assembly, hgvs, transcript)

            try:
                response = requests.get(request_url, headers=self._headers, timeout=self._timeout)
                response.raise_for_status()
                response = response.json()
                variant_coordinates = self._extract_variant_coordinates(response)
            except (requests.exceptions.RequestException, VariantValidatorDecodeException) as e:
                raise ValueError(f'Error processing {hgvs}', e)
        else:
            # The HGVS did not pass validation by a regular expression.
            # Please submit an issue to our GitHub tracker if you believe
            # the HGVS expression is correct.
            raise ValueError(f'Invalid HGVS string: {hgvs}')

        genotype = Genotype[genotype.upper()]
        return variant_coordinates, genotype

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
            raise VariantValidatorDecodeException(f'Contig {variant_data["chr"]} was not found in build {self._build.identifier}')

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
