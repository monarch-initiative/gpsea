# A module with classes that interact with Variant validator's REST API to fetch required data.
import re
import typing

import hpotk
import requests

from genophenocorr.model import VariantCoordinates, Genotype
from genophenocorr.model.genome import GenomeBuild, GenomicRegion, Strand, GRCh38
from ._api import VariantCoordinateFinder, T


class VVHgvsVariantCoordinateFinder(VariantCoordinateFinder[typing.Tuple[str, str]]):
    """
    `VVHgvsVariantCoordinateFinder` uses Variant Validator's REST API to build `VariantCoordinates` from an HGVS string.

    The finder takes a tuple with two strings:

    * HGVS `str` (e.g. `NM_005912.3:c.253A>G`)
    * genotype `str` (e.g. `heterozygous`)

    and extracts the variant coordinates from the response.

    URL: https://rest.variantvalidator.org/
    """
    TIME_OUT = 10  # TODO @ielis: please set this to a value you find reasonable
    select_options = ['all', 'raw', 'select', 'mane_select', 'mane', 'refseq_select']
    SELECT = select_options[2]

    def __init__(self, genome_build: GenomeBuild):
        self._build = hpotk.util.validate_instance(genome_build, GenomeBuild, 'genome_build')
        self._url = 'https://rest.variantvalidator.org/VariantValidator/variantvalidator/%s/%s/%s'
        self.hgvs_pattern = re.compile(r'^NM_\d+\.\d+:c.\d+(_\d+)?.*')

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
        transcript, _ = hgvs.split(':')
        print(f"URL: {self._url}")
        print(f"Values: {self._build.identifier_name}, {hgvs}")
        request_url = self._url % (self._build.identifier_name, hgvs, transcript)
        print(f"Request URL: {request_url}")
        headers = {'Content-type': 'application/json'}
        if self.hgvs_pattern.match(hgvs):
            try:
                response = requests.get(request_url, headers=headers, timeout=VVHgvsVariantCoordinateFinder.TIME_OUT)
                response.raise_for_status()
                response = response.json()
                print(response)
            except requests.exceptions.RequestException as e:
                print(f"Error: {e}")

            variant_coordinates = self._extract_variant_coordinates(response)

            genotype = Genotype[genotype.upper()]
            return variant_coordinates, genotype
        else:
            raise ValueError(f'Invalid HGVS string: {hgvs}')

    def _extract_variant_coordinates(self, response: typing.Dict) \
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
        variant_identifier = list(response.keys())[0]
        response = response[variant_identifier]

        selected_assembly = response['selected_assembly']
        variant_data = response['primary_assembly_loci'][selected_assembly.lower()]['vcf']
        strand = Strand.POSITIVE
        chrom = variant_data['chr']
        pos = int(variant_data['pos'])
        ref = variant_data['ref']
        alt = variant_data['alt']
        start, end, change_length = self._extract_0based_mutation_range_from_1based_pos(pos, ref, alt)

        variant_coordinates = VariantCoordinates(
            region=GenomicRegion(
                contig=self._build.contig_by_name(chrom),
                start=start,
                end=end,
                strand=strand,
            ),
            ref=ref,
            alt=alt,
            change_length=change_length,
        )
        return variant_coordinates

    def _extract_0based_mutation_range_from_1based_pos(
            self, pos: int, ref: str, alt: str
    ) -> typing.Tuple[int, int, int]:
        """
        :return: 0-based start, end and change_length
        """
        change_length = 1  # TODO: determine change length from ref and alt
        change_length = len(alt) - len(ref)
        end = pos
        start = pos - len(ref)

        return start, end, change_length
