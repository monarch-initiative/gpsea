# A module with classes that interact with Variant validator's REST API to fetch required data.
import typing

import hpotk

from genophenocorr.model import VariantCoordinates, Genotype
from genophenocorr.model.genome import GenomeBuild
from ._api import VariantCoordinateFinder


class VVHgvsVariantCoordinateFinder(VariantCoordinateFinder[typing.Tuple[str, str]]):
    """
    `VVHgvsVariantCoordinateFinder` uses Variant Validator's REST API to build `VariantCoordinates` from an HGVS string.

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
                    Strand.NEGATIVE if strand == -1
                    else Strand.POSITIVE if strand == 1
                    else None
                ),
            ),
            ref=transcript_consequences['given_ref'],
            alt=transcript_consequences['variant_allele'],
            change_length=response['end'] - response['start'] + 1,
        )
        return variant_coordinates
