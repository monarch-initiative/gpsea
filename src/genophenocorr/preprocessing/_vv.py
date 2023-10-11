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
        self._build = hpotk.util.validate_instance(genome_build, GenomeBuild, 'genome_build')
        # An example URL:
        #  https://rest.variantvalidator.org/VariantValidator/variantvalidator/hg38/NM_005912.3:c.253A>G/mane_select?content-type=application/json
        #self._url = None # TODO

    def find_coordinates(self, item: typing.Tuple[str, str]) -> typing.Tuple[VariantCoordinates, Genotype]:
        # TODO - implement
        #  `item` is a tuple described in class docstring.
        #  - construct a URL to query the VEP endpoint, use the placeholder `self._url`.
        #    We may want to include sanity check before query execution; we may not want to query VEP
        #    with rubbish values. Perhaps a regexp to check if the input looks roughly like a HGVS string?
        #     - has RefSeq transcript identifier (NM_...) with a version
        #     - represents annotation with respect to a coding sequence: `c.123C>T`, `c.100_110del`, etc.
        #       (not a protein (p.Met123Gly))
        #  - query the endpoint to get response
        #  - extract `VariantCoordinates` from the response
        #  - process the genotype string into `Genotype` enum member
        #  - return the results in a tuple
        raise NotImplementedError

    def _extract_variant_coordinates(self, response) -> typing.Optional[VariantCoordinates]:
        # TODO - implement
        #  `response` is a dict with JSON response.
        #  An example response is at `test_data/vep_response/hgvs_missense.json`, see `test_data/vep_response/README.md`
        #  for more info.
        #  Here we need to extract the relevant pieces and form `VariantCoordinates`.
        #  Note, you can get a `Contig` from `self._build`.
        #
        #  I prepared unit tests with examples of variants which we we are most likely to encounter.
        #  We have the following variant categories: missense, deletion, insertion, and duplication.
        #  If the tests in `_test_vep/TestVepHgvsVariantCoordinateFinder` pass then we're practically done! 😎
        return None
