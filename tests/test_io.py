import json

from genophenocorr.model import Cohort
from genophenocorr.io import GenophenocorrJSONEncoder, GenophenocorrJSONDecoder


def test_round_trip(suox_cohort: Cohort):
    dumped = json.dumps(suox_cohort, cls=GenophenocorrJSONEncoder, indent=2)
    decoded = json.loads(dumped, cls=GenophenocorrJSONDecoder)

    assert suox_cohort == decoded
