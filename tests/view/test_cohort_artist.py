import matplotlib.pyplot as plt

import pytest

from gpsea.model import Cohort
from gpsea.view import configure_default_cohort_artist, CohortArtist


class TestCohortArtist:
    @pytest.fixture(scope="class")
    def cohort_artist(self) -> CohortArtist:
        return configure_default_cohort_artist()

    @pytest.mark.skip("Run manually on demand")
    def test_draw_protein(
        self,
        suox_cohort: Cohort,
        suox_mane_protein_id: str,
        cohort_artist: CohortArtist,
    ):
        fig, ax = plt.subplots(figsize=(20, 20))
        cohort_artist.draw_protein(
            cohort=suox_cohort,
            protein_id=suox_mane_protein_id,
            ax=ax,
        )

        fig.savefig("protein.png")
