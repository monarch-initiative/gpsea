import matplotlib.pyplot as plt
import pytest

from genophenocorr.model import TranscriptCoordinates, ProteinMetadata, Cohort
from genophenocorr.view import ProteinVisualizer, ProteinVisualizable


class TestProteinVisualizer:

    @pytest.fixture
    def visualizer(self) -> ProteinVisualizer:
        return ProteinVisualizer()

    @pytest.fixture
    def visualizable(
            self,
            suox_mane_tx_coordinates: TranscriptCoordinates,
            suox_protein_metadata: ProteinMetadata,
            suox_cohort: Cohort,
    ) -> ProteinVisualizable:
        return ProteinVisualizable(
            tx_coordinates=suox_mane_tx_coordinates,
            protein_meta=suox_protein_metadata,
            cohort=suox_cohort,
        )

    @pytest.mark.skip('Run manually on demand')
    def test_protein_visualizer(
            self,
            visualizer: ProteinVisualizer,
            visualizable: ProteinVisualizable,
    ):
        fig, ax = plt.subplots(figsize=(20, 20))
        visualizer.draw_fig(
            pvis=visualizable,
            ax=ax,
        )

        fig.savefig('protein.png')
