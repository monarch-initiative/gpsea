import matplotlib.pyplot as plt
import pytest

from gpsea.model import TranscriptCoordinates, ProteinMetadata, Cohort
from gpsea.view import ProteinVisualizer, ProteinVisualizable, ProteinViewable


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

    @pytest.mark.skip('Run manually on demand')
    def test_protein_viewable(
        self, 
        suox_cohort: Cohort,
        visualizable: ProteinVisualizable,
    ):
        protein_viewable = ProteinViewable()
        view = protein_viewable.process(suox_cohort, visualizable)
        with open('protein_viewable.html', 'w') as fh:
            fh.write(view)        
