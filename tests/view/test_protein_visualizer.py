import io
import matplotlib.pyplot as plt
import pytest

from gpsea.model import TranscriptCoordinates, ProteinMetadata, Cohort
from gpsea.view import (
    GpseaReport,
    ProteinVisualizer,
    ProteinVisualizable,
    ProteinVariantViewer,
)


class TestProteinVisualizer:

    @pytest.fixture(scope="class")
    def visualizer(self) -> ProteinVisualizer:
        return ProteinVisualizer()

    @pytest.fixture(scope="class")
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

    @pytest.mark.skip("Run manually on demand")
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

        fig.savefig("protein.png")

    def test_protein_viewable(
        self,
        suox_cohort: Cohort,
        suox_protein_metadata: ProteinMetadata,
        suox_mane_tx_id: str,
    ):
        protein_viewable = ProteinVariantViewer(
            protein_metadata=suox_protein_metadata, tx_id=suox_mane_tx_id
        )
        report = protein_viewable.process(suox_cohort)
        assert isinstance(report, GpseaReport)

        buf = io.StringIO()
        report.write(buf)
        val = buf.getvalue()

        assert "gpsea-body" in val
