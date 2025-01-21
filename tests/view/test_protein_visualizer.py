import io
import matplotlib.pyplot as plt
import pytest

from gpsea.model import Cohort, ProteinMetadata
from gpsea.view import (
    GpseaReport,
    configure_default_protein_visualizer,
    BaseProteinVisualizer,
    ProteinVariantViewer,
)


class TestProteinVisualizer:
    @pytest.fixture(scope="class")
    def visualizer(self) -> BaseProteinVisualizer:
        return configure_default_protein_visualizer()

    @pytest.mark.skip("Run manually on demand")
    def test_protein_visualizer(
        self,
        visualizer: BaseProteinVisualizer,
        suox_cohort: Cohort,
        suox_protein_metadata: ProteinMetadata,
    ):
        fig, ax = plt.subplots(figsize=(20, 20))
        visualizer.draw_protein(
            cohort=suox_cohort,
            protein_metadata=suox_protein_metadata,
            ax=ax,
        )

        fig.savefig("protein.png")

    def test_protein_viewable(
        self,
        suox_cohort: Cohort,
        suox_protein_metadata: ProteinMetadata,
    ):
        protein_viewable = ProteinVariantViewer(
            protein_metadata=suox_protein_metadata,
        )
        report = protein_viewable.process(suox_cohort)
        assert isinstance(report, GpseaReport)

        buf = io.StringIO()
        report.write(buf)
        val = buf.getvalue()

        assert "gpsea-body" in val
