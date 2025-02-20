import matplotlib.pyplot as plt
import pytest

from gpsea.model import Cohort, ProteinMetadata
from gpsea.view import (
    configure_default_protein_visualizer,
    BaseProteinVisualizer,
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

    def test_drawing_with_empty_cohort(
        self,
        visualizer: BaseProteinVisualizer,
        suox_protein_metadata: ProteinMetadata,
    ):
        empty = Cohort(members=(), excluded_member_count=0)

        _, ax = plt.subplots(figsize=(20, 20))
        with pytest.raises(ValueError) as e:
            visualizer.draw_protein(
                cohort=empty,
                protein_metadata=suox_protein_metadata,
                ax=ax,
            )

        assert e.value.args == (
            'No variants annotated with respect to "NP_001027558.1" were found',
        )
