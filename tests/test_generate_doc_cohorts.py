import json
import os

import pytest

import hpotk

from ppktstore.registry import PhenopacketStoreRegistry, configure_phenopacket_registry

from gpsea.io import GpseaJSONEncoder
from gpsea.model import Cohort
from gpsea.preprocessing import load_phenopackets
from gpsea.preprocessing import CohortCreator, configure_caching_cohort_creator


@pytest.mark.skip("Run only to update cohorts")
class TestGenerateCohortsForDocumentation:

    PHENOPACKET_STORE_VERSION = "0.1.20"

    @pytest.fixture
    def cohort_creator(
        self,
        hpo: hpotk.MinimalOntology,
    ) -> CohortCreator:
        return configure_caching_cohort_creator(
            hpo=hpo,
        )

    @pytest.fixture
    def phenopacket_registry(self) -> PhenopacketStoreRegistry:
        return configure_phenopacket_registry()

    @pytest.fixture
    def fpath_doc_data_dir(
        self,
        fpath_project_dir: str,
    ) -> str:
        return os.path.join(fpath_project_dir, "docs", "cohort-data")

    @staticmethod
    def load_cohort(
        cohort_name: str,
        cohort_creator: CohortCreator,
        phenopacket_registry: PhenopacketStoreRegistry,
    ) -> Cohort:
        with phenopacket_registry.open_phenopacket_store(
            TestGenerateCohortsForDocumentation.PHENOPACKET_STORE_VERSION,
        ) as ps:
            phenopackets = tuple(ps.iter_cohort_phenopackets(cohort_name))

        cohort, qc = load_phenopackets(
            phenopackets=iter(phenopackets),
            cohort_creator=cohort_creator,
        )
        assert qc.is_ok(), "There must be no Q/C issues in the cohort"

        return cohort

    @staticmethod
    def dump_cohort(
        cohort_name: str,
        cohort: Cohort,
        fpath_doc_dir: str,
    ):
        fpath_cohort_file = os.path.join(
            fpath_doc_dir,
            f"{cohort_name}.{TestGenerateCohortsForDocumentation.PHENOPACKET_STORE_VERSION}.json",
        )
        with open(fpath_cohort_file, "w") as fh:
            json.dump(
                cohort, fh,
                cls=GpseaJSONEncoder,
                indent=2,
            )

    @pytest.mark.parametrize(
        "cohort_name",
        [
            "CYP21A2",
            "RERE",
            "TBX5",
        ]
    )
    def test_generate_cohorts(
        self,
        cohort_name: str,
        phenopacket_registry: PhenopacketStoreRegistry,
        cohort_creator: CohortCreator,
        fpath_doc_data_dir: str,
    ):
        """
        Generate cohorts used in the tutorial and documentation.
        """
        cohort = TestGenerateCohortsForDocumentation.load_cohort(
            cohort_name=cohort_name,
            cohort_creator=cohort_creator,
            phenopacket_registry=phenopacket_registry,
        )

        TestGenerateCohortsForDocumentation.dump_cohort(
            cohort_name=cohort_name,
            cohort=cohort,
            fpath_doc_dir=fpath_doc_data_dir,
        )
