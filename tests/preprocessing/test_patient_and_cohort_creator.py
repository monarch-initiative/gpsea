import os

import hpotk
import pytest

from genophenocorr.model.genome import GenomeBuild
from genophenocorr.preprocessing import PhenotypeCreator, FunctionalAnnotator, ImpreciseSvFunctionalAnnotator, VariantCoordinateFinder
from genophenocorr.preprocessing import VepFunctionalAnnotator, VarCachingFunctionalAnnotator, NoOpImpreciseSvFunctionalAnnotator, VVHgvsVariantCoordinateFinder
from genophenocorr.preprocessing import PhenopacketPatientCreator
from genophenocorr.preprocessing import CohortCreator, load_phenopacket_folder


class TestPhenopacketCohortCreator:

    @pytest.fixture
    def phenotype_creator(
        self,
        hpo: hpotk.MinimalOntology,
        validation_runner: hpotk.validate.ValidationRunner,
    ) -> PhenotypeCreator:
        return PhenotypeCreator(
            hpo=hpo,
            validator=validation_runner,
        )

    @pytest.fixture
    def functional_annotator(
        self,
        fpath_project_dir: str,
    ) -> FunctionalAnnotator:
        fpath_cache_dir = os.path.join(fpath_project_dir, '.genophenocorr_cache')
        fpath_variant_cache_dir = os.path.join(fpath_cache_dir, 'variant_cache')
        os.makedirs(fpath_variant_cache_dir, exist_ok=True)
        
        return VarCachingFunctionalAnnotator.with_cache_folder(
            fpath_cache_dir=fpath_variant_cache_dir,
            fallback=VepFunctionalAnnotator(
                timeout=20,
            ),
        )

    @pytest.fixture
    def imprecise_so_functional_annotator(self) -> ImpreciseSvFunctionalAnnotator:
        return NoOpImpreciseSvFunctionalAnnotator()

    @pytest.fixture
    def variant_coordinate_finder(
        self,
        genome_build: GenomeBuild,
    ) -> VariantCoordinateFinder:
        return VVHgvsVariantCoordinateFinder(
            genome_build=genome_build,
        )

    @pytest.fixture
    def phenopacket_cohort_creator(
        self,
        genome_build: GenomeBuild,
        phenotype_creator: PhenotypeCreator,
        functional_annotator: FunctionalAnnotator,
        imprecise_so_functional_annotator: ImpreciseSvFunctionalAnnotator,
        variant_coordinate_finder: VariantCoordinateFinder,
    ) -> CohortCreator:
        patient_creator = PhenopacketPatientCreator(
            build=genome_build,
            phenotype_creator=phenotype_creator,
            var_func_ann=functional_annotator,
            imprecise_sv_functional_annotator=imprecise_so_functional_annotator,
            hgvs_coordinate_finder=variant_coordinate_finder,
        )
        return CohortCreator(
            patient_creator=patient_creator,
        )

    @pytest.mark.skip('Skipping online test')
    def test_load_phenopacket(
        self, 
        fpath_project_dir: str,
        phenopacket_cohort_creator: CohortCreator,
    ):
        fpath_test_cohort = os.path.join(fpath_project_dir, 'docs', 'data', 'simple_cohort')
        cohort = load_phenopacket_folder(
            pp_directory=fpath_test_cohort, 
            cohort_creator=phenopacket_cohort_creator,
        )
        print(cohort)
