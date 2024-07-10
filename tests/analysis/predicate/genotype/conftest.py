import pytest

from genophenocorr.model import *
from genophenocorr.model.genome import *
from genophenocorr.preprocessing import ProteinMetadataService

from genophenocorr.analysis.predicate.genotype import ProteinPredicates


@pytest.fixture(scope="package")
def protein_metadata_service() -> ProteinMetadataService:
    response = ProteinMetadata(
        protein_id="pt:xyz",
        label="xyz_label",
        protein_features=(
            ProteinFeature.create(
                FeatureInfo(name="MOCK_REPEAT", region=Region(55, 80)),
                FeatureType.REPEAT,
            ),
            ProteinFeature.create(
                FeatureInfo(name="MOCK_DOMAIN", region=Region(30, 50)),
                FeatureType.DOMAIN,
            ),
        ),
    )
    return MockProteinMetadataService(response)


@pytest.fixture(scope="package")
def protein_predicates(
    protein_metadata_service: ProteinMetadataService,
) -> ProteinPredicates:
    return ProteinPredicates(protein_metadata_service)


class MockProteinMetadataService(ProteinMetadataService):

    def __init__(self, response: ProteinMetadata):
        self._response = response

    def annotate(self, protein_id: str) -> ProteinMetadata:
        return self._response


@pytest.fixture(scope="package")
def variant(genome_build: GenomeBuild) -> Variant:
    chr22 = genome_build.contig_by_name("chr22")
    assert chr22 is not None
    return Variant(
        var_coordinates=VariantCoordinates(
            region=GenomicRegion(
                contig=chr22,
                start=100,
                end=101,
                strand=Strand.POSITIVE,
            ),
            ref="C",
            alt="G",
            change_length=0,
        ),
        tx_annotations=(
            TranscriptAnnotation(
                gene_id="a_gene",
                tx_id="tx:xyz",
                hgvs_cdna=None,
                is_preferred=False,
                variant_effects=(
                    VariantEffect.MISSENSE_VARIANT,
                    VariantEffect.SPLICE_DONOR_VARIANT,
                ),
                affected_exons=(4,),
                protein_id="pt:xyz",
                protein_effect_coordinates=Region(40, 41),
            ),
        ),
        genotypes=Genotypes.single(SampleLabels("jim"), Genotype.HETEROZYGOUS),
    )
