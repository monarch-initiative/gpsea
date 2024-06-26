import typing

import pytest

from genophenocorr.analysis.predicate.genotype import VariantPredicates, ProteinPredicates
from genophenocorr.model import *
from genophenocorr.model.genome import *
from genophenocorr.preprocessing import ProteinMetadataService


class TestVariantPredicates:

    @pytest.mark.parametrize(
        'effect, expected',
        [
            (VariantEffect.MISSENSE_VARIANT, True),
            (VariantEffect.SPLICE_DONOR_VARIANT, True),
            (VariantEffect.STOP_GAINED, False),
            (VariantEffect.STOP_LOST, False),
        ]
    )
    def test_variant_effect_predicate(
            self,
            variant: Variant,
            effect: VariantEffect,
            expected: bool,
    ):
        predicate = VariantPredicates.variant_effect(effect, tx_id='tx:xyz')

        assert predicate.test(variant) == expected

    @pytest.mark.parametrize(
        'variant_key, expected',
        [
            ('22_101_101_C_G', True),
            ('22_101_101_C_T', False),
        ]
    )
    def test_variant_key_predicate(
            self,
            variant: Variant,
            variant_key: str,
            expected: bool,
    ):
        predicate = VariantPredicates.variant_key(variant_key)

        assert predicate.test(variant) == expected

    @pytest.mark.parametrize(
        'exon, expected',
        [
            (0, False),
            (4, True),
            (5, False),
        ]
    )
    def test_exon_predicate(
            self,
            variant: Variant,
            exon: int,
            expected: bool,
    ):
        predicate = VariantPredicates.exon(exon, tx_id='tx:xyz')

        assert predicate.test(variant) == expected

    @pytest.mark.parametrize(
        'tx_id, expected',
        [
            ('tx:xyz', True),
            ('other', False),
        ]
    )
    def test_transcript_predicate(
            self,
            variant: Variant,
            tx_id: str,
            expected: bool,
    ):
        predicate = VariantPredicates.transcript(tx_id)

        assert predicate.test(variant) == expected

    @pytest.mark.parametrize(
        'symbol, expected',
        [
            ('a_gene', True),
            ('b_gene', False),
        ]
    )
    def test_gene_predicate(
            self,
            variant: Variant,
            symbol: str,
            expected: bool,
    ):
        predicate = VariantPredicates.gene(symbol)

        assert predicate.test(variant) == expected


class TestProteinPredicates:

    @pytest.fixture(scope='class')
    def protein_metadata_service(self) -> ProteinMetadataService:
        response = ProteinMetadata(
                protein_id='pt:xyz',
                label='xyz_label',
                protein_features=(
                    ProteinFeature.create(
                        FeatureInfo(name='MOCK_REPEAT', region=Region(55, 80)),
                        FeatureType.REPEAT,
                    ),
                    ProteinFeature.create(
                        FeatureInfo(name='MOCK_DOMAIN', region=Region(30, 50)),
                        FeatureType.DOMAIN,
                    ),
                )
            )
        return MockProteinMetadataService(response)

    @pytest.fixture
    def protein_predicates(
            self,
            protein_metadata_service: ProteinMetadataService,
    ) -> ProteinPredicates:
        return ProteinPredicates(protein_metadata_service)

    @pytest.mark.parametrize(
        'feature_type, expected',
        [
            (FeatureType.DOMAIN, True),
            (FeatureType.REPEAT, False),
        ]
    )
    def test_protein_feature_type(
            self,
            protein_predicates: ProteinPredicates,
            variant: Variant,
            feature_type: FeatureType,
            expected: bool,
    ):
        predicate = protein_predicates.protein_feature_type(feature_type, tx_id='tx:xyz')

        assert predicate.test(variant) == expected

    @pytest.mark.parametrize(
        'feature_id, expected',
        [
            ('MOCK_DOMAIN', True),
            ('REAL_DOMAIN', False),
            ('MOCK_REPEAT', False),
        ]
    )
    def test_protein_feature_id(
            self,
            protein_predicates: ProteinPredicates,
            variant: Variant,
            feature_id: str,
            expected: bool,
    ):
        predicate = protein_predicates.protein_feature(feature_id, tx_id='tx:xyz')

        assert predicate.test(variant) == expected


class MockProteinMetadataService(ProteinMetadataService):

    def __init__(
            self,
            response: ProteinMetadata
    ):
        self._response = response

    def annotate(self, protein_id: str) -> ProteinMetadata:
        return self._response


@pytest.fixture
def variant(genome_build: GenomeBuild) -> Variant:
    chr22 = genome_build.contig_by_name('chr22')
    assert chr22 is not None
    return Variant(
        var_coordinates=VariantCoordinates(
            region=GenomicRegion(
                contig=chr22,
                start=100,
                end=101,
                strand=Strand.POSITIVE,
            ),
            ref='C',
            alt='G',
            change_length=0,
        ),
        tx_annotations=(
            TranscriptAnnotation(
                gene_id='a_gene',
                tx_id='tx:xyz',
                hgvs_cdna=None,
                is_preferred=False,
                variant_effects=(
                    VariantEffect.MISSENSE_VARIANT,
                    VariantEffect.SPLICE_DONOR_VARIANT,
                ),
                affected_exons=(4,),
                protein_id='pt:xyz',
                protein_effect_coordinates=Region(40, 41),
            ),
        ),
        genotypes=Genotypes.single(SampleLabels('jim'), Genotype.HETEROZYGOUS)
    )
