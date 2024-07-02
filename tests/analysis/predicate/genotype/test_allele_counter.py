import pytest

from genophenocorr.analysis.predicate.genotype import AlleleCounter, VariantPredicate
from genophenocorr.model import *
from genophenocorr.model.genome import *


@pytest.fixture(scope='module')
def sample_labels() -> SampleLabels:
    return SampleLabels('A_III-1', meta_label='PMID_10580070_A_III')

@pytest.fixture(scope='module')
def chr1(genome_build: GenomeBuild) -> Contig:
    contig = genome_build.contig_by_name('chr1')
    assert contig is not None
    return contig

@pytest.fixture(scope='module')
def chrX(genome_build: GenomeBuild) -> Contig:
    contig = genome_build.contig_by_name('chrX')
    assert contig is not None
    return contig

@pytest.fixture(scope='module')
def het_lmna(
        sample_labels: SampleLabels,
        chr1: Contig,
) -> Variant:
    return Variant(
        var_coordinates=VariantCoordinates(
            region=GenomicRegion(
                contig=chr1, start=156_137_755, end=156_137_756,
                strand=Strand.POSITIVE,
            ),
            ref='C', alt='A', change_length=0,
        ),
        tx_annotations=(
            TranscriptAnnotation(
                gene_id='LMNA', tx_id='NM_170707.4', hgvs_cdna='NM_170707.4:c.1698+13C>A', is_preferred=True,
                variant_effects=(VariantEffect.INTRON_VARIANT,), affected_exons=None, protein_id='NP_733821.1',
                protein_effect_coordinates=None,
            ),
            TranscriptAnnotation(
                gene_id='LMNA', tx_id='NM_005572.4', hgvs_cdna='NM_005572.4:c.1711C>A', is_preferred=False,
                variant_effects=(VariantEffect.MISSENSE_VARIANT,), affected_exons=(10,), protein_id='NP_005563.1',
                protein_effect_coordinates=Region(570, 571),
            ),
        ),
        genotypes=Genotypes.single(sample_labels, Genotype.HETEROZYGOUS),
    )


@pytest.fixture(scope='module')
def hom_alt_lmna(
        sample_labels: SampleLabels,
        chr1: Contig,
) -> Variant:
    return Variant(
        var_coordinates=VariantCoordinates(
            region=GenomicRegion(
                contig=chr1, start=156_134_852, end=156_134_853,
                strand=Strand.POSITIVE,
            ),
            ref='G', alt='A', change_length=0,
        ),
        tx_annotations=(
            TranscriptAnnotation(
                gene_id='LMNA', tx_id='NM_170707.4', hgvs_cdna='NM_170707.4:c.688G>A', is_preferred=True,
                variant_effects=(VariantEffect.MISSENSE_VARIANT,), affected_exons=(4,), protein_id='NP_733821.1',
                protein_effect_coordinates=Region(229, 230),
            ),
            TranscriptAnnotation(
                gene_id='LMNA', tx_id='NM_005572.4', hgvs_cdna='NM_005572.4:c.688G>A', is_preferred=False,
                variant_effects=(VariantEffect.MISSENSE_VARIANT,), affected_exons=(4,), protein_id='NP_005563.1',
                protein_effect_coordinates=Region(229, 230),
            ),
        ),
        genotypes=Genotypes.single(sample_labels, Genotype.HOMOZYGOUS_ALTERNATE),
    )

@pytest.fixture(scope='module')
def hemi_dmd(
        sample_labels: SampleLabels,
        chrX: Contig,
) -> Variant:
    return Variant(
        var_coordinates=VariantCoordinates(
            region=GenomicRegion(
                contig=chrX, start=31_180_436, end=31_180_437,
                strand=Strand.POSITIVE,
            ),
            ref='C', alt='T', change_length=0,
        ),
        tx_annotations=(
            TranscriptAnnotation(
                gene_id='DMD', tx_id='NM_000109.4', hgvs_cdna='NM_000109.4:c.9995G>A', is_preferred=False,
                variant_effects=(VariantEffect.MISSENSE_VARIANT,), affected_exons=(69,),
                protein_id='NP_000100.3', protein_effect_coordinates=Region(start=3331, end=3332)
            ),
            TranscriptAnnotation(
                gene_id='DMD', tx_id='NM_004006.3', hgvs_cdna='NM_004006.3:c.10019G>A', is_preferred=True,
                variant_effects=(VariantEffect.MISSENSE_VARIANT,), affected_exons=(69,),
                protein_id='NP_003997.2', protein_effect_coordinates=Region(start=3339, end=3340)
            ),
        ),
        genotypes=Genotypes.single(sample_labels, Genotype.HEMIZYGOUS),
    )

class TestAlleleCounter:

    @pytest.fixture
    def patient(
            self,
            sample_labels: SampleLabels,
            het_lmna: Variant,
            hom_alt_lmna: Variant,
            hemi_dmd: Variant,
    ) -> Patient:
        return Patient(
            labels=sample_labels,
            variants=(
                het_lmna,
                hom_alt_lmna,
                hemi_dmd,
            ),
            diseases=(),
            phenotypes=(),
        )

    @pytest.mark.parametrize(
        'variant_key, expected',
        [
            ('invariant', 0),
            ('1_156137756_156137756_C_A', 1),
            ('1_156134853_156134853_G_A', 2),
            ('X_31180437_31180437_C_T', 1),
        ]
    )
    def test_count_keys(
            self,
            patient: Patient,
            variant_key: str,
            expected: int,
    ):
        predicate = MockVariantKeyPredicate(variant_key)
        counter = AlleleCounter(predicate)

        assert counter.count(patient) == expected

    @pytest.mark.parametrize(
        'variant_effect, tx_id, expected',
        [
            (VariantEffect.START_LOST, 'NM_170707.4', 0),
            (VariantEffect.INTRON_VARIANT, 'NM_170707.4', 1),
            (VariantEffect.MISSENSE_VARIANT, 'NM_170707.4', 2),
            (VariantEffect.MISSENSE_VARIANT, 'NM_005572.4', 3),
        ]
    )
    def test_count_effects(
        self, 
        patient: Patient, 
        variant_effect: VariantEffect,
        tx_id: str,
        expected: int
    ):
        predicate = MockVariantEffectPredicate(variant_effect, tx_id)
        counter = AlleleCounter(predicate)
        
        assert counter.count(patient) == expected


class MockVariantKeyPredicate(VariantPredicate):

    def __init__(
            self,
            variant_key: str
    ):
        self._variant_key = variant_key

    def test(self, variant: Variant) -> bool:
        return variant.variant_coordinates.variant_key == self._variant_key

class MockVariantEffectPredicate(VariantPredicate):
    
    def __init__(self,
                 variant_effect: VariantEffect, tx_id: str) -> None:
        self._variant_effect = variant_effect
        self._tx_id = tx_id
        
    def test(self, variant: Variant) -> bool:
        tx_anno = variant.get_tx_anno_by_id(self._tx_id)
        if tx_anno is None:
            return False
        for effect in tx_anno.variant_effects:
            if effect == self._variant_effect:
                return True
        return False