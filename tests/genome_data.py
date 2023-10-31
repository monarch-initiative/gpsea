import typing

import pytest

from genophenocorr.model import TranscriptCoordinates, ProteinMetadata, Variant, VariantCoordinates, TranscriptAnnotation, Genotypes, VariantEffect
from genophenocorr.model.genome import Contig, GenomicRegion, Strand


@pytest.fixture
def toy_contig() -> Contig:
    return Contig('1', 'GB_ACC', 'REFSEQ_ACC', 'USCD_ACC', 2_000)


@pytest.fixture
def toy_tx_positive(toy_contig: Contig) -> TranscriptCoordinates:
    return TranscriptCoordinates('fake-tx-pos',
                                 GenomicRegion(toy_contig, 200, 800, Strand.POSITIVE),
                                 (
                                     GenomicRegion(toy_contig, 200, 250, Strand.POSITIVE),
                                     GenomicRegion(toy_contig, 300, 340, Strand.POSITIVE),
                                     GenomicRegion(toy_contig, 510, 540, Strand.POSITIVE),
                                     GenomicRegion(toy_contig, 650, 800, Strand.POSITIVE)
                                 ),
                                 310, 740)


@pytest.fixture
def toy_protein_positive() -> ProteinMetadata:
    return ProteinMetadata('fake-prot-pos', 'FAKE PROTEIN ON POSITIVE STRAND', ())


@pytest.fixture
def toy_tx_negative(toy_contig: Contig) -> TranscriptCoordinates:
    return TranscriptCoordinates('fake-tx-neg',
                                 GenomicRegion(toy_contig, 100, 600, Strand.NEGATIVE),
                                 (
                                     GenomicRegion(toy_contig, 100, 180, Strand.NEGATIVE),
                                     GenomicRegion(toy_contig, 210, 270, Strand.NEGATIVE),
                                     GenomicRegion(toy_contig, 400, 450, Strand.NEGATIVE),
                                     GenomicRegion(toy_contig, 550, 600, Strand.NEGATIVE)
                                 ),
                                 130, 430)


@pytest.fixture
def toy_protein_negative() -> ProteinMetadata:
    return ProteinMetadata('fake-prot-neg', 'FAKE PROTEIN ON NEGATIVE STRAND', ())


@pytest.fixture
def toy_variants(toy_contig: Contig) -> typing.Sequence[Variant]:
    return (
        Variant(
            VariantCoordinates(GenomicRegion(toy_contig, 325, 326, Strand.POSITIVE),'C', 'T', 0),
            (
                    TranscriptAnnotation('some-gene', 'fake-tx-pos', 'fake-tx-pos-hgvsc:v1',
                                         True, (VariantEffect.MISSENSE_VARIANT,),
                                         (2,),
                                         (), None),
                ), Genotypes.empty()
                ),
        Variant(VariantCoordinates(GenomicRegion(toy_contig, 530, 531, Strand.POSITIVE),
                                   'C', 'CCT', 2),
                (
                    TranscriptAnnotation('some-gene', 'fake-tx-pos', 'fake-tx-pos-hgvsc:v2',
                                         True, (VariantEffect.FRAMESHIFT_VARIANT,),
                                         (3,),
                                         (), None),
                ), Genotypes.empty()
                ),
        Variant(VariantCoordinates(GenomicRegion(toy_contig, 160, 161, Strand.NEGATIVE).with_strand(Strand.POSITIVE),
                                   'G', 'A', 0),
                (
                    TranscriptAnnotation('other-gene', 'fake-tx-neg', 'fake-tx-neg-hgvsc:v3',
                                         True, (VariantEffect.SYNONYMOUS_VARIANT,),
                                         (1,),
                                         (), None),
                ), Genotypes.empty()
                ),
        Variant(VariantCoordinates(GenomicRegion(toy_contig, 570, 574, Strand.NEGATIVE).with_strand(Strand.POSITIVE),
                                   'CTGA', 'C', -3),
                (
                    TranscriptAnnotation('other-gene', 'fake-tx-neg', 'fake-tx-neg-hgvsc:v4',
                                         True, (VariantEffect.THREE_PRIME_UTR_VARIANT,),
                                         (4,),
                                         (), None),
                ), Genotypes.empty()
                )

    )
