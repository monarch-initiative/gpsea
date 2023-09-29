import pytest

from ._genome import Contig, GenomicRegion, Strand


@pytest.fixture
def contig() -> Contig:
    return Contig('1', 'GB_BLA', 'NC_BLA', 'UCSC_BLA', 100)


def test_start_on_strand(contig):
    region = GenomicRegion(contig, 30, 80, Strand.POSITIVE)

    assert region.start_on_strand(Strand.POSITIVE) == 30
    assert region.start_on_strand(Strand.NEGATIVE) == 20


def test_end_on_strand(contig):
    region = GenomicRegion(contig, 30, 80, Strand.POSITIVE)

    assert region.end_on_strand(Strand.POSITIVE) == 80
    assert region.end_on_strand(Strand.NEGATIVE) == 70


def test_with_strand(contig):
    positive = GenomicRegion(contig, 30, 80, Strand.POSITIVE)
    negative = GenomicRegion(contig, 15, 35, Strand.NEGATIVE)

    assert positive.with_strand(Strand.POSITIVE) is positive
    assert negative.with_strand(Strand.NEGATIVE) is negative

    pos_on_neg = positive.with_strand(Strand.NEGATIVE)
    assert pos_on_neg.start == 20
    assert pos_on_neg.end == 70
    assert pos_on_neg.strand == Strand.NEGATIVE

    neg_on_pos = negative.with_strand(Strand.POSITIVE)
    assert neg_on_pos.start == 65
    assert neg_on_pos.end == 85
    assert neg_on_pos.strand == Strand.POSITIVE

