import pytest

from ._genome import Contig, GenomicRegion, Strand, Region


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


class TestOverlaps:

    @pytest.mark.parametrize("a_start, a_end, b_start, b_end, expected",
                             [
                                 (1, 5, 1, 5, True),
                                 (1, 5, 4, 8, True),
                                 (1, 5, 5, 8, False),
                                 (0, 0, 0, 0, True),  # Empty regions do overlap
                                 (0, 0, 0, 1, False),  # But the empty intervals on immediate boundaries do not
                                 (1, 1, 0, 1, False),
                             ])
    def test_overlap_region(self, a_start: int, a_end: int,
                            b_start: int, b_end: int,
                            expected: bool):
        a = Region(a_start, a_end)
        b = Region(b_start, b_end)

        assert a.overlaps_with(b) == expected
        assert b.overlaps_with(a) == expected  # The overlap is transitive!

    @pytest.mark.parametrize("a_strand, a_start, a_end, b_strand, b_start, b_end, expected",
                             [
                                 (Strand.POSITIVE, 1, 2, Strand.POSITIVE, 0, 1, False),
                                 (Strand.POSITIVE, 1, 2, Strand.POSITIVE, 1, 2, True),
                                 (Strand.POSITIVE, 1, 2, Strand.POSITIVE, 2, 3, False),

                                 (Strand.POSITIVE, 1, 3, Strand.POSITIVE, 0, 1, False),
                                 (Strand.POSITIVE, 1, 3, Strand.POSITIVE, 1, 2, True),
                                 (Strand.POSITIVE, 1, 3, Strand.POSITIVE, 2, 3, True),
                                 (Strand.POSITIVE, 1, 3, Strand.POSITIVE, 3, 4, False),

                                 # -----------------------------------------------
                                 #    POS -> 0 1 2 3 4 5
                                 #    NEG <- 5 4 3 2 1 0
                                 (Strand.POSITIVE, 1, 2, Strand.NEGATIVE, 0, 1, False),
                                 (Strand.POSITIVE, 0, 5, Strand.NEGATIVE, 0, 5, True),
                                 (Strand.POSITIVE, 1, 5, Strand.NEGATIVE, 1, 5, True),
                                 (Strand.POSITIVE, 2, 5, Strand.NEGATIVE, 2, 5, True),
                                 (Strand.POSITIVE, 3, 5, Strand.NEGATIVE, 3, 5, False),
                                 (Strand.POSITIVE, 2, 3, Strand.NEGATIVE, 2, 3, True),
                                 (Strand.POSITIVE, 3, 4, Strand.NEGATIVE, 3, 4, False),
                                 (Strand.POSITIVE, 1, 3, Strand.NEGATIVE, 2, 4, True),
                             ])
    def test_overlap_genomic_region(self, a_strand: Strand, a_start: int, a_end: int,
                                    b_strand: Strand, b_start: int, b_end: int,
                                    expected: bool):
        contig = Contig('1', 'GB', 'NC', 'BLA', 5)
        a = GenomicRegion(contig, a_start, a_end, a_strand)
        b = GenomicRegion(contig, b_start, b_end, b_strand)

        assert a.overlaps_with(b) == expected
        assert b.overlaps_with(a) == expected

    def test_overlap_genomic_region__other_contig(self):
        one = Contig('1', 'GB1', 'NC1', 'BLA1', 5)
        two = Contig('2', 'GB2', 'NC2', 'BLA2', 5)
        a = GenomicRegion(one, 0, 5, Strand.POSITIVE)
        b = GenomicRegion(two, 1, 4, Strand.POSITIVE)

        assert not a.overlaps_with(b)


class TestContains:

    @pytest.mark.parametrize("a_start, a_end, b_start, b_end, expected",
                             [
                                 (1, 3, 0, 1, False),
                                 (1, 3, 1, 2, True),
                                 (1, 3, 2, 3, True),
                                 (1, 3, 3, 4, False),
                                 (1, 3, 4, 5, False),

                                 # `b` is bigger or overlaps partially
                                 (1, 3, 0, 5, False),
                                 (1, 3, 0, 3, False),
                                 (1, 3, 0, 2, False),
                                 (1, 3, 2, 4, False),

                                 # An empty region is considered as contained if it lies on either boundary of self.
                                 (0, 1, 0, 0, True),
                                 (0, 1, 1, 1, True),

                             ])
    def test_contains_region(self, a_start: int, a_end: int,
                             b_start: int, b_end: int,
                             expected: bool):
        a = Region(a_start, a_end)
        b = Region(b_start, b_end)

        assert a.contains(b) == expected

    @pytest.mark.parametrize("a_strand, a_start, a_end, b_strand, b_start, b_end, expected",
                             [
                                 (Strand.POSITIVE, 1, 2, Strand.POSITIVE, 0, 1, False),
                                 (Strand.POSITIVE, 1, 2, Strand.POSITIVE, 1, 2, True),
                                 (Strand.POSITIVE, 1, 2, Strand.POSITIVE, 2, 3, False),

                                 (Strand.POSITIVE, 1, 3, Strand.POSITIVE, 0, 1, False),
                                 (Strand.POSITIVE, 1, 3, Strand.POSITIVE, 1, 1, True),
                                 (Strand.POSITIVE, 1, 3, Strand.POSITIVE, 2, 2, True),
                                 (Strand.POSITIVE, 1, 3, Strand.POSITIVE, 3, 4, False),

                                 # Now let's also tweak the strand
                                 # -----------------------------------------------
                                 #    POS -> 0 1 2 3 4 5
                                 #    NEG <- 5 4 3 2 1 0
                                 (Strand.POSITIVE, 1, 2, Strand.NEGATIVE, 0, 1, False),
                                 (Strand.POSITIVE, 0, 5, Strand.NEGATIVE, 0, 5, True),
                                 (Strand.POSITIVE, 1, 5, Strand.NEGATIVE, 1, 5, False),
                                 (Strand.POSITIVE, 2, 5, Strand.NEGATIVE, 2, 5, False),
                                 (Strand.POSITIVE, 2, 3, Strand.NEGATIVE, 2, 3, True),
                                 (Strand.POSITIVE, 3, 4, Strand.NEGATIVE, 3, 4, False),
                                 (Strand.POSITIVE, 1, 3, Strand.NEGATIVE, 2, 4, True),
                             ])
    def test_contains_genomic_region(self, a_strand: Strand, a_start: int, a_end: int,
                                     b_strand: Strand, b_start: int, b_end: int,
                                     expected: bool):
        contig = Contig('1', 'GB', 'NC', 'BLA', 5)
        a = GenomicRegion(contig, a_start, a_end, a_strand)
        b = GenomicRegion(contig, b_start, b_end, b_strand)

        assert a.contains(b) == expected

    def test_contains_genomic_region__other_contig(self):
        one = Contig('1', 'GB1', 'NC1', 'BLA1', 5)
        two = Contig('2', 'GB2', 'NC2', 'BLA2', 5)
        a = GenomicRegion(one, 0, 5, Strand.POSITIVE)
        b = GenomicRegion(two, 1, 4, Strand.POSITIVE)

        assert not a.contains(b)

    @pytest.mark.parametrize("a_start, a_end, pos, expected",
                             [
                                 (2, 3, 0, False),
                                 (2, 3, 1, False),
                                 (2, 3, 2, False),
                                 (2, 3, 3, True),
                                 (2, 3, 4, False),
                                 (2, 3, 5, False),

                                 # Empty region contains no position
                                 (3, 3, 2, False),
                                 (3, 3, 3, False),
                                 (3, 3, 4, False),

                                 # Negative coordinate does not raise
                                 (2, 3, -1, False),
                             ])
    def test_contains_pos(self, a_start, a_end, pos, expected):
        region = Region(a_start, a_end)

        assert region.contains_pos(pos) == expected


class TestDistanceTo:

    @pytest.mark.parametrize("a_start, a_end, b_start, b_end, expected",
                             [

                                 (0, 2, 1, 3, 0),  # Overlapping
                                 (0, 2, 2, 4, 0),  # Adjacent
                                 (0, 2, 4, 6, 2),
                             ])
    def test_distance_to_region(self, a_start: int, a_end: int,
                                b_start: int, b_end: int,
                                expected: int):
        a = Region(a_start, a_end)
        b = Region(b_start, b_end)

        assert a.distance_to(b) == expected
        assert b.distance_to(a) == -expected

    @pytest.mark.parametrize("a_strand, a_start, a_end, b_strand, b_start, b_end, expected",
                             [
                                 # -----------------------------------------------
                                 #    POS ->  0 1 2 3 4 5 6 7 8 9 10
                                 #    NEG <- 10 9 8 7 6 5 4 3 2 1 0
                                 (Strand.POSITIVE, 3, 5, Strand.POSITIVE, 0, 2, -1),
                                 (Strand.POSITIVE, 3, 5, Strand.NEGATIVE, 8, 10, -1),
                                 (Strand.POSITIVE, 3, 5, Strand.POSITIVE, 1, 3, 0),
                                 (Strand.POSITIVE, 3, 5, Strand.NEGATIVE, 7, 9, 0),
                                 (Strand.POSITIVE, 3, 5, Strand.POSITIVE, 2, 4, 0),
                                 (Strand.POSITIVE, 3, 5, Strand.NEGATIVE, 6, 8, 0),
                                 (Strand.POSITIVE, 3, 5, Strand.POSITIVE, 5, 7, 0),
                                 (Strand.POSITIVE, 3, 5, Strand.NEGATIVE, 3, 5, 0),
                                 (Strand.POSITIVE, 3, 5, Strand.POSITIVE, 6, 8, 1),
                                 (Strand.POSITIVE, 3, 5, Strand.NEGATIVE, 2, 4, 1),
                             ])
    def test_distance_to_genomic_region(self, a_strand: Strand, a_start: int, a_end: int,
                                        b_strand: Strand, b_start: int, b_end: int,
                                        expected: int):
        contig = Contig('1', 'GB', 'NC', 'BLA', 10)
        a = GenomicRegion(contig, a_start, a_end, a_strand)
        b = GenomicRegion(contig, b_start, b_end, b_strand)

        assert a.distance_to(b) == expected

    def test_distance_to_genomic_region__other_contig(self):
        one = Contig('1', 'GB1', 'NC1', 'BLA1', 5)
        two = Contig('2', 'GB2', 'NC2', 'BLA2', 5)
        a = GenomicRegion(one, 0, 5, Strand.POSITIVE)
        b = GenomicRegion(two, 1, 4, Strand.POSITIVE)

        with pytest.raises(ValueError) as e:
            a.distance_to(b)
        assert e.value.args == ('Cannot calculate distance between regions on different contigs: 1 <-> 2',)
