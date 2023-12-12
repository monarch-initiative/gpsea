import pytest

from ._tx import TranscriptCoordinates
from .genome import Contig, GenomicRegion, Strand


class TestTranscriptCoordinates:

    @pytest.fixture
    def contig(self) -> Contig:
        return Contig('a', 'b', 'c', 'd', 100)

    @pytest.fixture
    def simple_tx(self, contig: Contig) -> TranscriptCoordinates:
        return TranscriptCoordinates('some_id', region=GenomicRegion(contig, 20, 80, Strand.POSITIVE), exons=(
            GenomicRegion(contig, 20, 30, Strand.POSITIVE), GenomicRegion(contig, 50, 60, Strand.POSITIVE),
            GenomicRegion(contig, 70, 80, Strand.POSITIVE)), cds_start=25, cds_end=76)

    @pytest.fixture
    def noncoding_exons_tx(self, contig: Contig) -> TranscriptCoordinates:
        """
        Transcript coordinates with the same coding regions as above, but with extra non-coding exons for UTRs.
        """
        return TranscriptCoordinates('some_id', region=GenomicRegion(contig, 5, 95, Strand.POSITIVE), exons=(
            GenomicRegion(contig, 5, 10, Strand.POSITIVE), GenomicRegion(contig, 20, 30, Strand.POSITIVE),
            GenomicRegion(contig, 50, 60, Strand.POSITIVE), GenomicRegion(contig, 70, 80, Strand.POSITIVE),
            GenomicRegion(contig, 88, 95, Strand.POSITIVE),), cds_start=25, cds_end=76)

    def test_simple_tx_compute_n_codons(self, simple_tx: TranscriptCoordinates):
        assert simple_tx.get_coding_base_count() == 18
        assert simple_tx.get_codon_count() == 6

    def test_tx_w_noncoding_exons__compute_n_codons(self, noncoding_exons_tx: TranscriptCoordinates):
        assert noncoding_exons_tx.get_coding_base_count() == 18
        assert noncoding_exons_tx.get_codon_count() == 6

    def test_five_prime_utrs(self, simple_tx: TranscriptCoordinates, noncoding_exons_tx: TranscriptCoordinates):
        utrs = simple_tx.get_five_prime_utrs()
        assert len(utrs) == 1
        utr = utrs[0]
        assert utr.start == 20
        assert utr.end == 25

        utrs = noncoding_exons_tx.get_five_prime_utrs()
        assert len(utrs) == 2

        assert utrs[0].start == 5
        assert utrs[0].end == 10
        assert utrs[1].start == 20
        assert utrs[1].end == 25

    def test_three_prime_utrs(self, simple_tx: TranscriptCoordinates, noncoding_exons_tx: TranscriptCoordinates):
        utrs = simple_tx.get_three_prime_utrs()
        assert len(utrs) == 1
        utr = utrs[0]
        assert utr.start == 76
        assert utr.end == 80

        utrs = noncoding_exons_tx.get_three_prime_utrs()
        assert len(utrs) == 2
        assert utrs[0].start == 76
        assert utrs[0].end == 80
        assert utrs[1].start == 88
        assert utrs[1].end == 95
