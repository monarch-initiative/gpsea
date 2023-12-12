import pytest

from ._tx import TranscriptCoordinates
from .genome import Contig, GenomicRegion, Strand


class TestTranscriptCoordinates:

    @pytest.fixture
    def contig(self) -> Contig:
        return Contig('a', 'b', 'c', 'd', 100)

    def test_simple_tx_compute_n_codons(self, contig: Contig):
        tx = TranscriptCoordinates('some_id', region=GenomicRegion(contig, 20, 80, Strand.POSITIVE), exons=(
            GenomicRegion(contig, 20, 30, Strand.POSITIVE), GenomicRegion(contig, 50, 60, Strand.POSITIVE),
            GenomicRegion(contig, 70, 80, Strand.POSITIVE)), cds_start=25, cds_end=76)
        assert tx.get_codon_count() == 6

    def test_tx_w_noncoding_exons__compute_n_codons(self, contig: Contig):
        """
        Transcript coordinates with the same coding regions as above, but with extra non-coding exons for UTRs.
        """
        tx = TranscriptCoordinates('some_id', region=GenomicRegion(contig, 5, 95, Strand.POSITIVE), exons=(
            GenomicRegion(contig, 5, 10, Strand.POSITIVE), GenomicRegion(contig, 20, 30, Strand.POSITIVE),
            GenomicRegion(contig, 50, 60, Strand.POSITIVE), GenomicRegion(contig, 70, 80, Strand.POSITIVE),
            GenomicRegion(contig, 88, 95, Strand.POSITIVE),), cds_start=25, cds_end=76)
        assert tx.get_codon_count() == 6
