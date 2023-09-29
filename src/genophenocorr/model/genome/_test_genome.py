from ._genome import Contig, GenomicRegion


def test_start_on_strand():
    contig = Contig('1', 'GB_BLA', 'NC_BLA', 'UCSC_BLA', 100)
    region = GenomicRegion(contig, 30, 80, True)

    assert region.start_on_strand(True) == 30
    assert region.start_on_strand(False) == 20

    assert region.end_on_strand(True) == 80
    assert region.end_on_strand(False) == 70
