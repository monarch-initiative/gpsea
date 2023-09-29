import pytest

from ._builds import read_assembly_report, GRCH37, GRCH38


def test_read_assembly_report():
    build = read_assembly_report('GRCh37.p13', 'GCF_000001405.25_GRCh37.p13_assembly_report.tsv')
    assert build is not None


@pytest.mark.parametrize("name, genbank, refseq, ucsc, length",
                         [
                             ('1', 'CM000663.1', 'NC_000001.10', 'chr1', 249_250_621),
                             ('X', 'CM000685.1', 'NC_000023.10', 'chrX', 155_270_560),
                          ])
def test_hg19(name, genbank, refseq, ucsc, length):
    assert GRCH37.identifier == 'GRCh37.p13'
    contig = GRCH37.contig_by_name(name)
    assert contig.name == name
    assert contig.genbank_acc == genbank
    assert contig.refseq_name == refseq
    assert contig.ucsc_name == ucsc
    assert len(contig) == length


@pytest.mark.parametrize("name, genbank, refseq, ucsc, length",
                         [
                             ('1', 'CM000663.2', 'NC_000001.11', 'chr1', 248_956_422),
                             ('X', 'CM000685.2', 'NC_000023.11', 'chrX', 156_040_895),
                          ])
def test_hg38(name, genbank, refseq, ucsc, length):
    assert GRCH38.identifier == 'GRCh38.p13'
    contig = GRCH38.contig_by_name(name)
    assert contig.name == name
    assert contig.genbank_acc == genbank
    assert contig.refseq_name == refseq
    assert contig.ucsc_name == ucsc
    assert len(contig) == length
