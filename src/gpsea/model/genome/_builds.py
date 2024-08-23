import platform
import warnings

from ._genome import Contig, GenomeBuild, GenomeBuildIdentifier

major, minor, patch = platform.python_version_tuple()

if major == '3':
    import importlib.resources
    minor = int(minor)
    if minor < 9:
        # Versions 3.7, 3.8
        def resource_loader(path: str):
            return importlib.resources.open_text('gpsea.model.genome', path)
    else:
        def resource_loader(path: str):
            return importlib.resources.files('gpsea.model.genome').joinpath(path).open()

else:
    warnings.warn(f'Untested Python version v{major}.{minor}.{patch}')


def read_assembly_report(identifier: GenomeBuildIdentifier, path: str) -> GenomeBuild:
    contigs = []
    with resource_loader(path) as fp:
        for line in fp:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            name = fields[0]
            gb_acc = fields[4]
            refseq = fields[6]
            ucsc = fields[9]
            length = int(fields[8])
            contig = Contig(name, gb_acc, refseq, ucsc, length)
            contigs.append(contig)

    return GenomeBuild(identifier, contigs)


GRCh37 = read_assembly_report(GenomeBuildIdentifier('GRCh37', 'p13'), 'GCF_000001405.25_GRCh37.p13_assembly_report.tsv')
"""
The `GRCh37.p13` genomic build.
"""

GRCh38 = read_assembly_report(GenomeBuildIdentifier('GRCh38', 'p13'), 'GCF_000001405.39_GRCh38.p13_assembly_report.tsv')
"""
The `GRCh38.p13` genomic build.
"""
