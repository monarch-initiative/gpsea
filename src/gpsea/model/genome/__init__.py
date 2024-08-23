"""
`gpsea.model.genome` provides data classes to model :class:`GenomeBuild`, :class:`Contig`,
and (genomic) regions (:class:`Region`, :class:`GenomicRegion`).

The classes can do basic region arithmetics such as finding intersections, overlaps, and distances between regions.
Genomic regions are transposable - they can flip the coordinates between DNA strands.

The module provides *GRCh37.p13* and *GRCh38.p13*, the two most commonly used human genome builds.

The classes are largely a port of `Svart <https://github.com/exomiser/svart>`_ library.
"""

from ._builds import GRCh37, GRCh38
from ._genome import Contig, GenomeBuild, GenomeBuildIdentifier, Strand, Stranded, Transposable, GenomicRegion, Region
from ._genome import transpose_coordinate

__all__ = [
    "GenomeBuild", "Contig", "GenomeBuildIdentifier", "Region", "GenomicRegion",
    "Strand", "Stranded", "Transposable",
    "transpose_coordinate",
    "GRCh37", "GRCh38",
]
