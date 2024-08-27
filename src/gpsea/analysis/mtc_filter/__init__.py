"""
The `gpsea.analysis.mtc_filter` provides the strategies for reducing multiple testing burden
by pre-filtering the phenotype terms to test in advance.

See :ref:`MTC filters <mtc-filters>` section for more info.
"""

from ._impl import PhenotypeMtcFilter, PhenotypeMtcResult
from ._impl import UseAllTermsMtcFilter, SpecifiedTermsMtcFilter, HpoMtcFilter

__all__ = [
    'PhenotypeMtcFilter', 'PhenotypeMtcResult',
    'UseAllTermsMtcFilter', 'SpecifiedTermsMtcFilter', 'HpoMtcFilter',
]
