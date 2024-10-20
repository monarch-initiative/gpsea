"""
The `gpsea.analysis.temporal.endpoint` package provides endpoints
for comparing survivals between genotype groups.
"""

from ._impl import death, disease_onset, hpo_onset

__all__ = [
    "death",
    "disease_onset",
    "hpo_onset",
]
