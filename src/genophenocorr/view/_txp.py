import typing
import os
import math
import numpy as np
from matplotlib import pyplot as plt
from shapely.geometry.polygon import Polygon
from descartes import PolygonPatch
from genophenocorr.model import Variant, TranscriptCoordinates, ProteinMetadata


class VariantTranscriptProteinArtist:
    """
    `VariantTranscriptProteinArtist` creates a graphic to show distributions of variants across a provided transcript
    and its protein product.
    """

    def __init__(self, width):
        # Global variables for converting genomic coordinates to plot coordinates and vice versa
        self._margin_width = 20
        self._figure_size = 10
        self._min = 0
        self._max = width + 2*self._margin_width
        self._span = self._max - self._min
        self._factor = self._figure_size / self._span

    def draw_variants(self, variants: typing.Sequence[Variant],
                      tx: TranscriptCoordinates,
                      verbose: bool = True):
        # TODO - Peter implement
        if verbose:
            print("[INFO] Creating visualization ...")
        if not isinstance(variants, list):
            #raise ValueError(f"variants argument must be a list but was {type(variants)}")
            print(f"variants argument must be a list but was {type(variants)}")
            # TODO -- Should this be a set or a list?
            variants = list(variants)
        if len(variants) == 0:
            raise ValueError("variant objects was an empty list")
        if not isinstance(variants[0], Variant):
            raise ValueError(f"variant list must contain Variant objects but contained {type(variants[0])}")
        if not isinstance(tx, TranscriptCoordinates):
            raise ValueError(f"tx arguments must be of type TranscriptCoordinates but was {type(tx)}")
        figure_height = 6
        figure_width = 2 * figure_height
        if verbose:
            print(f"\t[INFO] Will plot transcripts {tx}")
            print(f"\t[INFO] Will plot {len(variants)} variants")
        # Create figure
        fig, ax = plt.subplots(1, figsize=(figure_width, figure_height))
        transcript_id = tx.id
        plot_title = f"Variants in {transcript_id}"
        ax.set_title(plot_title, loc='left', fontsize='x-large')
        ax.set_xlabel('TODO X Label', labelpad=12, fontsize=14.5)
        ax.set_ylabel('TODO Y Label', labelpad=12, fontsize=14.5)







