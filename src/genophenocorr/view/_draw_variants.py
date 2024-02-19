import typing
import os
from itertools import cycle

import numpy as np
import matplotlib.pyplot as plt
import hpotk
from hpotk.validate import ValidationRunner
from hpotk.validate import ObsoleteTermIdsValidator, PhenotypicAbnormalityValidator, AnnotationPropagationValidator

from genophenocorr.model import Cohort
from genophenocorr.preprocessing import configure_caching_patient_creator
from genophenocorr.preprocessing import load_phenopacket_folder
from genophenocorr.model.genome import GRCh38
from genophenocorr.preprocessing import VVTranscriptCoordinateService
from genophenocorr.preprocessing import UniprotProteinMetadataService


#  BASIC DRAWING METHODS
def draw_rectangle(start_x, start_y, end_x, end_y, line_color='black', fill_color=None, line_width=1.0):
    rect = plt.Rectangle((start_x, start_y), end_x - start_x, end_y - start_y, edgecolor=line_color,
                         fill=fill_color is not None, linewidth=line_width, facecolor=fill_color)
    plt.gca().add_patch(rect)


def draw_line(start_x, start_y, end_x, end_y, line_color='black', line_width=1.0):
    plt.plot([start_x, end_x], [start_y, end_y], color=line_color, linewidth=line_width)


def draw_circle(center_x, center_y, radius, line_color='black', fill_color=None, line_width=1.0):
    circle = plt.Circle((center_x, center_y), radius, edgecolor=line_color, fill=fill_color is not None,
                        linewidth=line_width, facecolor=fill_color)
    plt.gca().add_patch(circle)


def draw_string(text, x, y, ha, va, color='black', fontsize=12, rotation=0):
    plt.text(x, y, text, fontsize=fontsize, color=color, ha=ha, va=va, rotation=rotation)


class VariantsVisualizer:
    def __init__(self, transcript_id: str, protein_id: str, cohort: Cohort):
        self.transcript_id = transcript_id
        self.tx_coordinates = VVTranscriptCoordinateService(genome_build=GRCh38).fetch(self.transcript_id)
        self.protein_id = protein_id
        # TODO @ielis: in your notebook you wrote: "Genophenocorr knows that `NM_001032386.2` transcript corresponds to `NP_001027558.1` protein and it can fetch the corresponding metadata:"
        # but at the top you wrote that we should input it manually. Should I take it as a parameter or can i use genophenocorr to figure it out
        self.cohort = cohort
        self.variants = self.cohort.all_variants

        pms = UniprotProteinMetadataService()
        protein_metas = pms.annotate(protein_id)
        assert len(protein_metas) == 1
        self.protein_meta = protein_metas[0]

        self.protein_track_color = '#d3d3d3'
        self.marker_colors = ['red', 'green', 'yellow', 'orange', 'purple']
        self.feature_colors = ['blue', 'green', 'yellow', 'orange']
        self.exon_colors = cycle(['blue', 'lightblue'])
        self.exon_outline = 'black'
        self.axis_color = 'black'

    def _draw_marker(self, x, min_y, max_y, circle_radius, color):
        draw_line(x, min_y, x, max_y, line_color=self.protein_track_color, line_width=0.5)
        draw_circle(x, max_y, circle_radius, line_color=self.protein_track_color, fill_color=color, line_width=0.5)

    def _marker_dim(self, marker_count, protein_track_y_max, marker_length=0.02, marker_radius=0.0025):
        radius = marker_radius + np.sqrt(marker_count - 1) * marker_radius
        length = protein_track_y_max + marker_length + np.sqrt(marker_count - 1) * marker_length
        return radius, length

    def draw_fig(self, limits, protein_limits, markers):
        protein_track_x_min, protein_track_x_max = 0.15, 0.85
        protein_track_y_min, protein_track_y_max = 0.492, 0.508
        font_size = 12
        text_padding = 0.004

        plt.figure(figsize=(20, 20))

        max_x = max(np.max(limits), np.max(protein_limits), np.max(markers))

        # count marker occurrences and remove duplicates
        markers, marker_counts = np.unique(markers, return_counts=True)
        max_marker_count = np.max(marker_counts)

        # normalize into [0, 1], leaving some space on the sides
        preprocess = lambda x: (x / max_x) * (protein_track_x_max - protein_track_x_min) + protein_track_x_min
        protein_limits = preprocess(protein_limits)
        limits = preprocess(limits)
        markers = preprocess(markers)

        # draw the protein track
        draw_rectangle(protein_track_x_min, protein_track_y_min, protein_track_x_max, protein_track_y_max,
                       line_color=self.protein_track_color, fill_color=self.protein_track_color,line_width=2.0)
        # x_axis
        x_axis_y = protein_track_y_min - 0.02
        x_axis_min_x, x_axis_max_x = protein_track_x_min, protein_track_x_max
        big_tick_length, small_tick_length = 0.01, 0.005
        draw_line(x_axis_min_x, x_axis_y, x_axis_max_x, x_axis_y, line_color=self.axis_color, line_width=1.0)  # main line
        draw_line(x_axis_min_x, x_axis_y - big_tick_length, x_axis_min_x, x_axis_y, line_color=self.axis_color,
                  line_width=1.0)  # 0 tick
        draw_string("0", x_axis_min_x, x_axis_y - big_tick_length - text_padding, fontsize=font_size, ha='center',
                    va='top')
        draw_line(x_axis_max_x, x_axis_y - big_tick_length, x_axis_max_x, x_axis_y, line_color=self.axis_color,
                  line_width=1.0)  # max tick
        draw_string(str(max_x), x_axis_max_x, x_axis_y - big_tick_length - text_padding, fontsize=font_size,
                    ha='center', va='top')

        # y_axis
        y_axis_x = protein_track_x_min - 0.02
        y_axis_min_y = protein_track_y_max + 0.01
        _, y_axis_max_y = self._marker_dim(max_marker_count, protein_track_y_max)
        draw_line(y_axis_x, y_axis_min_y, y_axis_x, y_axis_max_y, line_color=self.axis_color, line_width=1.0)
        draw_line(y_axis_x - small_tick_length, y_axis_min_y, y_axis_x, y_axis_min_y, line_color=self.axis_color,
                  line_width=1.0)  # 0 tick
        draw_string("0", y_axis_x - small_tick_length - text_padding, y_axis_min_y, fontsize=font_size, ha='right',
                    va='center')
        draw_line(y_axis_x - small_tick_length, y_axis_max_y, y_axis_x, y_axis_max_y, line_color=self.axis_color,
                  line_width=1.0)  # max tick
        draw_string(str(max_marker_count), y_axis_x - small_tick_length - text_padding, y_axis_max_y,
                    fontsize=font_size, ha='right', va='center')
        draw_string("# Markers", y_axis_x - 0.05, (y_axis_min_y + y_axis_max_y) / 2, fontsize=font_size, ha='center',
                    va='center', rotation=90)  # x axis label

        # draw markers
        marker_y_min = protein_track_y_max
        for marker in markers:
            marker_count = marker_counts[np.where(markers == marker)[0][0]]
            cur_radius, cur_length = self._marker_dim(marker_count, protein_track_y_max)
            self._draw_marker(marker, marker_y_min, cur_length, cur_radius, np.random.choice(self.marker_colors))

        # draw the features
        exon_y_min, exon_y_max = 0.485, 0.515
        for exon_x_min, exon_x_max in limits:
            draw_rectangle(exon_x_min, exon_y_min, exon_x_max, exon_y_max, line_color=self.feature_outline_color,
                           fill_color=np.random.choice(self.feature_colors), line_width=1.0)

        # draw the protein
        protein_y_min, protein_y_max = 0.39, 0.43
        # iterate over pairs in protein_limits
        for protein_x_min, protein_x_max in [protein_limits[i:i + 2] for i in range(len(protein_limits) - 1)]:
            draw_rectangle(protein_x_min, protein_y_min, protein_x_max, protein_y_max, line_color=self.exon_outline_color,
                           fill_color=next(self.exon_colors), line_width=1.0)

        plt.xlim(0, 1)
        plt.ylim(0.3, 0.7)
        plt.gca().set_aspect('equal')
        plt.axis('off')
        plt.show()
