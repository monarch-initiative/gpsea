import random
import typing

from itertools import cycle

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np

from genophenocorr.model import VariantEffect

from ._protein_visualizable import ProteinVisualizable


def round_to_nearest_power_ten(x, base=None):
    if base is None:
        order_of_magnitude = np.floor(np.log10(np.abs(x)))
        base = 10 ** order_of_magnitude
        return (base * np.round(x / base)).astype(int), base
    else:
        return (base * np.round(x / base)).astype(int)


#  BASIC DRAWING METHODS
def draw_rectangle(
        ax: plt.Axes,
        start_x, start_y, end_x, end_y, line_color='black', fill_color=None, line_width=1.0,
):
    rect = plt.Rectangle((start_x, start_y), end_x - start_x, end_y - start_y, edgecolor=line_color,
                         fill=fill_color is not None, linewidth=line_width, facecolor=fill_color)
    ax.add_patch(rect)


def draw_line(
        ax: plt.Axes,
        start_x, start_y, end_x, end_y, line_color='black', line_width=1.0,
):
    ax.plot([start_x, end_x], [start_y, end_y], color=line_color, linewidth=line_width)


def draw_circle(
        ax: plt.Axes,
        center_x, center_y, radius, line_color='black', fill_color=None, line_width=1.0,
):
    circle = plt.Circle((center_x, center_y), radius, edgecolor=line_color, fill=fill_color is not None,
                        linewidth=line_width, facecolor=fill_color)
    ax.add_patch(circle)


def draw_string(
        ax: plt.Axes,
        text, x, y, ha, va, color='black', fontsize=12, rotation=0,
):
    ax.text(x, y, text, fontsize=fontsize, color=color, ha=ha, va=va, rotation=rotation)


def generate_ticks(apprx_n_ticks, min, max):
    tick_step_size, base = round_to_nearest_power_ten((max - min) / apprx_n_ticks)
    ticks = np.array(list(filter(
        lambda x: x < max,
        (min + i * tick_step_size for i in range(1, apprx_n_ticks + 1))
    ))).astype(int)
    return round_to_nearest_power_ten(ticks, base)


class ProteinVisualizer:
    """
    Draw a schema of a protein with variants of the cohort.
    """

    def __init__(self) -> None:
        self.protein_track_color = '#a9a9a9'
        self.transcript_track_color = '#a9a9a9'
        self.marker_colors = {
            VariantEffect.TRANSCRIPT_ABLATION: "#ff0000",
            VariantEffect.SPLICE_ACCEPTOR_VARIANT: "#00ff00",
            VariantEffect.SPLICE_DONOR_VARIANT: "#0000ff",
            VariantEffect.STOP_GAINED: "#ff00ff",
            VariantEffect.FRAMESHIFT_VARIANT: "#ffff00",
            VariantEffect.STOP_LOST: "#00ffff",
            VariantEffect.START_LOST: "#ff9900",
            VariantEffect.TRANSCRIPT_AMPLIFICATION: "#9900ff",
            VariantEffect.INFRAME_INSERTION: "#ff0099",
            VariantEffect.INFRAME_DELETION: "#99ff00",
            VariantEffect.MISSENSE_VARIANT: "#00ff99",
            VariantEffect.PROTEIN_ALTERING_VARIANT: "#990000",
            VariantEffect.SPLICE_REGION_VARIANT: "#009900",
            VariantEffect.SPLICE_DONOR_5TH_BASE_VARIANT: "#009999",
            VariantEffect.SPLICE_DONOR_REGION_VARIANT: "#990099",
            VariantEffect.SPLICE_POLYPYRIMIDINE_TRACT_VARIANT: "#999900",
            VariantEffect.INCOMPLETE_TERMINAL_CODON_VARIANT: "#999999",
            VariantEffect.START_RETAINED_VARIANT: "#ffcc00",
            VariantEffect.STOP_RETAINED_VARIANT: "#ccff00",
            VariantEffect.SYNONYMOUS_VARIANT: "#00ccff",
            VariantEffect.CODING_SEQUENCE_VARIANT: "#ff00cc",
            VariantEffect.MATURE_MIRNA_VARIANT: "#cc00ff",
            VariantEffect.FIVE_PRIME_UTR_VARIANT: "#ff6600",
            VariantEffect.THREE_PRIME_UTR_VARIANT: "#6600ff",
            VariantEffect.NON_CODING_TRANSCRIPT_EXON_VARIANT: "#ff3366",
            VariantEffect.INTRON_VARIANT: "#3366ff",
            VariantEffect.NMD_TRANSCRIPT_VARIANT: "#ffcc99",
            VariantEffect.NON_CODING_TRANSCRIPT_VARIANT: "#cc99ff",
            VariantEffect.UPSTREAM_GENE_VARIANT: "#ff6633",
            VariantEffect.DOWNSTREAM_GENE_VARIANT: "#6633ff",
            VariantEffect.TFBS_ABLATION: "#cc3300",
            VariantEffect.TFBS_AMPLIFICATION: "#ccff66",
            VariantEffect.TF_BINDING_SITE_VARIANT: "#66ccff",
            VariantEffect.REGULATORY_REGION_ABLATION: "#ff3366",
            VariantEffect.REGULATORY_REGION_AMPLIFICATION: "#3366ff",
            VariantEffect.FEATURE_ELONGATION: "#ffcc33",
            VariantEffect.REGULATORY_REGION_VARIANT: "#ccff33",
            VariantEffect.FEATURE_TRUNCATION: "#33ccff",
            VariantEffect.INTERGENIC_VARIANT: "#ff0033",
            VariantEffect.SEQUENCE_VARIANT: "#33ff00",
        }
        mycolors = [m for m in mcolors.CSS4_COLORS.keys() if "grey" not in m and "white" not in m]
        random.shuffle(mycolors)
        self._available_colors = mycolors
        # self._color_idx = 0
        self.feature_outline_color = 'black'
        self.exon_colors = cycle(['blue', 'lightblue'])
        self.exon_outline_color = 'black'
        self.axis_color = 'black'

    def _draw_marker(
            self,
            ax: plt.Axes,
            x_start, x_end, min_y, max_y, circle_radius, color,
    ):
        """
        Draw a lollipop representing a variant and the number of counts for the variant effect type
        currently putting marker in the middle of start and end, can change this later
        """
        x = (x_start + x_end) / 2
        draw_line(ax, x, min_y, x, max_y - circle_radius, line_color=self.protein_track_color, line_width=0.5)
        draw_circle(ax, x, max_y, circle_radius, line_color=self.protein_track_color, fill_color=color, line_width=0.5)

    @staticmethod
    def _marker_dim(marker_count, protein_track_y_max, marker_length=0.02, marker_radius=0.0025):
        radius = marker_radius + np.sqrt(marker_count - 1) * marker_radius
        length = protein_track_y_max + marker_length + np.sqrt(marker_count - 1) * marker_length
        return radius, length

    def draw_fig(
            self,
            pvis: ProteinVisualizable,
            ax: typing.Optional[plt.Axes] = None,
            labeling_method: typing.Literal['abbreviate', 'enumerate'] = 'abbreviate',
            legend_x: float = 0.87,
    ) -> typing.Optional[plt.Axes]:
        """
        Visualize the cohort variants on a protein diagram.

        By default, the legend is drawn to the right of the figure to avoid overlap between the variant markers
        and the legend.

        Args:
             pvis: :class:`ProteinVisualizable` with information about the transcript coordinates, protein metadata,
               and the cohort for plotting
             ax: a Matplotlib :class:`plt.Axes` to plot on or `None` if a new `Axes` should be created
             labeling_method: the strategy for generating labels.
               Valid values of labeling_method are `{'abbreviate', 'enumerate'}`
             legend_x: left x coordinate of the legend bounding box
        Returns:
            `None` if an :class:`plt.Axes` was provided via `ax` argument
            or an `Axes` created by the visualizer if `ax` was `None`.
        """
        if ax is None:
            should_return_ax = True
            _, ax = plt.subplots(figsize=(20, 20))
        else:
            should_return_ax = False
        min_aa_pos = 1
        max_aa_pos = pvis.protein_length
        feature_limits = list()
        for i in range(len(pvis.protein_feature_ends)):
            feature_limits.append((pvis.protein_feature_starts[i], pvis.protein_feature_ends[i]))
        feature_limits = np.array(feature_limits)
        color_idx = 0
        protein_feature_colors = dict()
        seen_features = set()
        for feature_name in pvis.protein_feature_names:
            if feature_name in seen_features:
                continue
            seen_features.add(feature_name)
            color = self._available_colors[color_idx]
            protein_feature_colors[feature_name] = color
            color_idx += 1
        feature_colors = [protein_feature_colors[ft] for ft in pvis.protein_feature_names]
        feature_names = pvis.protein_feature_names
        variant_locations = pvis.variant_locations
        # The following has the variant positions (in variant_locations_counted_absolute)
        # and the number of occurrences of each position (in marker counts)
        variant_locations_counted_absolute = pvis.variant_locations_counted_absolute
        marker_counts = pvis.marker_counts
        variant_effect_colors = []
        for vl in variant_locations_counted_absolute:
            i = np.where(variant_locations == vl)[0][0]  # find index of unique variant loc in all locs to find effect
            effect = pvis.variant_effects[i]
            variant_effect_colors.append(self.marker_colors[effect])

        protein_track_x_min, protein_track_x_max = 0.15, 0.85
        protein_track_y_min, protein_track_y_max = 0.492, 0.508
        font_size = 12
        text_padding = 0.004

        max_marker_count = np.max(marker_counts)

        x_ticks = generate_ticks(apprx_n_ticks=6, min=min_aa_pos, max=max_aa_pos)
        y_ticks = generate_ticks(apprx_n_ticks=5, min=0, max=max_marker_count)

        # normalize into [0, 1], leaving some space on the sides
        def preprocess(
                absolute: np.ndarray,
                min_absolute=min_aa_pos, max_absolute=max_aa_pos,
                min_relative=protein_track_x_min, max_relative=protein_track_x_max,
        ) -> np.ndarray:
            if absolute.any() > max_absolute:
                # should never happen
                print(f"[ERROR _variant_drawer] x={absolute} but max={max_absolute}")
            shifted_to_0_1 = ((absolute - min_absolute) / (max_absolute - min_absolute))
            relative_scale = (max_relative - min_relative)
            return shifted_to_0_1 * relative_scale + min_relative

        feature_limits_relative = preprocess(feature_limits)
        variant_locations_relative = preprocess(variant_locations_counted_absolute)
        x_ticks_relative = preprocess(x_ticks)

        # draw the tracks
        draw_rectangle(
            ax,
            protein_track_x_min, protein_track_y_min, protein_track_x_max, protein_track_y_max,
            line_color=self.protein_track_color, fill_color=self.protein_track_color, line_width=2.0,
        )
        # x_axis
        x_axis_y = protein_track_y_min - 0.02
        x_axis_min_x, x_axis_max_x = protein_track_x_min, protein_track_x_max
        big_tick_length, small_tick_length = 0.015, 0.005
        draw_line(  # main line
            ax,
            x_axis_min_x, x_axis_y, x_axis_max_x, x_axis_y,
            line_color=self.axis_color, line_width=1.0,
        )
        draw_string(
            ax,
            str(min_aa_pos), x_axis_min_x, x_axis_y - big_tick_length - text_padding,
            fontsize=font_size, ha='center', va='top',
        )
        draw_string(
            ax,
            str(max_aa_pos), x_axis_max_x, x_axis_y - big_tick_length - text_padding,
            fontsize=font_size, ha='center', va='top',
        )
        # draw x ticks
        draw_line(  # max tick
            ax,
            x_axis_max_x, x_axis_y - big_tick_length, x_axis_max_x, x_axis_y,
            line_color=self.axis_color, line_width=1.0,
        )
        draw_line(  # minimum tick
            ax,
            x_axis_min_x, x_axis_y - big_tick_length, x_axis_min_x, x_axis_y,
            line_color=self.axis_color, line_width=1.0,
        )
        for x_tick_relative, x_tick_absolute in zip(x_ticks_relative, x_ticks):
            draw_line(
                ax,
                x_tick_relative, x_axis_y - small_tick_length, x_tick_relative, x_axis_y,
                line_color=self.axis_color, line_width=1.0,
            )
            draw_string(
                ax,
                str(x_tick_absolute), x_tick_relative, x_axis_y - small_tick_length - text_padding,
                fontsize=font_size, ha='center', va='top',
            )

        # y_axis
        y_axis_x = protein_track_x_min - 0.02
        y_axis_min_y = protein_track_y_max + 0.01
        _, y_axis_max_y = ProteinVisualizer._marker_dim(max_marker_count, protein_track_y_max)
        y_ticks_relative = preprocess(
            y_ticks,
            min_absolute=0, max_absolute=max_marker_count, min_relative=y_axis_min_y, max_relative=y_axis_max_y,
        )
        draw_line(
            ax,
            y_axis_x, y_axis_min_y, y_axis_x, y_axis_max_y,
            line_color=self.axis_color, line_width=1.0,
        )
        draw_string(
            ax,
            "0", y_axis_x - small_tick_length - text_padding, y_axis_min_y,
            fontsize=font_size, ha='right', va='center',
        )
        draw_string(
            ax,
            str(max_marker_count), y_axis_x - small_tick_length - text_padding, y_axis_max_y,
            fontsize=font_size, ha='right', va='center',
        )
        draw_string(  # x axis label
            ax,
            "# Variants", y_axis_x - 0.05, (y_axis_min_y + y_axis_max_y) / 2,
            fontsize=font_size, ha='center', va='center', rotation=90,
        )
        # draw y ticks
        draw_line(  # 0 tick
            ax,
            y_axis_x - small_tick_length, y_axis_min_y, y_axis_x, y_axis_min_y,
            line_color=self.axis_color, line_width=1.0,
        )
        draw_line(  # max tick
            ax,
            y_axis_x - small_tick_length, y_axis_max_y, y_axis_x, y_axis_max_y,
            line_color=self.axis_color, line_width=1.0,
        )
        for y_tick_relative, y_tick_absolute in zip(y_ticks_relative, y_ticks):
            draw_line(
                ax,
                y_axis_x - small_tick_length, y_tick_relative, y_axis_x, y_tick_relative,
                line_color=self.axis_color, line_width=1.0,
            )
            draw_string(
                ax,
                str(y_tick_absolute), y_axis_x - small_tick_length - text_padding, y_tick_relative,
                fontsize=font_size, ha='right', va='center',
            )

        # draw variants
        marker_y_min = protein_track_y_max
        for marker, marker_color in zip(variant_locations_relative, variant_effect_colors):
            marker_count = marker_counts[np.where(variant_locations_relative == marker)[0][0]]
            cur_radius, cur_length = ProteinVisualizer._marker_dim(marker_count, protein_track_y_max)
            x_start, x_end = marker, marker #  WAS  marker[0], marker[1]
            self._draw_marker(ax, x_start, x_end, marker_y_min, cur_length, cur_radius, marker_color)

        # draw the features (protein track)
        feature_y_min, feature_y_max = 0.485, 0.515

        unique_feature_names = list(set(feature_names))

        if labeling_method == 'enumerate':
            ascii_A = 65
            labels = {fn: chr(ascii_A + i) for i, fn in enumerate(unique_feature_names)}

        elif labeling_method == 'abbreviate':
            labels = {fn: fn[0:5] for fn in unique_feature_names}
        else:
            raise ValueError(f'Unsupported labeling method {labeling_method}')

        for feature_x, feature_color, feature_name in zip(feature_limits_relative, feature_colors, feature_names):
            feature_x_min, feature_x_max = feature_x
            draw_rectangle(
                ax,
                feature_x_min, feature_y_min, feature_x_max, feature_y_max,
                line_color=self.feature_outline_color, fill_color=feature_color, line_width=1.0,
            )
            if (feature_x_max - feature_x_min) <= 0.03:  # too small to display name
                draw_string(
                    ax, labels[feature_name],
                    0.05 * (feature_x_max - feature_x_min) + feature_x_min,
                    0.55 * (feature_y_max - feature_y_min) + feature_y_min,
                    ha="left", va="center", rotation=90, color='black', fontsize=8,
                )
            elif (feature_x_max - feature_x_min) <= 0.005:  # too small even to draw vertical string
                # TODO @ielis: How to display name here?
                pass
            else:
                draw_string(
                    ax, labels[feature_name],
                    0.2 * (feature_x_max - feature_x_min) + feature_x_min,
                    0.4 * (feature_y_max - feature_y_min) + feature_y_min,
                    ha="left", va="center", color='black',
                )

        # draw legend
        n_unique_features = len(unique_feature_names)
        color_box_x_dim = 0.01
        color_box_y_dim = 0.01
        if labeling_method == 'abbreviate':
            color_box_x_dim *= 3.5
        row_spacing = 0.005
        legend_width = 0.2 + color_box_x_dim
        legend_max_y = 0.75
        legend_min_y = legend_max_y - (n_unique_features + 1) * row_spacing - n_unique_features * color_box_y_dim
        legend_min_x = legend_x
        legend_max_x = legend_min_x + legend_width

        # legend box
        draw_rectangle(ax, legend_min_x, legend_min_y, legend_max_x, legend_max_y, 'black')
        for i, feature_name in enumerate(unique_feature_names):
            # colored box
            color_box_min_x = legend_min_x + row_spacing
            color_box_max_x = color_box_min_x + color_box_x_dim
            color_box_max_y = legend_max_y - (i + 1) * row_spacing - i * color_box_y_dim
            color_box_min_y = color_box_max_y - color_box_y_dim
            draw_rectangle(
                ax, color_box_min_x, color_box_min_y, color_box_max_x, color_box_max_y,
                line_color='black', fill_color=protein_feature_colors[feature_name],
            )
            # label in colored box (same as on x axis in the boxes)
            draw_string(
                ax, labels[feature_name], color_box_min_x + 0.002, color_box_min_y + 0.005,
                fontsize=10, ha="left", va="center", color='black',
            )
            # full feature name
            draw_string(
                ax, feature_name,
                color_box_max_x + 0.005, color_box_min_y + 0.005,
                fontsize=12, ha="left", va="center", color='black',
            )

        ax.set(
            xlim=(0, max(1.0, legend_x + legend_width + 0.02)),
            ylim=(0.3, 0.75),
            aspect='equal',
            axis='off',
            title=f'transcript: {pvis.transcript_id}, '
                  f'protein: {pvis.protein_id}, '
                  f'protein name {pvis.protein_metadata.label}',
        )

        if should_return_ax:
            return ax
