import random
from typing import Union, Tuple, List, Optional, Literal
from dataclasses import dataclass

from itertools import cycle

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np

from genophenocorr.model import VariantEffect

from ._protein_visualizable import ProteinVisualizable


class ProteinVisualizer:
    """
    Draw a schema of a protein with variants of the cohort.
    """

    def __init__(self) -> None:
        self.marker_colors = {
            VariantEffect.TRANSCRIPT_ABLATION: "#ff0000",
            VariantEffect.SPLICE_ACCEPTOR_VARIANT: "#00ff00",
            VariantEffect.SPLICE_DONOR_VARIANT: "#ff0099",
            VariantEffect.STOP_GAINED: "#ffcc00",
            VariantEffect.FRAMESHIFT_VARIANT: "#990099",
            VariantEffect.STOP_LOST: "#00ffff",
            VariantEffect.START_LOST: "#ff9900",
            VariantEffect.TRANSCRIPT_AMPLIFICATION: "#9900ff",
            VariantEffect.INFRAME_INSERTION: "#0000ff",
            VariantEffect.INFRAME_DELETION: "#990000",
            VariantEffect.MISSENSE_VARIANT: "#ff0033",
            VariantEffect.PROTEIN_ALTERING_VARIANT: "#99ff00",
            VariantEffect.SPLICE_REGION_VARIANT: "#009900",
            VariantEffect.SPLICE_DONOR_5TH_BASE_VARIANT: "#009999",
            VariantEffect.SPLICE_DONOR_REGION_VARIANT: "#ffff00",
            VariantEffect.SPLICE_POLYPYRIMIDINE_TRACT_VARIANT: "#999900",
            VariantEffect.INCOMPLETE_TERMINAL_CODON_VARIANT: "#999999",
            VariantEffect.START_RETAINED_VARIANT: "#00ff99",
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
            VariantEffect.INTERGENIC_VARIANT: "#ff00ff",
            VariantEffect.SEQUENCE_VARIANT: "#33ff00",
        }
        # colors
        mycolors = [m for m in mcolors.CSS4_COLORS.keys() if "grey" not in m and "white" not in m]
        random.seed(42)
        random.shuffle(mycolors)
        self._available_colors = mycolors
        self.feature_outline_color = 'black'
        self.exon_colors = cycle(['blue', 'lightblue'])
        self.exon_outline_color = 'black'
        self.axis_color = 'black'
        self.protein_track_color = '#a9a9a9'
        self.transcript_track_color = '#a9a9a9'
        # plot options
        self.protein_track_x_min = 0.15
        self.protein_track_x_max = 0.85
        self.protein_track_y_min = 0.492
        self.protein_track_y_max = 0.508
        self.font_size = 12
        self.text_padding = 0.004
        self.feature_y_min = 0.485
        self.feature_y_max = 0.515
        # legend
        self.color_box_x_dim = 0.01
        self.color_box_y_dim = 0.01
        self.row_spacing = 0.005
        self.legend1_max_y = 0.75
        self.color_circle_radius = 0.005
        self.legend1_min_x = 0.87
        self.legend2_min_x = 0.1
        self.legend2_max_x = 0.3
        self.legend2_max_y = 0.75

    def draw_fig(
            self,
            pvis: ProteinVisualizable,
            ax: Optional[plt.Axes] = None,
            labeling_method: Literal['abbreviate', 'enumerate'] = 'abbreviate'
    ) -> Optional[plt.Axes]:
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
        Returns:
            `None` if an :class:`plt.Axes` was provided via `ax` argument
            or an `Axes` created by the visualizer if `ax` was `None`.
        """
        if ax is None:
            should_return_ax = True
            _, ax = plt.subplots(figsize=(20, 20))
        else:
            should_return_ax = False

        # STATE
        feature_handler = DrawableProteinFeatureHandler(pvis, labeling_method, self._available_colors)

        marker_counts, variant_locations_counted_absolute, variant_effect_colors = (
            generate_variant_markers(pvis, self.marker_colors))

        x_ticks = generate_ticks(apprx_n_ticks=6, min=1, max=pvis.protein_length)
        y_ticks = generate_ticks(apprx_n_ticks=5, min=0, max=np.max(marker_counts))

        # normalize into [0, 1], leaving some space on the sides
        for f in feature_handler.features:
            cur_limits = f.min_pos_abs, f.max_pos_abs
            f.min_pos_plotting, f.max_pos_plotting = translate_to_ax_coordinates(
                np.array(cur_limits), min_absolute=1, max_absolute=pvis.protein_length,
                min_relative=self.protein_track_x_min, max_relative=self.protein_track_x_max
            )

        variant_locations_relative = translate_to_ax_coordinates(variant_locations_counted_absolute,
                                                                 min_absolute=1, max_absolute=pvis.protein_length,
                                                                 min_relative=self.protein_track_x_min,
                                                                 max_relative=self.protein_track_x_max, clip=True)
        x_ticks_relative = translate_to_ax_coordinates(x_ticks, min_absolute=1,
                                                       max_absolute=pvis.protein_length,
                                                       min_relative=self.protein_track_x_min,
                                                       max_relative=self.protein_track_x_max)

        # PLOTTING
        draw_axes(ax,
                  x_ticks, x_ticks_relative, y_ticks,
                  np.max(marker_counts), 1, pvis.protein_length,
                  self.protein_track_x_min, self.protein_track_x_max,
                  self.protein_track_y_min, self.protein_track_y_max,
                  self.font_size, self.text_padding,
                  self.axis_color, self.protein_track_color
                  )

        draw_variants(ax, variant_locations_relative, variant_effect_colors, marker_counts, self.protein_track_y_max)

        feature_handler.draw_features(ax, self.feature_y_min, self.feature_y_max, self.feature_outline_color)

        legend1_width = draw_legends(ax, feature_handler, pvis,
                                     self.color_box_x_dim, self.color_box_y_dim, self.color_circle_radius,
                                     self.row_spacing,
                                     self.legend1_min_x, self.legend1_max_y,
                                     self.legend2_min_x, self.legend2_max_y,
                                     variant_effect_colors, self.marker_colors,
                                     labeling_method, )

        ax.set(
            xlim=(0, max(1.0, self.legend1_min_x + legend1_width + 0.02)),
            ylim=(0.3, 0.75),
            aspect='equal',
            title=f'{pvis.protein_metadata.label}\ntranscript: {pvis.transcript_id}, '
                  f'protein: {pvis.protein_id}',
        )

        ax.axis('off')

        if should_return_ax:
            return ax


@dataclass(slots=True)
class DrawableProteinFeature:
    name: str
    min_pos_abs: Union[int, float]
    max_pos_abs: Union[int, float]
    label: str
    color: str
    min_pos_plotting: float
    max_pos_plotting: float

    def draw(self, ax: plt.Axes, y_min: float, y_max: float, feature_outline_color: str):
        draw_rectangle(
            ax,
            self.min_pos_plotting, y_min, self.max_pos_plotting, y_max,
            line_color=feature_outline_color, fill_color=self.color, line_width=1.0,
        )
        # too small to display horizontally, so display vertically
        if (self.max_pos_plotting - self.min_pos_plotting) <= 0.03:
            draw_string(
                ax, self.label,
                0.05 * (self.max_pos_plotting - self.min_pos_plotting) + self.min_pos_plotting,
                0.55 * (y_max - y_min) + y_min,
                ha="left", va="center", rotation=90, color='black', fontsize=8,
            )
        else:
            draw_string(
                ax, self.label,
                0.2 * (self.max_pos_plotting - self.min_pos_plotting) + self.min_pos_plotting,
                0.4 * (y_max - y_min) + y_min,
                ha="left", va="center", color='black',
            )


class DrawableProteinFeatureHandler:
    def __init__(self, pvis: ProteinVisualizable, labeling_method: str, colors: List[str]):
        self.pvis = pvis
        self.labeling_method = labeling_method
        self._available_colors = colors

        self.cleaned_unique_feature_names, self.mapping_all2cleaned, self.labels, self.colors = self._generate_labels()
        self.features = self._generate_features()

    def _generate_labels(self):
        # aggregate similar feature names into one category
        unique_feature_names = list(set(self.pvis.protein_feature_names))
        cleaned_unique_feature_names = set()
        mapping_all2cleaned = dict()
        for feature_name in unique_feature_names:
            # remove digits from feature name
            cleaned_feature_name = str(''.join(char for char in feature_name if not char.isdigit()))
            cleaned_unique_feature_names.add(cleaned_feature_name)
            mapping_all2cleaned[feature_name] = cleaned_feature_name
        # generate labels for features
        if self.labeling_method == 'enumerate':
            ascii_capital_a = ord('A')
            labels = {fn: chr(ascii_capital_a + i) for i, fn in enumerate(cleaned_unique_feature_names)}

        elif self.labeling_method == 'abbreviate':
            labels = {feature_name: feature_name[0:5] for feature_name in cleaned_unique_feature_names}
        else:
            raise ValueError(f'Unsupported labeling method {self.labeling_method}')

        cleaned_unique_feature_names = list(cleaned_unique_feature_names)
        cleaned_unique_feature_names.sort()

        colors = assign_colors(cleaned_unique_feature_names, self._available_colors)

        return cleaned_unique_feature_names, mapping_all2cleaned, labels, colors

    def _generate_features(self):
        features = list()
        for i in range(len(self.pvis.protein_feature_ends)):
            features.append(DrawableProteinFeature(
                min_pos_abs=self.pvis.protein_feature_starts[i],
                max_pos_abs=self.pvis.protein_feature_ends[i],
                name=self.pvis.protein_feature_names[i],
                label=self.labels[self.mapping_all2cleaned[self.pvis.protein_feature_names[i]]],
                color=self.colors[self.mapping_all2cleaned[self.pvis.protein_feature_names[i]]],
                min_pos_plotting=-1.0, max_pos_plotting=-1.0  # will be set later in draw_fig(), when plot limits known
            ))

        return features

    def draw_features(self, ax: plt.Axes, y_min: float, y_max: float, feature_outline_color: str):
        for f in self.features:
            f.draw(ax, y_min, y_max, feature_outline_color)


@dataclass(slots=True)
class DrawableProteinVariant:
    name: str
    pos_abs: Union[int, float]
    label: str
    color: str
    pos_plotting: float
    count: int

    def draw(self, ax: plt.Axes, y_min: float, y_max: float, feature_outline_color: str):
        # TODO: implement
        pass


class DrawableProteinVariantHandler:
    def __init__(self):
        # TODO: implement
        pass

    def draw_variants(self):
        # TODO: implement
        pass


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


def marker_dim(marker_count, protein_track_y_max, marker_length=0.02, marker_radius=0.0025):
    radius = marker_radius + np.sqrt(marker_count - 1) * marker_radius
    length = protein_track_y_max + marker_length + np.sqrt(marker_count - 1) * marker_length
    return radius, length


def round_to_nearest_power_ten(x, base=None):
    if base is None:
        order_of_magnitude = np.floor(np.log10(np.abs(x)))
        base = 10 ** order_of_magnitude
        return (base * np.round(x / base)).astype(int), base
    else:
        return (base * np.round(x / base)).astype(int)


def generate_ticks(apprx_n_ticks, min, max):
    tick_step_size, base = round_to_nearest_power_ten((max - min) / apprx_n_ticks)
    ticks = np.array(list(filter(
        lambda x: x < max,
        (min + i * tick_step_size for i in range(1, apprx_n_ticks + 1))
    ))).astype(int)
    return round_to_nearest_power_ten(ticks, base)


def draw_axes(ax, x_ticks, x_ticks_relative, y_ticks, max_marker_count,
              min_aa_pos, max_aa_pos,
              protein_track_x_min: float, protein_track_x_max: float,
              protein_track_y_min: float, protein_track_y_max: float,
              font_size: int, text_padding: float,
              axis_color: str, protein_track_color: str
              ):
    # draw the tracks
    draw_rectangle(
        ax,
        protein_track_x_min, protein_track_y_min, protein_track_x_max, protein_track_y_max,
        line_color=protein_track_color, fill_color=protein_track_color, line_width=2.0,
    )
    # x_axis
    x_axis_y = protein_track_y_min - 0.02
    x_axis_min_x, x_axis_max_x = protein_track_x_min, protein_track_x_max
    big_tick_length, small_tick_length = 0.015, 0.005
    draw_line(  # main line
        ax,
        x_axis_min_x, x_axis_y, x_axis_max_x, x_axis_y,
        line_color=axis_color, line_width=1.0,
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
        line_color=axis_color, line_width=1.0,
    )
    draw_line(  # minimum tick
        ax,
        x_axis_min_x, x_axis_y - big_tick_length, x_axis_min_x, x_axis_y,
        line_color=axis_color, line_width=1.0,
    )
    for x_tick_relative, x_tick_absolute in zip(x_ticks_relative, x_ticks):
        draw_line(
            ax,
            x_tick_relative, x_axis_y - small_tick_length, x_tick_relative, x_axis_y,
            line_color=axis_color, line_width=1.0,
        )
        draw_string(
            ax,
            str(x_tick_absolute), x_tick_relative, x_axis_y - small_tick_length - text_padding,
            fontsize=font_size, ha='center', va='top',
        )

        # y_axis
        y_axis_x = protein_track_x_min - 0.02
        y_axis_min_y = protein_track_y_max + 0.01
        _, y_axis_max_y = marker_dim(max_marker_count, protein_track_y_max)
        y_ticks_relative = translate_to_ax_coordinates(
            y_ticks,
            min_absolute=0, max_absolute=max_marker_count, min_relative=y_axis_min_y, max_relative=y_axis_max_y,
        )
        draw_line(
            ax,
            y_axis_x, y_axis_min_y, y_axis_x, y_axis_max_y,
            line_color=axis_color, line_width=1.0,
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
            line_color=axis_color, line_width=1.0,
        )
        draw_line(  # max tick
            ax,
            y_axis_x - small_tick_length, y_axis_max_y, y_axis_x, y_axis_max_y,
            line_color=axis_color, line_width=1.0,
        )
        for y_tick_relative, y_tick_absolute in zip(y_ticks_relative, y_ticks):
            draw_line(
                ax,
                y_axis_x - small_tick_length, y_tick_relative, y_axis_x, y_tick_relative,
                line_color=axis_color, line_width=1.0,
            )
            draw_string(
                ax,
                str(y_tick_absolute), y_axis_x - small_tick_length - text_padding, y_tick_relative,
                fontsize=font_size, ha='right', va='center',
            )


def translate_to_ax_coordinates(
        absolute: np.ndarray,
        min_absolute, max_absolute,
        min_relative, max_relative,
        clip: bool = False,
) -> np.ndarray:
    if clip:
        # Put the coordinate of an item located at or after the stop codon to the location of the last AA
        absolute = np.minimum(absolute, max_absolute)
        # Put the coordinate of an item located at or before the start codon to the location of the first AA
        absolute = np.maximum(absolute, min_absolute)
    shifted_to_0_1 = ((absolute - min_absolute) / (max_absolute - min_absolute))
    relative_scale = (max_relative - min_relative)
    return shifted_to_0_1 * relative_scale + min_relative


def draw_marker(
        ax: plt.Axes,
        x_start, x_end, min_y, max_y, circle_radius, fill_color, stem_color='#a9a9a9',
):
    """
    Draw a lollipop representing a variant and the number of counts for the variant effect type
    currently putting marker in the middle of start and end, can change this later
    """
    x = (x_start + x_end) / 2
    draw_line(ax, x, min_y, x, max_y - circle_radius, line_color=stem_color, line_width=0.5)
    draw_circle(ax, x, max_y, circle_radius, line_color=stem_color, fill_color=fill_color, line_width=0.5)


def draw_variants(ax: plt.Axes, variant_locations_relative, variant_effect_colors, marker_counts, protein_track_y_max):
    for marker, marker_color in zip(variant_locations_relative, variant_effect_colors):
        marker_count = marker_counts[np.where(variant_locations_relative == marker)[0][0]]
        cur_radius, cur_length = marker_dim(marker_count, protein_track_y_max)
        x_start, x_end = marker, marker  # WAS  marker[0], marker[1]
        draw_marker(ax, x_start, x_end, protein_track_y_max, cur_length, cur_radius, marker_color)


def generate_variant_markers(pvis, marker_colors):
    variant_locations = pvis.variant_locations
    # The following has the variant positions (in variant_locations_counted_absolute)
    # and the number of occurrences of each position (in marker counts)
    variant_locations_counted_absolute = pvis.variant_locations_counted_absolute
    marker_counts = pvis.marker_counts
    variant_effect_colors = []
    for vl in variant_locations_counted_absolute:
        i = np.where(variant_locations == vl)[0][0]  # find index of unique variant loc in all locs to find effect
        effect = pvis.variant_effects[i]
        variant_effect_colors.append(marker_colors[effect])

    return marker_counts, variant_locations_counted_absolute, variant_effect_colors


def assign_colors(names, available_colors):
    num_colors = 0
    _colors = dict()
    seen = set()
    for n in names:
        if n in seen:
            continue
        seen.add(n)
        color = available_colors[num_colors]
        _colors[n] = color
        num_colors += 1

    return _colors


def draw_legends(ax: plt.Axes, feature_handler, pvis,
                 color_box_x_dim, color_box_y_dim, color_circle_radius, row_spacing,
                 legend1_min_x, legend1_max_y,
                 legend2_min_x, legend2_max_y,
                 variant_effect_colors, marker_colors,
                 labeling_method,
                 ):
    # draw legend 1 for protein features
    n_unique_features = len(feature_handler.cleaned_unique_feature_names)
    if labeling_method == 'abbreviate':
        color_box_x_dim *= 3.5
    legend1_width = 0.2 + color_box_x_dim
    legend1_min_y = legend1_max_y - (n_unique_features + 1) * row_spacing - n_unique_features * color_box_y_dim
    legend1_max_x = legend1_min_x + legend1_width

    # legend box
    draw_rectangle(ax, legend1_min_x, legend1_min_y, legend1_max_x, legend1_max_y, 'black')
    for i, feature_name in enumerate(feature_handler.cleaned_unique_feature_names):
        # colored box
        color_box_min_x = legend1_min_x + row_spacing
        color_box_max_x = color_box_min_x + color_box_x_dim
        color_box_max_y = legend1_max_y - (i + 1) * row_spacing - i * color_box_y_dim
        color_box_min_y = color_box_max_y - color_box_y_dim
        draw_rectangle(
            ax, color_box_min_x, color_box_min_y, color_box_max_x, color_box_max_y,
            line_color='black', fill_color=feature_handler.colors[feature_name],
        )
        # label in colored box (same as on x axis in the boxes)
        draw_string(
            ax, feature_handler.labels[feature_name], color_box_min_x + 0.002, color_box_min_y + 0.005,
            fontsize=10, ha="left", va="center", color='black',
        )
        # full feature name
        draw_string(
            ax, feature_name,
            color_box_max_x + 0.005, color_box_min_y + 0.005,
            fontsize=12, ha="left", va="center", color='black',
        )

    # draw legend 2 for variant effects
    unique_variant_effects = list(set(pvis.variant_effects))
    n_unique_effects = len(set(variant_effect_colors))
    legend2_min_y = legend2_max_y - (
            n_unique_effects + 1) * row_spacing - n_unique_effects * 2 * color_circle_radius
    legend2_max_x = legend2_min_x + 0.2
    draw_rectangle(ax, legend2_min_x, legend2_min_y, legend2_max_x, legend2_max_y, 'black')
    for i, variant_effect in enumerate(unique_variant_effects):
        colored_circle_x = legend2_min_x + row_spacing + color_circle_radius
        colored_circle_y = legend2_max_y - (i + 1) * row_spacing - i * 2 * color_circle_radius - color_circle_radius
        draw_circle(
            ax, colored_circle_x, colored_circle_y, color_circle_radius,
            line_color='black', fill_color=marker_colors[variant_effect],
        )
        draw_string(
            ax, variant_effect,
            colored_circle_x + 2 * color_circle_radius, colored_circle_y,
            fontsize=12, ha="left", va="center", color='black', )

    return legend1_width
