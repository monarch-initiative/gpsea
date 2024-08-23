import typing
from collections import defaultdict
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from matplotlib.lines import Line2D
import typing
from gpsea.model import Variant, TranscriptCoordinates, ProteinMetadata


def get_interpolated_location_in_protein(aa_pos_in_prot:int,
                                        prot_len:int,
                                        plot_width:int=1200,
                                        margin_size:int=20) -> typing.Tuple[int,int]:
    used_plot_width = plot_width - 2*margin_size
    rel_pos_in_prot = aa_pos_in_prot / prot_len
    interpolated_x = margin_size + used_plot_width * rel_pos_in_prot
    return interpolated_x



class VariantTranscriptVisualizer:
    """
    `VariantTranscriptProteinArtist` creates a graphic to show distributions of variants across a provided transcript
    and its protein product.
    """

    # noinspection PyMethodMayBeStatic
    def draw_variants(self, variants: typing.Iterable[Variant],
                    tx: TranscriptCoordinates,
                    protein: ProteinMetadata):
        title = f"{protein.protein_id} ({protein.label})"
        fig, ax = plt.subplots(1, figsize=(10, 10))
        protein_domains = set()
        THRESHOLD = 2
        BOTTOM_MARGIN = 20
        amino_acid_len = tx.get_codon_count()
        # draw a box that is ten aax tall, where aax is the dimension of one amino acid
        prot_start = get_interpolated_location_in_protein(1, amino_acid_len)
        prot_end = get_interpolated_location_in_protein(amino_acid_len, amino_acid_len)
        box_height = 10/amino_acid_len
        prot_width = prot_end - prot_start + 1
        protein_height = prot_width/20
        #rect = Rectangle((prot_start, BOTTOM_MARGIN), prot_width, protein_height)
        rect = Rectangle((50, 30), prot_width, 40, edgecolor='black', facecolor='white', alpha=0.5)
        ax.add_patch(rect)
        for variant in variants:
            tx_annot_list = [txa for txa in variant.tx_annotations if txa.transcript_id == tx.identifier]
            if len(tx_annot_list) == 0:
                raise ValueError(f"Could not find variant transcript annotations for {tx.identifier}")
            if len(tx_annot_list) > 1:
                # Should never happen
                raise ValueError(f"Found mutiple variant transcript annotations for {tx.identifier}")
            tx_annot = tx_annot_list[0]
            hgvs = tx_annot.hgvsc_id.split(":")
            if len(hgvs) > 1:
                hgvs_cdna = hgvs[1]
            else:
                hgvs_cdna = hgvs
            variant_effects = tx_annot.variant_effects
            if len(variant_effects) > 1:
                var_effect = "MULTIPLE"
            elif len(variant_effects) == 0:
                var_effect = "UNKNOWN"
            else:
                var_effect = variant_effects[0].name
            for p in tx_annot.protein_affected:
                for f in p.domains():
                    protein_domains.add(f.info)
            if tx_annot.protein_effect_location is None:
                print(f"skipping {hgvs_cdna} because no protein effect location was found")
                continue
            # for plotting, we will just show the start position; note there is also end
            pe_loc = tx_annot.protein_effect_location.start
            prot_variant_location = get_interpolated_location_in_protein(pe_loc, amino_acid_len)
            x = [prot_variant_location, prot_variant_location]
            ybase = BOTTOM_MARGIN + protein_height
            ytop = ybase + 10
            y = [ybase, ytop]
            line = Line2D(x, y)
            ax.add_line(line)
            # TODO show one box per mutation or otherwise represent frequency
            color_list = ['red', 'blue', 'green', 'orange', 'brown', 'yellow', 'purple']
            feature_to_color_index_d = defaultdict(int) # default is red
            for feature in protein_domains:
                name = feature.name
                if name not in feature_to_color_index_d:
                    new_idx = len(feature_to_color_index_d) % len(color_list)
                    feature_to_color_index_d[name] = new_idx
                color_idx = feature_to_color_index_d.get(name)
                box_color = color_list[color_idx]
                start = feature.start
                end = feature.end
                #print(name, start, end, box_color)
                box_height = 10/amino_acid_len
                prot_width = prot_end - prot_start + 1
                protein_height = prot_width/20
                #rect = Rectangle((prot_start, BOTTOM_MARGIN), prot_width, protein_height)
                box_start = get_interpolated_location_in_protein(start,  amino_acid_len)
                box_end = get_interpolated_location_in_protein(end, amino_acid_len)
                box_width = box_end - box_start
                rect = Rectangle((box_start, 30), box_width, 40, edgecolor='black', facecolor=box_color, alpha=0.5)
                ax.add_patch(rect)
        ax.set_title(title, loc='left', fontsize='x-large')
        ax.set_xlabel('X Label', labelpad=12, fontsize=14.5)
        ax.set_ylabel('Y label', labelpad=12, fontsize=14.5)
        ax.set_xlim(-10, 1200)
        ax.set_ylim(0,120)
