from genophenocorr.view import VariantTranscriptProteinArtist


from .genome_data import *


def test_draw_variants_on_transcript_pos(toy_variants: typing.Sequence[Variant],
                                         toy_tx_positive: TranscriptCoordinates,
                                         toy_protein_positive: ProteinMetadata):
    artist = VariantTranscriptProteinArtist()
    output = artist.draw_variants(toy_variants, toy_tx_positive, toy_protein_positive)
    print(output)


def test_draw_variants_on_transcript_neg(toy_variants: typing.Sequence[Variant],
                                         toy_tx_negative: TranscriptCoordinates,
                                         toy_protein_negative: ProteinMetadata):
    artist = VariantTranscriptProteinArtist()
    output = artist.draw_variants(toy_variants, toy_tx_negative, toy_protein_negative)
    print(output)
