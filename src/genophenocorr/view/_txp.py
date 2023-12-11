import typing

from genophenocorr.model import Variant, TranscriptCoordinates, ProteinMetadata


class VariantTranscriptProteinArtist:
    """
    `VariantTranscriptProteinArtist` creates a graphic to show distributions of variants across a provided transcript
    and its protein product.
    """

    def draw_variants(self, variants: typing.Iterable[Variant],
                      tx: TranscriptCoordinates,
                      protein: ProteinMetadata):
        # TODO - Peter implement
        pass
