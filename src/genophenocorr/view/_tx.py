import typing

from genophenocorr.model import Variant
from genophenocorr.model.genome import TranscriptCoordinates


class VariantTranscriptArtist:
    """
    `VariantTranscriptArtist` creates a graphic to show distributions of variants across a provided transcript.
    """

    def draw_variants(self, tx: TranscriptCoordinates, variants: typing.Sequence[Variant]):
        # TODO - implement
        pass
