import typing
import warnings

from gpsea.model import (
    Cohort,
    ProteinMetadata,
    TranscriptAnnotation,
    TranscriptCoordinates,
    Variant,
    VariantEffect,
)
import numpy as np

from gpsea.model.genome._genome import Region


class ProteinVisualizable:

    def __init__(
        self,
        tx_coordinates: typing.Optional[TranscriptCoordinates],
        protein_meta: ProteinMetadata,
        cohort: Cohort,
    ):
        # TODO[v1.0.0] - Remove from the public API!
        warnings.warn(
            "ProteinVisualizable was deprecated and will be remove in 1.0.0",
            DeprecationWarning,
            stacklevel=2,
        )
        self._cohort = cohort
        self._protein_meta = protein_meta

        # store the annotations for the target transcript
        transcript_annotations = ProteinVisualizable._get_tx_anns(
            cohort.all_variants(), protein_meta.protein_id,
        )
        
        variant_regions_on_protein: typing.List[Region] = list()
        self._variant_effect = list()
        for tx_ann in transcript_annotations:
            if tx_ann is None:
                continue
            variant_effects = tx_ann.variant_effects
            if len(variant_effects) == 0:
                continue
            prot_eff_loc = tx_ann.protein_effect_location
            if prot_eff_loc is not None:
                variant_regions_on_protein.append(prot_eff_loc)
                self._variant_effect.append(variant_effects[0])

        self._protein_feature_names = list()
        self._protein_feature_types = list()
        self._protein_feature_starts = list()
        self._protein_feature_ends = list()
        for feature in protein_meta.protein_features:
            self._protein_feature_names.append(feature.info.name)
            self._protein_feature_types.append(feature.feature_type.lower())
            self._protein_feature_starts.append(feature.info.start)
            self._protein_feature_ends.append(feature.info.end)

        self._variant_locations = np.array(
            [item.start for item in variant_regions_on_protein]
        )

        # variant_locations = (variant_locations * 3) - 2 + min_exon_limit  # to convert from codons to bases
        # variant_effects = np.array([(ann.variant_effects[0]) for ann in tx_anns])
        # count marker occurrences and remove duplicates
        self._variant_locations_counted_absolute, self._marker_counts = np.unique(
            self._variant_locations,
            axis=0,
            return_counts=True,
        )

        if protein_meta.protein_length > 0:
            self._protein_length = protein_meta.protein_length
        else:
            raise ValueError(
                f"Unable to get protein_length for {protein_meta.label}. "
                "Consider deleting cache and trying again. Also check accession numbers"
            )

    @staticmethod
    def _get_tx_anns(
        variants: typing.Iterable[Variant],
        protein_id: str,
    ) -> typing.Sequence[typing.Optional[TranscriptAnnotation]]:
        """
        By default, the API returns transcript annotations for many transcripts.
        We would like to store the annotations only for our protein of interest (protein_id)
        """
        tx_anns = []
        for v in variants:
            tx_ann = None
            for ann in v.tx_annotations:
                if ann.protein_id is not None and ann.protein_id == protein_id:
                    tx_ann = ann
                    break
            tx_anns.append(tx_ann)

        return tx_anns

    @property
    def cohort(self) -> Cohort:
        return self._cohort

    @property
    def protein_id(self) -> str:
        return self._protein_meta.protein_id

    @property
    def protein_metadata(self) -> ProteinMetadata:
        return self._protein_meta

    @property
    def protein_feature_starts(self) -> typing.Sequence[int]:
        return self._protein_feature_starts

    @property
    def protein_feature_ends(self) -> typing.Sequence[int]:
        return self._protein_feature_ends

    @property
    def protein_feature_types(self) -> typing.Sequence[str]:
        return self._protein_feature_types

    @property
    def protein_length(self) -> int:
        return self._protein_length
        # TODO: ignoring the code below because we check the length in __init__ method.
        #  Remove the commented code if there are no issues.
        # """
        # We try to parse the protein length from the UniProt API. If it worked, then self._protein_length is greater than zero.
        # If it did not work, then we take the maximum value of the protein features and variants and add 30 amino acids to it so that
        # the display will show something reasonable. We also print an error.
        # """
        # if self._protein_length > 0:
        #     return self._protein_length
        # else:
        #     print("There was some problem parsing the protein length; here we estimate the length based on features and variants")
        #     max_feat = max(self._protein_feature_ends)
        #     max_var = max(self._variant_locations)
        #     return  max(max_feat, max_var) + 30

    @property
    def protein_feature_names(self) -> typing.Sequence[str]:
        return self._protein_feature_names

    @property
    def variant_effects(self) -> typing.Sequence[VariantEffect]:
        return self._variant_effect

    @property
    def variant_locations_counted_absolute(self) -> np.ndarray:
        return self._variant_locations_counted_absolute

    @property
    def marker_counts(self) -> np.ndarray:
        # `self._marker_counts` should be a 1D array with non-negative ints (counts).
        return self._marker_counts

    @property
    def variant_locations(self) -> np.ndarray:
        """
        Get an array with 0-based start coordinates of aminoacids that overlap with the variant regions.

        Returns
          a 1D array with ints containing a coordinate for each variant
        """
        return self._variant_locations
