from genophenocorr.model import Cohort, ProteinMetadata, TranscriptCoordinates, VariantEffect, FeatureType
import numpy as np



class ProteinVisualizable:
    def __init__(self, tx_coordinates: TranscriptCoordinates, protein_meta: ProteinMetadata, cohort: Cohort) -> None:
        self._tx_id = tx_coordinates.identifier
        self._protein_id = protein_meta.protein_id
        self._protein_features = list()
        for feature in protein_meta.protein_features:
            self._protein_features.append(feature)
        variants = cohort.all_variants()
         # store the annotations for the correct transcript
        transcript_annotations = ProteinVisualizable._get_tx_anns(variants, self._tx_id)
        self._variant_pos = list()
        self._variant_effect = list()
        for tannot in transcript_annotations:
            variant_effects = tannot.variant_effects
            if variant_effects is None or len(variant_effects) == 0:
                continue
            prot_eff_loc = tannot.protein_effect_location
            if prot_eff_loc is not None:
                self._variant_pos.append([prot_eff_loc.start, prot_eff_loc.end])
                self._variant_effect.append(variant_effects[0])
        self._protein_feature_names = list()
        self._protein_feature_starts = list()
        self._protein_feature_ends = list()
        for feature in protein_meta.protein_features:
            self._protein_feature_names.append(feature.info.name)
            self._protein_feature_starts.append(feature.info.start)
            self._protein_feature_ends.append(feature.info.end)

        variant_locations = list(item[0] for item in self._variant_pos)
        self._variant_locations = np.array(variant_locations)

        #variant_locations = (variant_locations * 3) - 2 + min_exon_limit  # to convert from codons to bases
        #variant_effects = np.array([(ann.variant_effects[0]) for ann in tx_anns])
        # count marker occurrences and remove duplicates
        variant_locations_counted_absolute, marker_counts = np.unique(variant_locations, axis=0, return_counts=True)
        self._variant_locations_counted_absolute = variant_locations_counted_absolute
        self._marker_counts = marker_counts
        if protein_meta.protein_length > 0:
            self._protein_length = protein_meta.protein_length
        else:
            raise ValueError(f"Unable to get protein_length for {protein_meta.label}. Consider deleting cache and trying again. Also check accession numbers")

    @staticmethod
    def _get_protein_meta(protein_meta: ProteinMetadata):
        protein_id = protein_meta.protein_id
        for pm in protein_meta:
            if pm._id == protein_id:
                return pm
    
    @staticmethod
    def _get_tx_anns(variants, tx_id):
        """
        By default, the API returns transcript annotations for many transcripts. 
        We would like to store the annotations only for our transcript of interest (tx_id)
        """
        tx_anns = []
        for i, v in enumerate(variants):
            tx_ann = None
            for ann in v.tx_annotations:
                if ann.transcript_id == tx_id:
                    tx_ann = ann
                    break
            if tx_ann is None:
                raise ValueError(f'The transcript annotation for {tx_id} was not found!')
            else:
                tx_anns.append(tx_ann)

        return tx_anns
    
    @property
    def transcript_id(self):
        return self._tx_id
    
    @property
    def protein_id(self):
        return self._protein_id
    
    @property
    def protein_feature_starts(self):
        return self._protein_feature_starts
    @property
    def protein_feature_ends(self):
        return self._protein_feature_ends
    
    @property
    def protein_length(self):
        """
        We try to parse the protein length from the UniProt API. If it worked, then self._protein_length is greater than zero.
        If it did not work, then we take the maximum value of the protein features and variants and add 30 amino acids to it so that
        the display will show something reasonable. We also print an error.
        """
        if self._protein_length > 0:
            return self._protein_length
        else:
            print("There was some problem parsing the protein length; here we estimate the length based on features and variants")
            max_feat = max(self._protein_feature_ends)
            max_var = max(self._variant_locations)
            return  max(max_feat, max_var) + 30
    
    @property
    def protein_feature_names(self):
        return self._protein_feature_names
    
    @property
    def variant_effects(self):
        return self._variant_effect
    
    @property
    def variant_locations_counted_absolute(self):
        return self._variant_locations_counted_absolute
    @property
    def marker_counts(self):
        return self._marker_counts
    
    @property
    def variant_locations(self):
        return self._variant_locations
   