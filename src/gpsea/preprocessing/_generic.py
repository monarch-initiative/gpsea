import typing

from gpsea.model import ImpreciseSvInfo, TranscriptAnnotation, VariantEffect, VariantClass

from ._api import ImpreciseSvFunctionalAnnotator, GeneCoordinateService


class DefaultImpreciseSvFunctionalAnnotator(ImpreciseSvFunctionalAnnotator):

    def __init__(
        self,
        gene_coordinate_service: GeneCoordinateService,
    ):
        self._coordinate_service = gene_coordinate_service

    def annotate(self, item: ImpreciseSvInfo) -> typing.Sequence[TranscriptAnnotation]:
        gene_id = item.gene_id
        tx_coordinates = self._coordinate_service.fetch_for_gene(gene_id)
        
        tx_annotations = []
        for txc in tx_coordinates:
            is_preferred = False if txc.is_preferred is None else txc.is_preferred
            variant_effects = self._map_to_variant_effects(item.variant_class)
            
            # We assume all exons are affected
            affected_exons = range(len(txc.exons))  
            annotation = TranscriptAnnotation(
                gene_id=gene_id,
                tx_id=txc.identifier,
                hgvs_cdna=None,  # We can't provide this here at this time.
                is_preferred=is_preferred,
                variant_effects=variant_effects,
                affected_exons=affected_exons,
                protein_id=None,
                hgvsp=None,
                protein_effect_coordinates=None,
            )
            tx_annotations.append(annotation)
        
        return tuple(tx_annotations)
    
    def _map_to_variant_effects(
        self, 
        variant_class: str,
    ) -> typing.Sequence[VariantEffect]:
        if variant_class == VariantClass.DEL:
            return (VariantEffect.TRANSCRIPT_ABLATION,)
        elif variant_class == VariantClass.DUP:
            return (VariantEffect.TRANSCRIPT_AMPLIFICATION,)
        else:
            # This mapping is most likely incomplete.
            # Please open a ticket if support 
            # for a variant class should be added.
            raise ValueError(f'Unsupported variant class {variant_class}')
