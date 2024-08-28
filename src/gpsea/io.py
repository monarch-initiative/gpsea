import typing

import hpotk

from json import JSONDecoder, JSONEncoder

from gpsea.model import *
from gpsea.model.genome import *


class GpseaJSONEncoder(JSONEncoder):
    """
    `GpseaJSONEncoder` encodes gpsea's types into a JSON message.

    The encoder is supposed to be used along with Python's `json` module via the `cls` parameter of :func:`json.dump`
    or :func:`json.dumps`.
    """

    def default(self, o):
        if isinstance(o, Variant):
            return {
                'variant_info': o.variant_info,
                'tx_annotations': o.tx_annotations,
                'genotypes': o.genotypes,
            }
        elif isinstance(o, VariantInfo):
            return {
                'variant_coordinates': o.variant_coordinates,
                'sv_info': o.sv_info,
            }
        elif isinstance(o, VariantCoordinates):
            return {
                'region': o.region,
                'ref': o.ref,
                'alt': o.alt,
                'change_length': o.change_length,
            }
        elif isinstance(o, ImpreciseSvInfo):
            return {
                'structural_type': o.structural_type.value,
                'variant_class': o.variant_class,
                'gene_id': o.gene_id,
                'gene_symbol': o.gene_symbol
            }
        elif isinstance(o, Region):
            val = {
                'start': o.start,
                'end': o.end,
            }
            if isinstance(o, GenomicRegion):
                val['contig'] = o.contig
                val['strand'] = o.strand

            return val
        elif isinstance(o, Contig):
            return {
                'name': o.name,
                'genbank_acc': o.genbank_acc,
                'refseq_name': o.refseq_name,
                'ucsc_name': o.ucsc_name,
                'length': len(o),
            }
        elif isinstance(o, TranscriptAnnotation):
            return {
                'gene_id': o.gene_id,
                'transcript_id': o.transcript_id,
                'hgvs_cdna': o.hgvs_cdna,
                'is_preferred': o.is_preferred,
                'variant_effects': o.variant_effects,
                'overlapping_exons': o.overlapping_exons,
                'protein_id': o.protein_id,
                'hgvsp': o.hgvsp,
                'protein_effect_location': o.protein_effect_location,
            }
        elif isinstance(o, Genotypes):
            samples = []
            genotypes = []
            for s, g in o:
                samples.append(s)
                genotypes.append(g)
            return {
                'samples': samples,
                'genotypes': genotypes,
            }
        elif isinstance(o, SampleLabels):
            return {
                'label': o.label,
                'meta_label': o.meta_label,
            }
        elif isinstance(o, (Sex, Genotype, VariantEffect, Strand, VariantClass)):
            # enums
            return o.name
        elif isinstance(o, Phenotype):
            return {
                'term_id': o.identifier.value,
                'is_present': o.is_present,
            }
        elif isinstance(o, Disease):
            return {
                'term_id': o.identifier.value,
                'name': o.name,
                'is_observed': o.is_present,
            }
        elif isinstance(o, Patient):
            return {
                'labels': o.labels,
                'sex': o.sex,
                'phenotypes': o.phenotypes,
                'diseases': o.diseases,
                'variants': o.variants,
            }
        elif isinstance(o, Cohort):
            return {
                'members': sorted(o.all_patients, key=lambda patient: patient.labels),
                'excluded_patient_count': o.get_excluded_count(),
            }
        elif isinstance(o, TranscriptCoordinates):
            return {
                'identifier': o.identifier,
                'region': o.region,
                'exons': o.exons,
                'cds_start': o.cds_start,
                'cds_end': o.cds_end,
            }
        elif isinstance(o, ProteinMetadata):
            return {
                'protein_id': o.protein_id,
                'label': o.label,
                'protein_features': o.protein_features,
                'protein_length': o.protein_length,
            }
        elif isinstance(o, ProteinFeature):
            return {
                'info': o.info,
                'feature_type': o.feature_type.name,
            }
        elif isinstance(o, FeatureInfo):
            return {
                'name': o.name,
                'region': o.region,
            }
        else:
            return super().default(o)


_VARIANT_FIELDS = ('variant_info', 'tx_annotations', 'genotypes')
_VARIANT_INFO_FIELDS = ('variant_coordinates', 'sv_info')
_IMPRECISE_SV_INFO_FIELDS = ('structural_type', 'variant_class', 'gene_id', 'gene_symbol')
_VARIANT_COORDINATES_FIELDS = ('region', 'ref', 'alt', 'change_length')
_REGION_FIELDS = ('start', 'end')
_GENOMIC_REGION_FIELDS = ('contig', 'start', 'end', 'strand')
_CONTIG_FIELDS = ('name', 'genbank_acc', 'refseq_name', 'ucsc_name', 'length')
_SAMPLE_LABELS_FIELDS = ('label', 'meta_label')
_GENOTYPES_FIELDS = ('samples', 'genotypes')
_TX_ANNOTATION_FIELDS = (
    'gene_id', 'transcript_id', 'hgvs_cdna', 'is_preferred', 'variant_effects',
    'overlapping_exons', 'protein_id', 'protein_effect_location',
)
_TX_COORDINATES = ('identifier', 'region', 'exons', 'cds_start', 'cds_end')
_PROTEIN_METADATA = ('protein_id', 'label', 'protein_features', 'protein_length')
_PROTEIN_FEATURE = ('info', 'feature_type')
_FEATURE_INFO = ('name', 'region')
_PHENOTYPE_FIELDS = ('term_id', 'is_present')
_DISEASE_FIELDS = ('term_id', 'name', 'is_observed')
_PATIENT_FIELDS = ('labels', 'sex', 'phenotypes', 'diseases', 'variants')
_COHORT_FIELDS = ('members', 'excluded_patient_count')


class GpseaJSONDecoder(JSONDecoder):
    """
    `GpseaJSONDecoder` decodes gpsea's types from a JSON message.

    The decoder is supposed to be used along with Python's `json` module via the `cls` parameter of :func:`json.load`
    or :func:`json.loads`.
    """

    def __init__(
            self,
            *args: typing.Any, **kwargs: typing.Any,
    ):
        super().__init__(
            *args,
            object_hook=GpseaJSONDecoder.object_hook,
            **{k: v for k, v in kwargs.items() if k != 'object_hook'},
        )

    @staticmethod
    def _has_all_fields(
            obj: typing.Dict[str, typing.Any],
            query: typing.Iterable[str],
    ) -> bool:
        return all(key in obj for key in query)
    
    @staticmethod
    def _has_any_field(
            obj: typing.Dict[str, typing.Any],
            query: typing.Iterable[str],
    ) -> bool:
        return any(key in obj for key in query)

    @staticmethod
    def object_hook(obj: typing.Dict[typing.Any, typing.Any]) -> typing.Any:
        if GpseaJSONDecoder._has_all_fields(obj, _VARIANT_FIELDS):
            return Variant(
                variant_info=obj['variant_info'],
                tx_annotations=obj['tx_annotations'],
                genotypes=obj['genotypes'],
            )
        elif GpseaJSONDecoder._has_any_field(obj, _VARIANT_INFO_FIELDS):
            return VariantInfo(
                variant_coordinates=obj['variant_coordinates'] if 'variant_coordinates' in obj else None,
                sv_info=obj['sv_info'] if 'sv_info' in obj else None,
            )
        elif GpseaJSONDecoder._has_all_fields(obj, _VARIANT_COORDINATES_FIELDS):
            return VariantCoordinates(
                region=obj['region'],
                ref=obj['ref'],
                alt=obj['alt'],
                change_length=obj['change_length'],
            )
        elif GpseaJSONDecoder._has_all_fields(obj, _IMPRECISE_SV_INFO_FIELDS):
            return ImpreciseSvInfo(
                structural_type=hpotk.TermId.from_curie(obj['structural_type']),
                variant_class=VariantClass[obj['variant_class']],
                gene_id=obj['gene_id'],
                gene_symbol=obj['gene_symbol'],
            )
        elif GpseaJSONDecoder._has_all_fields(obj, _REGION_FIELDS):
            if GpseaJSONDecoder._has_all_fields(obj, _GENOMIC_REGION_FIELDS):
                return GenomicRegion(
                    contig=obj['contig'],
                    start=obj['start'],
                    end=obj['end'],
                    strand=Strand[obj['strand']],
                )
            else:
                return Region(
                    start=obj['start'],
                    end=obj['end'],
                )
        elif GpseaJSONDecoder._has_all_fields(obj, _CONTIG_FIELDS):
            return Contig(
                name=obj['name'],
                gb_acc=obj['genbank_acc'],
                refseq_name=obj['refseq_name'],
                ucsc_name=obj['ucsc_name'],
                length=obj['length'],
            )
        elif GpseaJSONDecoder._has_all_fields(obj, _SAMPLE_LABELS_FIELDS):
            return SampleLabels(
                label=obj['label'],
                meta_label=obj['meta_label'],
            )
        elif GpseaJSONDecoder._has_all_fields(obj, _GENOTYPES_FIELDS):
            return Genotypes(
                samples=obj['samples'],
                genotypes=(Genotype[gt] for gt in obj['genotypes']),
            )
        elif GpseaJSONDecoder._has_all_fields(obj, _TX_ANNOTATION_FIELDS):
            return TranscriptAnnotation(
                gene_id=obj['gene_id'],
                tx_id=obj['transcript_id'],
                hgvs_cdna=obj['hgvs_cdna'],
                is_preferred=obj['is_preferred'],
                variant_effects=(VariantEffect[ve] for ve in obj['variant_effects']),
                affected_exons=obj['overlapping_exons'],
                protein_id=obj['protein_id'],
                hgvsp=obj['hgvsp'],
                protein_effect_coordinates=obj['protein_effect_location'],
            )
        elif GpseaJSONDecoder._has_all_fields(obj, _PHENOTYPE_FIELDS):
            return Phenotype(
                term_id=hpotk.TermId.from_curie(obj['term_id']),
                is_observed=obj['is_present'],
            )
        elif GpseaJSONDecoder._has_all_fields(obj, _DISEASE_FIELDS):
            return Disease(
                term_id=hpotk.TermId.from_curie(obj['term_id']),
                name=obj['name'],
                is_observed=obj['is_observed'],
            )
        elif GpseaJSONDecoder._has_all_fields(obj, _PATIENT_FIELDS):
            return Patient(
                labels=obj['labels'],
                sex=Sex[obj['sex']],
                phenotypes=obj['phenotypes'],
                diseases=obj['diseases'],
                variants=obj['variants'],
            )
        elif GpseaJSONDecoder._has_all_fields(obj, _TX_COORDINATES):
            return TranscriptCoordinates(
                identifier=obj['identifier'],
                region=obj['region'],
                exons=obj['exons'],
                cds_start=obj['cds_start'],
                cds_end=obj['cds_end'],
            )
        elif GpseaJSONDecoder._has_all_fields(obj, _PROTEIN_METADATA):
            return ProteinMetadata(
                protein_id=obj['protein_id'],
                label=obj['label'],
                protein_features=obj['protein_features'],
                protein_length=obj['protein_length'],
            )
        elif GpseaJSONDecoder._has_all_fields(obj, _PROTEIN_FEATURE):
            return ProteinFeature.create(
                info=obj['info'], feature_type=FeatureType[obj['feature_type']],
            )
        elif GpseaJSONDecoder._has_all_fields(obj, _FEATURE_INFO):
            return FeatureInfo(
                name=obj['name'], region=obj['region'],
            )
        elif GpseaJSONDecoder._has_all_fields(obj, _COHORT_FIELDS):
            return Cohort(
                members=obj['members'],
                excluded_member_count=obj['excluded_patient_count'],
            )
        else:
            return obj
