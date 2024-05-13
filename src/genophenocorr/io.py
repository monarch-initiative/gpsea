import typing
from json import JSONDecoder, JSONEncoder

import hpotk

from genophenocorr.model import *
from genophenocorr.model.genome import *


class GenophenocorrJSONEncoder(JSONEncoder):
    """
    `GenophenocorrJSONEncoder` encodes genophenocorr's types into a JSON message.

    The encoder is supposed to be used along with Python's `json` module via the `cls` parameter of :func:`json.dump`
    or :func:`json.dumps`.
    """

    def default(self, o):
        if isinstance(o, Variant):
            return {
                'variant_coordinates': o.variant_coordinates,
                'tx_annotations': o.tx_annotations,
                'genotypes': o.genotypes,
            }
        elif isinstance(o, VariantCoordinates):
            return {
                'region': o.region,
                'ref': o.ref,
                'alt': o.alt,
                'change_length': o.change_length,
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
        elif isinstance(o, (Genotype, VariantEffect, Strand)):
            return o.name
        elif isinstance(o, Phenotype):
            return {
                'term_id': o.identifier.value,
                'name': o.name,
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
                'phenotypes': o.phenotypes,
                'diseases': o.diseases,
                'variants': o.variants,
            }
        elif isinstance(o, Cohort):
            return {
                'members': sorted(o.all_patients, key=lambda patient: patient.labels),
                'excluded_patient_count': o.get_excluded_count(),
            }
        else:
            return super().default(o)


_VARIANT_FIELDS = ('variant_coordinates', 'tx_annotations', 'genotypes')
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
_PHENOTYPE_FIELDS = ('term_id', 'name', 'is_present')
_DISEASE_FIELDS = ('term_id', 'name', 'is_observed')
_PATIENT_FIELDS = ('labels', 'phenotypes', 'diseases', 'variants')
_COHORT_FIELDS = ('members', 'excluded_patient_count')


class GenophenocorrJSONDecoder(JSONDecoder):
    """
    `GenophenocorrJSONDecoder` decodes genophenocorr's types from a JSON message.

    The decoder is supposed to be used along with Python's `json` module via the `cls` parameter of :func:`json.load`
    or :func:`json.loads`.
    """

    def __init__(
            self,
            *args: typing.Any, **kwargs: typing.Any,
    ):
        super().__init__(
            *args,
            object_hook=GenophenocorrJSONDecoder.object_hook,
            **{k: v for k, v in kwargs.items() if k != 'object_hook'},
        )

    @staticmethod
    def _has_all_fields(
            obj: typing.Dict[str, typing.Any],
            query: typing.Iterable[str],
    ) -> bool:
        return all(key in obj for key in query)

    @staticmethod
    def object_hook(obj: typing.Dict[typing.Any, typing.Any]) -> typing.Any:
        if GenophenocorrJSONDecoder._has_all_fields(obj, _VARIANT_FIELDS):
            return Variant(
                var_coordinates=obj['variant_coordinates'],
                tx_annotations=obj['tx_annotations'],
                genotypes=obj['genotypes'],
            )
        elif GenophenocorrJSONDecoder._has_all_fields(obj, _VARIANT_COORDINATES_FIELDS):
            return VariantCoordinates(
                region=obj['region'],
                ref=obj['ref'],
                alt=obj['alt'],
                change_length=obj['change_length'],
            )
        elif GenophenocorrJSONDecoder._has_all_fields(obj, _REGION_FIELDS):
            if GenophenocorrJSONDecoder._has_all_fields(obj, _GENOMIC_REGION_FIELDS):
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
        elif GenophenocorrJSONDecoder._has_all_fields(obj, _CONTIG_FIELDS):
            return Contig(
                name=obj['name'],
                gb_acc=obj['genbank_acc'],
                refseq_name=obj['refseq_name'],
                ucsc_name=obj['ucsc_name'],
                length=obj['length'],
            )
        elif GenophenocorrJSONDecoder._has_all_fields(obj, _SAMPLE_LABELS_FIELDS):
            return SampleLabels(
                label=obj['label'],
                meta_label=obj['meta_label'],
            )
        elif GenophenocorrJSONDecoder._has_all_fields(obj, _GENOTYPES_FIELDS):
            return Genotypes(
                samples=obj['samples'],
                genotypes=(Genotype[gt] for gt in obj['genotypes']),
            )
        elif GenophenocorrJSONDecoder._has_all_fields(obj, _TX_ANNOTATION_FIELDS):
            return TranscriptAnnotation(
                gene_id=obj['gene_id'],
                tx_id=obj['transcript_id'],
                hgvs_cdna=obj['hgvs_cdna'],
                is_preferred=obj['is_preferred'],
                variant_effects=(VariantEffect[ve] for ve in obj['variant_effects']),
                affected_exons=obj['overlapping_exons'],
                protein_id=obj['protein_id'],
                protein_effect_coordinates=obj['protein_effect_location'],
            )
        elif GenophenocorrJSONDecoder._has_all_fields(obj, _PHENOTYPE_FIELDS):
            return Phenotype(
                term_id=hpotk.TermId.from_curie(obj['term_id']),
                name=obj['name'],
                is_observed=obj['is_present'],
            )
        elif GenophenocorrJSONDecoder._has_all_fields(obj, _DISEASE_FIELDS):
            return Disease(
                term_id=hpotk.TermId.from_curie(obj['term_id']),
                name=obj['name'],
                is_observed=obj['is_observed'],
            )
        elif GenophenocorrJSONDecoder._has_all_fields(obj, _PATIENT_FIELDS):
            return Patient(
                labels=obj['labels'],
                phenotypes=obj['phenotypes'],
                diseases=obj['diseases'],
                variants=obj['variants'],
            )
        elif GenophenocorrJSONDecoder._has_all_fields(obj, _COHORT_FIELDS):
            return Cohort(
                members=obj['members'],
                excluded_member_count=obj['excluded_patient_count'],
            )
        else:
            return obj
