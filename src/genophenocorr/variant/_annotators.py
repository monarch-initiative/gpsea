import logging
import os
import pickle
import re
import typing

import hpotk
import requests
# pyright: reportGeneralTypeIssues=false
from phenopackets import GenomicInterpretation

from ._models import FunctionalAnnotator, VariantCoordinateFinder
from ._variant_data import VariantCoordinates, TranscriptAnnotation, Variant


class PhenopacketVariantCoordinateFinder(VariantCoordinateFinder[GenomicInterpretation]):

    def __init__(self):
        self._logger = logging.getLogger(__name__)

    def find_coordinates(self, item: GenomicInterpretation) -> VariantCoordinates:
        # TODO(ielis&lnrekerle) - we need to write tests for this logic
        if not isinstance(item, GenomicInterpretation):
            raise ValueError(f"item must be a Phenopacket GenomicInterpretation but was type {type(item)}")
        chrom, ref, alt, genotype = None, None, None, None
        start, end = 0, 0
        variant_descriptor = item.variant_interpretation.variation_descriptor
        if len(variant_descriptor.vcf_record.chrom) == 0 and len(
                variant_descriptor.variation.copy_number.allele.sequence_location.sequence_id) != 0:
            ref = 'N'
            start = int(
                variant_descriptor.variation.copy_number.allele.sequence_location.sequence_interval.start_number.value)
            end = int(
                variant_descriptor.variation.copy_number.allele.sequence_location.sequence_interval.end_number.value)
            number = int(variant_descriptor.variation.copy_number.number.value)
            if number > 2:
                alt = '<DUP>'
            else:
                alt = '<DEL>'
            chrom = re.findall(r'NC_0000(\d{2}).\d\d',
                               variant_descriptor.variation.copy_number.allele.sequence_location.sequence_id)[0]
            if chrom.startswith('0'):
                chrom = str(int(chrom))
            elif chrom == '23':
                chrom = 'X'
            elif chrom == '24':
                chrom = 'Y'
        elif len(variant_descriptor.vcf_record.chrom) != 0 and len(
                variant_descriptor.variation.copy_number.allele.sequence_location.sequence_id) == 0:
            ref = variant_descriptor.vcf_record.ref
            alt = variant_descriptor.vcf_record.alt
            start = int(variant_descriptor.vcf_record.pos) - 1
            end = int(variant_descriptor.vcf_record.pos) + abs(len(alt) - len(ref))
            chrom = variant_descriptor.vcf_record.chrom[3:]
        genotype = variant_descriptor.allelic_state.label

        if any(field is None for field in (chrom, ref, alt, genotype)):
            raise ValueError(f'Cannot determine variant coordinate from genomic interpretation {item}')

        return VariantCoordinates(chrom, start, end, ref, alt, len(alt) - len(ref), genotype)


def verify_start_end_coordinates(vc: VariantCoordinates):
    """
    Converts the 0-based VariantCoordinates to ones that will be interpreted correctly by VEP
    
    Example - an insertion/duplication of G after the given G at coordinate 3:
    1 2 3 4 5
    A C G T A

    0-based: 2 3 G GG       1-based: 3 G GG         VEP: 4 3 - G 
    """

    chrom = vc.chrom
    end = vc.end
    start = vc.start + 1
    alt = vc.alt
    if vc.is_structural():
        alt = vc.alt[1:-1]
        # TODO: Verify <INS> are working correctly
    else:
        if len(vc.ref) == 0 or len(vc.alt) == 0:
            raise ValueError(f'Trimmed alleles are not yet supported!')
        else:
            if len(vc.ref) == 1 and len(vc.alt) != 1:
                # INS/DUP
                start = start + 1  # we must "trim"
                alt = vc.alt[1:]
            # 100 AC AGT
            # MNV

    return f'{chrom}:{start}-{end}/{alt}'


class VepFunctionalAnnotator(FunctionalAnnotator):

    def __init__(self):
        self._logging = logging.getLogger(__name__)
        self._url = 'https://rest.ensembl.org/vep/human/region/%s?LoF=1&canonical=1&domains=1&hgvs=1' \
                    '&mutfunc=1&numbers=1&protein=1&refseq=1&transcript_version=1&variant_class=1'

    def annotate(self, variant_coordinates: VariantCoordinates) -> Variant:
        variant = self._query_vep(variant_coordinates)
        variant_id = variant.get('id')
        variant_class = variant.get('variant_class')
        # TODO - which one is the MANE transcript?
        canon_tx = None
        transcript_list = []
        for trans in variant.get('transcript_consequences'):
            trans_id = trans.get('transcript_id')
            if not trans_id.startswith('NM'):
                self._logging.info(f'We will not be using this transcript: {trans_id}')
                continue
            if trans.get('canonical') == 1:
                canon_tx = trans_id
            hgvsc_id = trans.get('hgvsc')
            consequences = trans.get('consequence_terms')
            gene_name = trans.get('gene_symbol')
            protein_id = trans.get('protein_id')
            protein_effect_start = trans.get('protein_start')
            protein_effect_end = trans.get('protein_end')
            if protein_effect_start is None and protein_effect_end is not None:
                protein_effect_start = 1
            if protein_effect_end is not None:
                protein_effect_end = int(protein_effect_end)
            if protein_effect_start is not None:
                protein_effect_start = int(protein_effect_start)
            exons_effected = trans.get('exon')
            if exons_effected is not None:
                exons_effected = exons_effected.split('/')[0].split('-')
                if len(exons_effected) == 2:
                    exons_effected = range(int(exons_effected[0]), int(exons_effected[1]) + 1)
                exons_effected = [int(x) for x in exons_effected]
            transcript_list.append(
                TranscriptAnnotation(gene_name,
                                     trans_id,
                                     hgvsc_id,
                                     consequences,
                                     exons_effected,
                                     protein_id,
                                     protein_effect_start,
                                     protein_effect_end)
            )

        return Variant(variant_id, variant_class, variant_coordinates, canon_tx, transcript_list,
                       variant_coordinates.genotype)

    def _query_vep(self, variant_coordinates) -> dict:
        api_url = self._url % (verify_start_end_coordinates(variant_coordinates))
        r = requests.get(api_url, headers={'Content-Type': 'application/json'})
        if not r.ok:
            self._logging.error(f"Expected a result but got an Error for variant: {variant_coordinates.as_string()}")
            r.raise_for_status()
        results = r.json()
        if not isinstance(results, list):
            self._logging.error(results.get('error'))
            raise ConnectionError(f"Expected a result but got an Error. See log for details.")
        if len(results) > 1:
            self._logging.error([result.id for result in results])
            raise ValueError(f"Expected only one variant per request but received {len(results)} different variants.")
        return results[0]


class VariantAnnotationCache:

    def __init__(self, datadir: str):
        if not os.path.isdir(datadir):
            raise ValueError(f'datadir {datadir} must be an existing directory')
        self._datadir = datadir

    def get_annotations(self, variant_coordinates: VariantCoordinates) -> typing.Optional[Variant]:
        # TODO - EXPLAIN - removed `directory` param. A fine example of class state instead of parameter.
        fname = f'{variant_coordinates.as_string()}.pickle'
        fpath = os.path.join(self._datadir, fname)
        if os.path.exists(fpath):
            with open(fpath, 'rb') as fh:
                return pickle.load(fh)
        else:
            return None

    def store_annotations(self, variant_coordinates: VariantCoordinates, annotation: Variant):
        fname = f'{variant_coordinates.as_string()}.pickle'
        fpath = os.path.join(self._datadir, fname)
        with open(fpath, 'wb') as f:
            pickle.dump(annotation, f)


class CachingFunctionalAnnotator(FunctionalAnnotator):

    def __init__(self, cache: VariantAnnotationCache, fallback: FunctionalAnnotator):
        # TODO - EXPLAIN - removed output_directory. The caching annotator does not care about the cache internals.
        self._cache = hpotk.util.validate_instance(cache, VariantAnnotationCache, 'cache')
        self._fallback = hpotk.util.validate_instance(fallback, FunctionalAnnotator, 'fallback')

    def annotate(self, variant_coordinates: VariantCoordinates) -> Variant:
        hpotk.util.validate_instance(variant_coordinates, VariantCoordinates, 'variant_coordinates')
        annotations = self._cache.get_annotations(variant_coordinates)
        if annotations is not None:
            # we had cached annotations
            return annotations
        else:
            ann = self._fallback.annotate(variant_coordinates)
            self._cache.store_annotations(variant_coordinates, ann)
            return ann
