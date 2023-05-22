import abc
import logging
import pickle
import typing
import re
import requests
import os
# pyright: reportGeneralTypeIssues=false
from phenopackets import GenomicInterpretation
from ._models import FunctionalAnnotator, VariantCoordinateFinder
from ._variant_data import VariantCoordinates, TranscriptAnnotation, Variant

OUTPUT_DIR = os.getcwd()

class PhenopacketVariantCoordinateFinder(VariantCoordinateFinder):

    def __init__(self):
        self._logger = logging.getLogger(__name__)

    def find_coordinates(self, genomic_interp):
        if not isinstance(genomic_interp, GenomicInterpretation):
            raise ValueError(f"item must be a Phenopacket GenomicInterpretation but was type {type(genomic_interp)}")
        chrom, ref, alt, genotype = '', '', '', None
        start, end = 0, 0
        variant_descriptor = genomic_interp.variant_interpretation.variation_descriptor
        if len(variant_descriptor.vcf_record.chrom) == 0 and len(variant_descriptor.variation.copy_number.allele.sequence_location.sequence_id) != 0:
            ref = 'N'
            start = int(variant_descriptor.variation.copy_number.allele.sequence_location.sequence_interval.start_number.value)
            end = int(variant_descriptor.variation.copy_number.allele.sequence_location.sequence_interval.end_number.value)
            number = int(variant_descriptor.variation.copy_number.number.value)
            if number > 2:
                alt = '<DUP>'
            else:
                alt = '<DEL>'
            chrom = re.findall(r'NC_0000(\d{2}).\d\d', variant_descriptor.variation.copy_number.allele.sequence_location.sequence_id)[0]
            if chrom.startswith('0'):
                chrom = str(int(chrom))
            elif chrom == '23':
                chrom = 'X'
            elif chrom == '24':
                chrom = 'Y'
        elif len(variant_descriptor.vcf_record.chrom) != 0 and len(variant_descriptor.variation.copy_number.allele.sequence_location.sequence_id) == 0:
            ref = variant_descriptor.vcf_record.ref
            alt = variant_descriptor.vcf_record.alt
            start = int(variant_descriptor.vcf_record.pos) - 1
            end = int(variant_descriptor.vcf_record.pos) + abs(len(alt) - len(ref))
            chrom = variant_descriptor.vcf_record.chrom[3:]
        genotype = variant_descriptor.allelic_state.label
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
        #TODO: Verify <INS> are working correctly
    else:
        if len(vc.ref) == 0 or len(vc.alt) == 0:
            raise ValueError(f'Trimmed alleles are not yet supported!')
        if len(vc.ref) == 1 and len(vc.alt) != 1:
            # INS/DUP
            start = start + 1  # we must "trim"
            end = end - 1
            alt = vc.alt[1:]
            # 100 AC AGT
            # MNV
            

    return f'{chrom}:{start}-{end}/{alt}'


class VepFunctionalAnnotator(FunctionalAnnotator):

    def __init__(self):
        self._logging = logging.getLogger(__name__)
        self._url = 'https://rest.ensembl.org/vep/human/region/%s?LoF=1&canonical=1&domains=1&hgvs=1&mutfunc=1&numbers=1&protein=1&refseq=1&transcript_version=1&variant_class=1'

    def annotate(self, variant_coordinates: VariantCoordinates) -> Variant:
        api_url = self._url % (verify_start_end_coordinates(variant_coordinates))
        r = requests.get(api_url, headers={'Content-Type':'application/json'})
        if not r.ok:
            self._logging.error(f"Expected a result but got an Error for variant: {variant_coordinates.as_string()}")
            r.raise_for_status()
        results = r.json()
        if not isinstance(results, list):
            self._logging.error(results.get('error'))
            raise ConnectionError(f"Expected a result but got an Error. See log for details.")
        if len(results) > 1:
            self._logging.error([re.id for re in results])
            raise ValueError(f"Expected only one variant per request but received {len(results)} different variants.")
        variant = results[0]
        variant_id = variant.get('id')
        variant_class = variant.get('variant_class')
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
            transcript_list.append(TranscriptAnnotation(gene_name, trans_id, hgvsc_id, consequences, exons_effected,
                protein_id, protein_effect_start, protein_effect_end))
        return Variant(variant_id, variant_class, variant_coordinates, canon_tx, transcript_list, variant_coordinates.genotype)


class VariantAnnotationCache:

    def get_annotations(self, variant_coordinates: VariantCoordinates, directory) -> typing.Optional[Variant]:
        if os.path.exists(f'{directory}/{variant_coordinates.as_string()}.pickle'):
            with open(f'{directory}/{variant_coordinates.as_string()}.pickle', 'rb') as f:
                variant = pickle.load(f)
        else:
            variant = None
        return variant

    def store_annotations(self, variant_coordinates: VariantCoordinates, annotation: Variant, directory):
        with open(f'{directory}/{variant_coordinates.as_string()}.pickle', 'wb') as f:
            pickle.dump(annotation, f)
        return None


class CachingFunctionalAnnotator(FunctionalAnnotator):

    def __init__(self, output_directory: str, cache: VariantAnnotationCache, fallback: FunctionalAnnotator):
        if not os.path.exists(output_directory):
            raise ValueError(f"Invalid path: {output_directory}")
        directory = os.path.join(output_directory, 'annotations')
        os.makedirs(directory, exist_ok=True)
        self._output = directory
        if not isinstance(cache, VariantAnnotationCache):
            raise ValueError(f"cache must be type VariantAnnotationCache but was type {type(cache)}")
        self._cache = cache
        if not isinstance(fallback, FunctionalAnnotator):
            raise ValueError(f"fallback must be type FunctionalAnnotator but was type {type(fallback)}")
        self._fallback = fallback

    def annotate(self, variant_coordinates: VariantCoordinates) -> Variant:
        if not isinstance(variant_coordinates, VariantCoordinates):
            raise ValueError(f"variant_coordinates must be type VariantCoordinates but was type {type(variant_coordinates)}")
        annotations = self._cache.get_annotations(variant_coordinates, self._output)
        if annotations is not None:
            # we had cached annotations
            return annotations
        else:
            ann = self._fallback.annotate(variant_coordinates)
            self._cache.store_annotations(variant_coordinates, ann, self._output)
            return ann