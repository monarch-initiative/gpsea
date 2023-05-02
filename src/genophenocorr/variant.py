import abc
import logging
import typing
import re
import requests
import os
from phenopackets import Phenopacket

URL = 'https://rest.ensembl.org/vep/human/hgvs/%s:%s?protein=1&variant_class=1&numbers=1&LoF=1&canonical=1&domains=1&hgvs=1&mutfunc=1&refseq=1&transcript_version=1'
STRUCT_URL = 'https://rest.ensembl.org/vep/human/region/%s:%s-%s:%s/%s?LoF=1&canonical=1&domains=1&hgvs=1&mutfunc=1&numbers=1&protein=1&refseq=1&transcript_version=1&variant_class=1&transcript_id=%s'
OUTPUT_DIR = os.getcwd()


class VariantCoordinates:
    """A representation of coordinates of sequence and symbolic variants.

    The breakend variants are not supported.
    """

    def __init__(self, chrom, start, end, ref, alt, change_length):
        # TODO(lnrekerle) - instance/type check
        # TODO - id?
        self._chrom = chrom
        self._start = start
        self._end = end
        self._ref = ref
        self._alt = alt
        self._change_length = change_length

    @property
    def chrom(self) -> str:
        """
        Get label of the chromosome/contig where the variant is located.
        """
        return self._chrom

    @property
    def start(self) -> int:
        """
        Get 0-based start coordinate (excluded) of the ref allele.
        """
        return self._start

    @property
    def end(self) -> int:
        """
        Get 0-based end coordinate (included) of the ref allele.
        """
        return self._end

    @property
    def ref(self) -> str:
        """
        Get reference allele (e.g. "A", "N"). The allele may be an empty string.
        """
        return self._ref

    @property
    def alt(self) -> str:
        """
        Get alternate allele (e.g. "A", "GG", "<DEL>"). The allele may be an empty string for sequence variants.
        The symbolic alternate allele follow the VCF notation and use the `<` and `>` characters
        (e.g. "<DEL>", "<INS:ME:SINE>").
        """
        return self._alt

    @property
    def change_length(self) -> int:
        """
        Get the change between the ref and alt alleles due to the variant presence. SNVs lead to change length of zero,
        deletions and insertions/duplications lead to negative and positive change lengths, respectively.
        """
        return self._change_length

    def is_structural(self) -> bool:
        """
        Return `True` if the variant coordinates use structural variant notation
        (e.g. `chr5  101 . N <DEL> .  .  SVTYPE=DEL;END=120;SVLEN=-10`)
        as opposed to the sequence/literal notation (`chr5  101 . NACGTACGTAC N`).
        """
        return len(self._alt) != 0 and self._alt.startswith('<') and self._alt.endswith('>')

    def __len__(self):
        """
        Get the number of bases on the ref allele that are affected by the variant.
        """
        return self._end - self._start

    def __eq__(self, other) -> bool:
        return isinstance(other, VariantCoordinates) \
            and self.alt == other.alt \
            and self.ref == other.ref \
            and self.chrom == other.chrom \
            and self.start == other.start \
            and self.end == other.end \
            and self.change_length == other.change_length

    def __str__(self) -> str:
        return f"VariantCoordinates(chrom={self.chrom}, " \
            f"start={self.start}, end={self.end}, " \
            f"ref={self.ref}, alt={self.alt}, " \
            f"change_length={self.change_length} "

    def __repr__(self) -> str:
        return str(self)

    def __hash__(self) -> int:
        return hash((self._chrom, self._start, self._end, self._ref, self._alt, self._change_length))


class VariantCoordinateFinder(metaclass=abc.ABCMeta):
    @abc.abstractmethod
    def find_coordinates(self, item) -> VariantCoordinates:
        pass


class PhenopacketVariantCoordinateFinder(VariantCoordinateFinder):

    def __init__(self):
        self._logger = logging.getLogger(__name__)

    def find_coordinates(self, item):
        if not isinstance(item, Phenopacket):
            raise ValueError(f"item must be a Phenopacket but was type {type(item)}")
        if len(item.interpretations) != 1:
            raise ValueError(f'Phenopackets with not exactly one interpretation are not yet supported.')
        genomic_interpretations = item.interpretations[0].diagnosis.genomicInterpretations
        variants_list = []
        # TODO(lnrekerle) - prepare a unit test for testing PhenopacketVariantCoordinateFinder
        # TODO(ielis & lnrekerle) - debug the code below
        for gi in genomic_interpretations:
            chrom, ref, alt, hgvs = '', '', '', ''
            start, end = 0, 0
            variant_descriptor = gi.variant_interpretation.variation_descriptor
            if len(variant_descriptor.vcf_record.chrom) == 0 and len(variant_descriptor.variation.copy_number.sequence_location.sequence_id) != 0:
                ref = '-'
                start = variant_descriptor.variation.copy_number.allele.sequence_location.sequence_interval.start_number.value
                end = variant_descriptor.variation.copy_number.allele.sequence_location.sequence_interval.start_number.value
                number = variant_descriptor.variation.copy_number.number.value
                if number > 2:
                    alt = '<DUP>'
                else:
                    alt = '<DEL>'
                chrom = re.findall(r'NC_0000(\d{2}).\d\d', variant_descriptor.variation.copy_number.sequence_location.sequence_id)[0]
                if chrom.startswith('0'):
                    chrom = str(int(chrom))
                elif chrom == '23':
                    chrom = 'X'
                elif chrom == '24':
                    chrom = 'Y'
            elif len(variant_descriptor.vcf_record.chrom) != 0 and len(variant_descriptor.variation.copy_number.sequence_location.sequence_id) == 0:
                expressions = variant_descriptor.expressions
                hgvs = ''
                for ex in expressions:
                    if ex.syntax == 'hgvs.g':
                        hgvs = ex.value.split(':')[1]
                ref = variant_descriptor.vcf_record.ref
                alt = variant_descriptor.vcf_record.alt
                #TODO: Verify these Starts and Ends are correct
                if len(ref) == len(alt) and len(ref) == 1:
                    start = int(variant_descriptor.vcf_record.pos) - 1
                    end = int(variant_descriptor.vcf_record.pos)
                elif len(ref) > len(alt):
                    start = int(variant_descriptor.vcf_record.pos) + 1
                    end = int(variant_descriptor.vcf_record.pos) + len(ref)
                elif len(ref) < len(alt):
                    if hgvs.endswith('dup'):
                        start = int(variant_descriptor.vcf_record.pos) + 1
                        end = int(variant_descriptor.vcf_record.pos) + len(alt) + 1
                    else:
                        start = int(variant_descriptor.vcf_record.pos)
                        end = int(variant_descriptor.vcf_record.pos) + len(alt)
                chrom = variant_descriptor.vcf_record.chrom[3:]
            variants_list.append(VariantCoordinates(chrom, start, end, ref, alt, len(alt) - len(ref)))

        # TODO(lnrekerle) - test if it works for the strange dup in ankyrin (phenopacket)
        return variants_list


class TranscriptAnnotation:
    """
    Class that represents results of the functional annotation of a variant with respect to single transcript of a gene.
    """

    def __init__(self, gene_id: str,
                 tx_id: str,
                 variant_effects,
                 affected_exons: list[int],
                 affected_protein: str,
                 protein_effect_start: int,
                 protein_effect_end: int):
        self._gene_id = gene_id
        self._tx_id = tx_id
        self._variant_effects = variant_effects
        self._affected_exons = affected_exons
        self._affected_protein = affected_protein
        self._protein_effect_location = [protein_effect_start, protein_effect_end]

    @property
    def gene_id(self) -> str:
        """
        Get the gene symbol (e.g. SURF1)
        """
        return self._gene_id

    @property
    def transcript_id(self) -> str:
        """
        Get the transcript identifier (e.g. NM_123456.7)
        """
        return self._tx_id

    @property
    def variant_effects(self):
        """
        Get a sequence of variant effects.
        """
        return self._variant_effects

    @property
    def overlapping_exons(self) -> list[int]:
        """
        Get a sequence of IDs of the exons that overlap with the variant.
        """
        return self._affected_exons

    @property
    def protein_affected(self) -> str:
        return self._affected_protein

    @property
    def protein_effect_location(self) -> list[int]:
        return self.protein_effect_location

class Variant:

    # NM_123456.7:c.1243T>C

    # Attributes:
    #  - id
    #  - list of TranscriptAnnotations

    def __init__(self, var_id, var_class, var_coordinates: VariantCoordinates, tx_annotations, genotype):
        self._id = var_id
        self._var_coordinates = var_coordinates
        self._var_class = var_class
        self._tx_annotations = tx_annotations
        self._genotype = genotype  # This is optional

    @property
    def variant_coordinates(self) -> VariantCoordinates:
        return self._var_coordinates

    @property
    def variant_string(self):
        return self._id

    @property
    def genotype(self):
        return self._genotype

    @property
    def tx_annotations(self) -> typing.Sequence[TranscriptAnnotation]:
        return self._tx_annotations

    @property
    def variant_class(self):
        return self._var_class


class FunctionalAnnotator(metaclass=abc.ABCMeta):

    @abc.abstractmethod
    def annotate(self, variant_coordinates: VariantCoordinates) -> typing.Sequence[TranscriptAnnotation]:
        pass


class VepFunctionalAnnotator(FunctionalAnnotator):

    def __init__(self):
        self._logging = logging.getLogger(__name__)
        self._hgvs_url = 'https://rest.ensembl.org/vep/human/hgvs/%s:%s?protein=1&variant_class=1&numbers=1&LoF=1&canonical=1&domains=1&hgvs=1&mutfunc=1&refseq=1&transcript_version=1'
        self._other_url = 'https://rest.ensembl.org/vep/human/region/%s:%s-%s:1/%s?LoF=1&canonical=1&domains=1&hgvs=1&mutfunc=1&numbers=1&protein=1&refseq=1&transcript_version=1&variant_class=1'

    def annotate(self, variant_coordinates: VariantCoordinates) -> typing.Sequence[Variant]:
        if variant_coordinates.hgvs is None:
            api_url = self._other_url % (variant_coordinates.chrom, variant_coordinates.start, variant_coordinates.end, variant_coordinates.ref)
        else:
            api_url = self._hgvs_url % (variant_coordinates.chrom, variant_coordinates.hgvs)
        r = requests.get(api_url, headers={'Content-Type':'application/json'})
        if not r.ok:
            r.raise_for_status()
        results = r.json()
        variant_list = []
        for variants in results:
            variant_id = results.get('id')
            variant_class = results.get('variant_class')
            transcript_list = []
            for trans in variants.get('transcript_consequences'):
                trans_id = trans.get('transcript_id')
                if not trans_id.startswith('NM'):
                    self._logging.info(f'We will not be using this transcript: {trans_id}')
                    continue
                consequences = trans.get('consequence_terms')
                gene_name = trans.get('gene_symbol')
                protein_id = trans.get('protein_id')
                protein_effect_start = trans.get('protein_start')
                protein_effect_end = trans.get('protein_end')
                exons_effected = trans.get('exon').split('/')[0].split('-')
                if len(exons_effected) == 2:
                    exons_effected = range(int(exons_effected[0]), int(exons_effected[1]) + 1)
                exons_effected = [int(x) for x in exons_effected]
                transcript_list.append(TranscriptAnnotation(gene_name, trans_id, consequences, exons_effected,
                protein_id, int(protein_effect_start), int(protein_effect_end)))
            variant_list.append(Variant(variant_id, variant_class, variant_coordinates, transcript_list, None))
        return variant_list




class VariantAnnotationCache:

    # TODO - implement the persistence strategy (e.g. pickle)

    def get_annotations(self, variant_coordinates: VariantCoordinates) -> typing.Optional[typing.Sequence[TranscriptAnnotation]]:
        # TODO - implement
        return None

    def store_annotations(self, variant_coordinates: VariantCoordinates, annotations: typing.Sequence[TranscriptAnnotation]):
        # TODO - implement
        pass


class CachingFunctionalAnnotator(FunctionalAnnotator):

    def __init__(self, cache: VariantAnnotationCache, fallback: FunctionalAnnotator):
        # TODO(lnrekerle) - instance/type check
        self._cache = cache
        self._fallback = fallback

    def annotate(self, variant_coordinates: VariantCoordinates) -> typing.Sequence[TranscriptAnnotation]:
        annotations = self._cache.get_annotations(variant_coordinates)
        if annotations is not None:
            # we had cached annotations
            return annotations
        else:
            ann = self._fallback.annotate(variant_coordinates)
            self._cache.store_annotations(variant_coordinates, ann)
            return ann