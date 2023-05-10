import abc
import logging
import pickle
import typing
import re
import requests
import os
from phenopackets import GenomicInterpretation

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

    def as_string(self) -> str:
        return f"{self.chrom}_{self.start}_{self.end}_{self.ref}_{self.alt}"

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

    def find_coordinates(self, genomic_interp):
        if not isinstance(genomic_interp, GenomicInterpretation):
            raise ValueError(f"item must be a Phenopacket GenomicInterpretation but was type {type(genomic_interp)}")
        chrom, ref, alt = '', '', ''
        start, end = 0, 0
        variant_descriptor = genomic_interp.variant_interpretation.variation_descriptor
        if len(variant_descriptor.vcf_record.chrom) == 0 and len(variant_descriptor.variation.copy_number.allele.sequence_location.sequence_id) != 0:
            ref = '-'
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
            start = int(variant_descriptor.vcf_record.pos) 
            end = int(variant_descriptor.vcf_record.pos) + abs(len(alt) - len(ref)) 
            chrom = variant_descriptor.vcf_record.chrom[3:]
        return VariantCoordinates(chrom, start, end, ref, alt, len(alt) - len(ref))


class TranscriptAnnotation:
    """
    Class that represents results of the functional annotation of a variant with respect to single transcript of a gene.
    """

    def __init__(self, gene_id: str,
                 tx_id: str,
                 hgvsc: typing.Optional[str],
                 variant_effects,
                 affected_exons: typing.Optional[list[int]],
                 affected_protein: str,
                 protein_effect_start: typing.Optional[int],
                 protein_effect_end: typing.Optional[int]):
        self._gene_id = gene_id
        self._tx_id = tx_id
        self._hgvsc_id = hgvsc
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
    def hgvsc_id(self):
        """
        Get the HGVS "coding-DNA" representation of the variant (e.g. NM_123456.7:c.9876G>T)
        """
        return self._hgvsc_id

    @property
    def variant_effects(self):
        """
        Get a sequence of variant effects. 
        Definitions of these can be found at: http://www.sequenceontology.org/
        """
        return self._variant_effects

    @property
    def overlapping_exons(self):
        """
        Get a sequence of IDs of the exons that overlap with the variant.
        """
        return self._affected_exons

    @property
    def protein_affected(self) -> str:
        """
        Get the protein ID that is affected by the alteration of this transcript (e.g. NP_037407.4)
        """
        return self._affected_protein

    @property
    def protein_effect_location(self) -> list[int]:
        """
        Get the start and end position on the protein sequence that the variant effects. (e.g. [1234, 1235])
        """
        return self.protein_effect_location

    def __str__(self) -> str:
        return f"TranscriptAnnotation(gene_id:{self.gene_id}," \
            f"transcript_id:{self.transcript_id}," \
            f"hgvsc_id:{self.hgvsc_id}," \
            f"variant_effects:{self.variant_effects}," \
            f"overlapping_exons:{self.overlapping_exons}," \
            f"protein_affected:{self.protein_affected}," \
            f"protein_effect_location:{self.protein_effect_location})"

    def __eq__(self, other) -> bool:
        return isinstance(other, TranscriptAnnotation) \
            and self.gene_id == other.gene_id \
            and self.hgvsc_id == other.hgvsc_id \
            and self.transcript_id == other.transcript_id \
            and self.variant_effects == other.variant_effects \
            and self.overlapping_exons == other.overlapping_exons \
            and self.protein_affected == other.protein_affected \
            and self.protein_effect_location == other.protein_effect_location

    def __repr__(self) -> str:
        return str(self)

    def __hash__(self) -> int:
        return hash((self.gene_id, self.hgvsc_id, self.transcript_id, self.overlapping_exons, self.variant_effects, self.protein_affected, self.protein_effect_location))

class Variant:
    """
    Class that represents results of the functional annotation of a variant with all included transcripts.
    """

    def __init__(self, var_id, var_class, var_coordinates: VariantCoordinates, tx_annotations, genotype):
        self._id = var_id
        self._var_coordinates = var_coordinates
        self._var_class = var_class
        self._tx_annotations = tx_annotations
        self._genotype = genotype  # This is optional

    @property
    def variant_coordinates(self) -> VariantCoordinates:
        """
        A representation of coordinates of a sequence and symbolic variant.
        """
        return self._var_coordinates

    @property
    def variant_string(self):
        """
        A readable representation of the variant's coordinates. 
        Format - "Chromosome_Start_Reference/Alternative" or
                 "Chromosome_Start_StructuralType"
        """
        return self._id

    @property
    def genotype(self):
        """
        Optional parameter. Required for recessive tests. 
        Possible values: Heterozygous, Homozygous, Hemizygous
        """
        return self._genotype

    @property
    def tx_annotations(self) -> typing.Sequence[TranscriptAnnotation]:
        """
        A collection of TranscriptAnnotations that each represent results of the functional annotation 
        of a variant with respect to single transcript of a gene.
        """
        return self._tx_annotations

    @property
    def variant_class(self):
        "The variant class. (e.g. Duplication, SNV, Deletion, etc.)"
        return self._var_class

    def __str__(self) -> str:
        return f"Variant(variant_coordinates:{str(self.variant_coordinates)}," \
            f"variant_string:{self.variant_string}," \
            f"genotype:{self.genotype}," \
            f"tx_annotations:{self.tx_annotations}," \
            f"variant_class:{self.variant_class})"


    def __eq__(self, other) -> bool:
        return isinstance(other, Variant) \
            and self.variant_string == other.variant_string \
            and self.genotype == other.genotype \
            and self.variant_class == other.variant_class \
            and self.variant_coordinates == other.variant_coordinates 
            #and self.tx_annotations == other.tx_annotations

    def __repr__(self) -> str:
        return str(self)

    def __hash__(self) -> int:
        return hash((self.variant_coordinates, self.variant_string, self.variant_class, self.genotype, self.tx_annotations))


class FunctionalAnnotator(metaclass=abc.ABCMeta):

    @abc.abstractmethod
    def annotate(self, variant_coordinates: VariantCoordinates) -> Variant:
        pass



class VepFunctionalAnnotator(FunctionalAnnotator):

    def __init__(self):
        self._logging = logging.getLogger(__name__)
        self._url = 'https://rest.ensembl.org/vep/human/region/%s?LoF=1&canonical=1&domains=1&hgvs=1&mutfunc=1&numbers=1&protein=1&refseq=1&transcript_version=1&variant_class=1'

    def _verify_start_end_coordinates(self, variant_coordinates: VariantCoordinates):
        if len(variant_coordinates.alt) > len(variant_coordinates.ref) and variant_coordinates.alt != '<DEL>' and variant_coordinates.alt != '<DUP>':
            return f"{variant_coordinates.chrom}:{variant_coordinates.start + 1}-{variant_coordinates.start}/{variant_coordinates.alt[1:]}"
        else:
            return f"{variant_coordinates.chrom}:{variant_coordinates.start}-{variant_coordinates.end}/{variant_coordinates.alt.replace('<', '').replace('>', '')}"
        

    def annotate(self, variant_coordinates: VariantCoordinates) -> Variant:

        api_url = self._url % (self._verify_start_end_coordinates(variant_coordinates))
        r = requests.get(api_url, headers={'Content-Type':'application/json'})
        if not r.ok:
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
        transcript_list = []
        for trans in variant.get('transcript_consequences'):
            trans_id = trans.get('transcript_id')
            if not trans_id.startswith('NM'):
                self._logging.info(f'We will not be using this transcript: {trans_id}')
                continue
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
        return Variant(variant_id, variant_class, variant_coordinates, transcript_list, None)


class VariantAnnotationCache:

    # TODO - implement the persistence strategy (e.g. pickle)

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