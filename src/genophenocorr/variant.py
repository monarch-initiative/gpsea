import abc
import typing
import re
import requests
import pandas as pd
import pickle
import os

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
        Return `True` if the variant coordinates describe a structural variant.
        """
        return len(self._alt) != 0 and self._alt.startswith('<') and self._alt.endswith('>')

    def __len__(self):
        """
        Get the number of bases on the ref allele that are affected by the variant.
        """
        return self._end - self._start

    # TODO - eq, hash, repr, str


class TranscriptAnnotation:
    """
    Class that represents results of the functional annotation of a variant with respect to single transcript of a gene.
    """

    def __init__(self, gene_id: str,
                 tx_id: str,
                 variant_effects,
                 affected_exons: list[int]):
        self._gene_id = gene_id
        self._tx_id = tx_id
        self._variant_effects = variant_effects
        self._affected_exons = affected_exons

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


class FunctionalAnnotator(metaclass=abc.ABCMeta):

    @abc.abstractmethod
    def annotate(self, variant_coordinates: VariantCoordinates) -> typing.Sequence[TranscriptAnnotation]:
        pass


class VepFunctionalAnnotator(FunctionalAnnotator):

    def __init__(self):
        self._seq_url = 'https://rest.ensembl.org/vep/human/hgvs/%s:%s?protein=1&variant_class=1&numbers=1&LoF=1&canonical=1&domains=1&hgvs=1&mutfunc=1&refseq=1&transcript_version=1'
        self._structural_url = 'https://rest.ensembl.org/vep/human/region/%s:%s-%s:%s/%s?LoF=1&canonical=1&domains=1&hgvs=1&mutfunc=1&numbers=1&protein=1&refseq=1&transcript_version=1&variant_class=1&transcript_id=%s'

    def annotate(self, variant_coordinates: VariantCoordinates) -> typing.Sequence[TranscriptAnnotation]:
        if variant_coordinates.is_structural():
            # TODO(lnrekerle) implement:
            #  - parametrize & call the structural URL
            #  - process the results into a sequence of TranscriptAnnotation
            #  - return the results
            return []
        else:
            # TODO(lnrekerle) implement:
            #  - parametrize & call the sequence URL
            #  - process the results into a sequence of TranscriptAnnotation
            #  - return the results
            return []


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


class Variant:

    # NM_123456.7:c.1243T>C

    # Attributes:
    #  - id
    #  - list of TranscriptAnnotations

    def __init__(self, var_id, var_coordinates: VariantCoordinates, tx_annotations, genotype):
        self._id = var_id
        self._var_coordinates = var_coordinates
        self._tx_annotations = tx_annotations
        self._genotype = genotype  # This is optional
        # self._transcript = transcript
        # self._genoInterp = genoInterp
        # self._pickled_dir = pickled_dir
        # self._varInterp = genoInterp.variant_interpretation.variation_descriptor
        # self._structural = False
        # self._variant_json = self.__find_variant()
        # self.__set_canon(transcript)

    def __find_variant(self):
        self._allelic_state = self._varInterp.allelic_state.label
        if len(self._varInterp.structural_type.id) != 0:
            self._structural = True
            return self.__structure_variant()
        else:
            chrom = self._varInterp.vcf_record.chrom
            contig = re.sub(r'[^0-9MXY]', '', chrom)
            if len(contig) == 0 or (contig.isdigit() and (int(contig) == 0 or int(contig) >= 24)):
                ## Chromosome can only be values 1-23, X, Y, or M
                raise ValueError(f"Chromosome unknown: {chrom}")
            HGVS = self._varInterp.expressions[1]
            if self._pickled_dir is not None:
                if not os.path.isdir(os.path.join(OUTPUT_DIR, self._pickled_dir)) and self._pickled_dir is not None:
                    os.mkdir(os.path.join(OUTPUT_DIR, self._pickled_dir))
                if os.path.isfile(os.path.join(OUTPUT_DIR, f'{self._pickled_dir}/{HGVS}.pkl')) and self._pickled_dir is not None:
                    with open(os.path.join(OUTPUT_DIR, f'{self._pickled_dir}/{HGVS}.pkl'), 'rb') as f:
                        varJson = pickle.load(f)
                        return varJson
            var_input = HGVS.value.split(':')[1]
            api_url = URL % (chrom, var_input)
            r = requests.get(api_url, headers={'Content-Type':'application/json'})
            if not r.ok:
                r.raise_for_status()
            varJson = r.json()
            if self._pickled_dir is not None:
                with open(os.path.join(OUTPUT_DIR, f'{self._pickled_dir}/{HGVS}.pkl'), 'wb') as f:
                    pickle.dump(varJson[0], f)
            return varJson[0]

    def __structure_variant(self):
        chrom = self._varInterp.description.split('q')[0]
        start = self._varInterp.variation.copy_number.allele.sequence_location.sequence_interval.start_number.value
        end = self._varInterp.variation.copy_number.allele.sequence_location.sequence_interval.end_number.value
        if self._varInterp.variation.copy_number.number.value == 1:
            copy_type = 'DEL'
            times = 1
            finalConsequence = 'copy_number_decrease'
        elif self._varInterp.variation.copy_number.number.value >= 3:
            copy_type = 'DUP'
            times = self._varInterp.variation.copy_number.number.value - 2
            finalConsequence = 'copy_number_increase'
        if self._pickled_dir is not None:
            if not os.path.isdir(os.path.join(OUTPUT_DIR, self._pickled_dir)):
                os.mkdir(os.path.join(OUTPUT_DIR, self._pickled_dir))
            if os.path.isfile(os.path.join(OUTPUT_DIR, f'{self._pickled_dir}/{self._varInterp.description}.pkl')) and self._pickled_dir is not None:
                with open(os.path.join(OUTPUT_DIR, f'{self._pickled_dir}/{self._varInterp.description}.pkl'), 'rb') as f:
                    varJson = pickle.load(f)
                    return varJson
        api_url = STRUCT_URL % (chrom,start,end,times,copy_type,self._transcript)
        r = requests.get(api_url, headers={ "Content-Type" : "application/json"})
        if not r.ok:
            r.raise_for_status()
        results = r.json()
        varJson = results[0]
        if varJson.get('transcript_consequences') is None:
            varJson['transcript_consequences'] = [{
                'canonical': 1,
                'consequence_terms':[],
                'gene_symbol': self._varInterp.gene_context.symbol
            }]
        varJson['transcript_consequences'][0]['consequence_terms'].append('copy_number_change')
        varJson['transcript_consequences'][0]['consequence_terms'].append(finalConsequence)
        if self._pickled_dir is not None:
            with open(os.path.join(OUTPUT_DIR, f'{self._pickled_dir}/{self._varInterp.description}.pkl'), 'wb') as f:
                pickle.dump(varJson, f)
        return varJson

    # @property
    # def variant_json(self):
    #     return self._variant_json

    @property
    def variant_coordinates(self) -> VariantCoordinates:
        return self._var_coordinates

    @property
    def variant_string(self):
        # __str__
        if self._structural:
            return self._varInterp.description
        else:
            return self.variant_json.get('id')

    @property
    def genotype(self):
        return self._genotype
        # if self._allelic_state is not None:
        #     return self._allelic_state
        # else:
        #     return None

    @property
    def tx_annotations(self) -> typing.Sequence[TranscriptAnnotation]:
        return self._tx_annotations

    # @property
    # def genomic_start(self):
    #     return self._variant_json.get('start')
    #
    # @property
    # def genomic_end(self):
    #     return self._variant_json.get('end')

    # TODO - decide if we want to store SNV, SV, or any other classes
    #  perhaps make it optional
    # @property
    # def variant_class(self):
    #     return self._variant_json.get('variant_class')

    @property
    def variant_types(self):
        return self._canonical.get('consequence_terms')
         
    @property
    def transcript(self):
        return self._transcript

    @property
    def gene_name(self):
        return self._canonical.get('gene_symbol')

    @property
    def effected_protein(self):
        return self._canonical.get('protein_id')

    @property
    def protein_effect_location(self):
        return self._canonical.get('protein_end')

    @property
    def protein_effect(self):
        return self._canonical.get('hgvsp').split(':')[1]

    @property
    def exons_effected(self):
        if self._canonical.get('exon') is not None:
            if '-' in self._canonical.get('exon'):
                exon_range = self._canonical.get('exon').split('/')[0].split('-')
                exons = [x for x in range(int(exon_range[0]), int(exon_range[1]) + 1)]
            else:
                exons = [int(self._canonical.get('exon').split('/')[0])]
            return exons
        else:
            return []


    def list_all_variant_effects_df(self):
        allSeries = []
        varID = self._variant_json.get('id')
        varClass = self._variant_json.get('variant_class')
        for trans in self._variant_json.get('transcript_consequences'):
            transID = trans.get('transcript_id')
            if not transID.startswith('NM'):
                continue
            consequence = trans.get('consequence_terms')
            geneName = trans.get('gene_symbol')
            proteinID = trans.get('protein_id')
            exonEffected = trans.get('exon').split('/')[0]
            proteinEffect = trans.get('hgvsp').split(':')[1]
            siftScore = trans.get('sift_score')
            polyphenScore = trans.get('polyphen_score')
            if trans.get('canonical') == 1:
                canon = True
            else:
                canon = False
            allSeries.append(pd.Series([consequence, geneName, proteinID, proteinEffect, exonEffected, siftScore, polyphenScore, canon], name = transID, index=['Variant Type', 'Gene Name', 'Protein ID', 'Effect on Protein', 'Exon Effected', 'Sift Score', 'Polyphen Score', 'Current Transcript']))
        if len(allSeries) > 1:
            results = pd.concat(allSeries, axis=1)
            results = results.transpose()
            results.Name = varID + ' - Class: ' + varClass
        else:
            results = allSeries[0]
            results.name = varID + ' - Class: ' + varClass
        return results 
            
    def __set_canon(self, transcript = None):
        self._canonical = None
        if transcript is not None:
            self._transcript = transcript
            if self._variant_json.get('transcript_consequences') is not None:
                for trans in self._variant_json.get('transcript_consequences'):
                    if trans.get('transcript_id') == transcript:
                        self._canonical = trans
                if self._canonical is None:
                    self._canonical = self._variant_json.get('transcript_consequences')[0]
            else:
                raise ValueError(f"{self.variant_string} has no transcript consequences:\n {self._variant_json}")
        else:
            if self._variant_json.get('transcript_consequences') is not None:
                for trans in self._variant_json.get('transcript_consequences'):
                    if trans.get('canonical') == 1:
                        self._transcript = trans.get('transcript_id')
                        self._canonical = trans
                if self._canonical is None:
                    self._canonical = self._variant_json.get('transcript_consequences')[0]
            else:
                raise ValueError(f"{self.variant_string} has no transcript consequences:\n {self._variant_json}")
        return None