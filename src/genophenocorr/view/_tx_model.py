import requests
import abc
import os
import json
from typing import Dict



# e.g., https://rest.variantvalidator.org/VariantValidator/tools/gene2transcripts/NM_000518.4
URL_SCHEME = "rest.variantvalidator.org/VariantValidator/tools/gene2transcripts/%s?content-type=application%%2Fjson"
NOT_AVAILABLE_INT = -42
NOT_AVAILABLE_STR = "n/a"


class VvTranscriptModel:
    """This class encapsulates the results of a call to VariantValidator to retrieve information about a transcript.

    The JSON that is returned by VariantValidator includes

        - current_name: e.g., protein tyrosine phosphatase non-receptor type 11
        - current_symbol: e.g., PTPN11
        - hgnc: e.g., HGNC:9644
        - requested_symbol: NM_002834.5 - this should be the same as the transcript we request
        - transcripts: list of obects with annotations that we represent using the class VvTranscript

    :param vval_json_data: A JSON object retrieved from VariantValidator
    :type vval_json_data: Dict
    """
    def __init__(self, vval_json_data) -> None:
        self._name = vval_json_data.get("current_name", NOT_AVAILABLE_STR)
        self._symbol = vval_json_data.get("current_symbol", NOT_AVAILABLE_STR)
        self._hgnc = vval_json_data.get("hgnc", NOT_AVAILABLE_STR)
        self._requested_symbol = vval_json_data.get("requested_symbol", NOT_AVAILABLE_STR)
        transcripts = vval_json_data.get("transcripts")
        self._transcripts = []
        if transcripts is None:
            raise ValueError("Could not retrieve transcripts")
        for transcript in transcripts:
            self._transcripts.append(VvTranscript(json_data=transcript))
        print(f"[INFO] retrieved {len(self._transcripts)} transcript models from VariantValidator")


    @property
    def name(self):
        return self._name

    @property
    def symbol(self):
        return self._symbol

    @property
    def hgnc(self):
        return self._hgnc

    @property
    def requested_symbol(self):
        return self._requested_symbol

    def get_transcripts(self):
        return self._transcripts

    def get_most_relevant_transcript(self):
        """
        return the higest priority transcript
        """
        if len(self._transcripts) == 0:
            raise ValueError("Attempt to sort empty transcript list")
        sorted_tr = sorted(self._transcripts, key=lambda tr: tr.get_priority_level())
        return sorted_tr[0]

    def __str__(self):
        n_transcripts = len(self._transcripts)
        return f"[VvTranscriptModel] {self._symbol}: {self._name} ({self._hgnc}); requested {self._requested_symbol}. Transcripts: n={n_transcripts}"



class VvTranscript:
    def __init__(self, json_data) -> None:
        if not "annotations" in json_data:
            # Should never happen
            raise ValueError("Malformed VariantValidator response - could not find annotations in JSON; these were the keys {';'.join(json_data.keys())}")
        annotations = json_data.get("annotations")
        self._coding_start = json_data.get("coding_start", NOT_AVAILABLE_INT)
        self._coding_end = json_data.get("coding_end", NOT_AVAILABLE_INT)
        self._reference = json_data.get("reference", NOT_AVAILABLE_STR)
        self._translation = json_data.get("translation", NOT_AVAILABLE_STR)
        self._length = json_data.get("length", NOT_AVAILABLE_INT)
        self._chromosome = annotations.get("chromosome", NOT_AVAILABLE_STR)
        self._ensembl_select = annotations.get("ensembl_select", False)
        self._mane_plus_clinical = annotations.get("mane_plus_clinical", False)
        self._mane_select = annotations.get("mane_select", False)
        self._refseq_select = annotations.get("refseq_select", False)
        db_xref = annotations.get("db_xref", {}) # we expect a dictionary with keys such as 'ncbigene'
        # db_xref v={'CCDS': 'CCDS7753.1', 'ensemblgene': None, 'hgnc': 'HGNC:4827', 'ncbigene': '3043', 'select': 'MANE'}
        # todo add optional entries from db xref
        if not 'genomic_spans' in json_data:
            # Should never happen
            raise ValueError("Malformed VariantValidator response - could not find genomic_spans in JSON; these were the keys {';'.join(json_data.keys())}")
        from pprint import pprint
        # biosequence is a string such as  "NC_000012.11"
        # to get hg38, we want  "NC_0000??.12", i.e., version 12
        genomic_spans = json_data.get('genomic_spans')
        if genomic_spans is None:
            raise ValueError("Could not retrieve genomic spans from VariantValidator JSON data")
        biosequence_list = [biosequence_id for biosequence_id in genomic_spans if biosequence_id.endswith(".12")]
        if len(biosequence_list) == 0:
            print(f"Could not retrieve hg38 genomic span from VariantValidator JSON data")
        elif len(biosequence_list) > 1:
            raise ValueError(f"Retrieved multiple hg38 genomic span elements from VariantValidator JSON data")
        else:
            hg38_biosequence_id = biosequence_list[0]
            seq_d = genomic_spans.get(hg38_biosequence_id)
            self._start_position = seq_d.get("start_position", NOT_AVAILABLE_INT)
            self._end_position = seq_d.get("end_position", NOT_AVAILABLE_INT)
            self._orientation = seq_d.get("orientation", NOT_AVAILABLE_INT)
            self._total_exons = seq_d.get("total_exons", NOT_AVAILABLE_INT)
            exon_structure = seq_d.get("exon_structure")
            if exon_structure is None:
                raise ValueError(f"Could not retrieve exon_structure from VariantValidator JSON data")
            self._exon_list = []
            for exon in exon_structure:
                self._exon_list.append(VvExon(json_data=exon))
            print(f"Extracted {len(self._exon_list)} exons")


    @property
    def reference(self):
        """TODO, can this be renamed reference_sequence_accession_id
        """
        return self._reference


    def get_exon_list(self):
        return self._exon_list

    def get_priority_level(self):
        """return a numerical priority level, with 1=mostt relevant and 5=less relevant

        The priority is based on the assessment mane_plus_clinical > mane_select > ensembl_select > refseq_select > other

        :returns: prioritz level from 1 (highest) to 5 (lowest)
        :rtype: int
        """
        if self._mane_plus_clinical:
            return 1
        elif self._mane_select:
            return 2
        elif self._ensembl_select:
            return 3
        elif self._refseq_select:
            return 4
        else:
            return 5


    def __str__(self):
        priority_levels = []
        if self._mane_plus_clinical:
            priority_levels.append("MANE-plus clinical")
        elif self._mane_select:
            priority_levels.append("MANE select")
        elif self._ensembl_select:
            priority_levels.append("ENSEMBL select")
        elif self._refseq_select:
            priority_levels.append("RefSeq select")
        if len(priority_levels) == 0:
            priority_levels.append("no priority level")
        priority = ', '.join(priority_levels)
        return f"[VvTranscript] {self._reference}-CDS: {self._coding_start}-{self._coding_end}, chrom: {self._chromosome} length {self._length} nt transl {self._translation}: {priority}"


class VvExon:

    def __init__(self, json_data) -> None:
        """Represents data about one exon extracted from the exon_structure array of VariantValidator JSON

        Contents of the JSON dictionary

            - cigar: e.g., "179="
            - exon_number
            - genomic_start
            - genomic_end
            - transcript_start
            - transcript_end

        :param json_data: JSON dictionary of data about an exon
        :type json_data: Dict[str,str]
        """
        self._cigar = json_data.get("cigar", NOT_AVAILABLE_STR)
        self._exon_number = json_data.get("exon_number")
        self._genomic_start = json_data.get("genomic_start")
        self._genomic_end  = json_data.get("genomic_end")
        self._transcript_start  = json_data.get("transcript_start")
        self._transcript_end  = json_data.get("transcript_end")

    @property
    def exon_number(self):
        return self._exon_number

    @property
    def genomic_start(self):
        return self._genomic_start

    @property
    def genomic_end(self):
        return self._genomic_end

    @property
    def transcript_start(self):
        return self._transcript_start

    @property
    def transcript_end(self):
        return self._transcript_end

    def __str__(self):
        return f"[VvExon] exon {self.exon_number}: {self.transcript_start}-{self._transcript_end}"



class TxRetriever(metaclass=abc.ABCMeta):

    @abc.abstractmethod
    def retrieve(self, transcript:str):
        """
        Retrieve a model of the transcript. Subclasses either call VariantValidator or another service, or load a JSON model from file to enable testing.

        :param transcript: transcript identifier
        :type transcript: str
        """
        pass

    def json_to_model(self, json_data):
        if not isinstance(json_data, dict):
            raise ValueError(f"Expecting json_data argument to be a Python dictionary but got {type(json_data)}")
        # for testing
        current_name = json_data.get("current_name", NOT_AVAILABLE_STR)
        current_symbol = json_data.get("current_symbol", NOT_AVAILABLE_STR)
        hgnc = json_data.get("hgnc", NOT_AVAILABLE_STR)
        requested_symbol = json_data.get("requested_symbol", NOT_AVAILABLE_STR)

        return json_data


class TxJsonFileRetriever(TxRetriever):

    def __init__(self, directory_location) -> None:
        super().__init__()
        if not os.path.isdir(directory_location):
            raise FileNotFoundError(f"Could not find directory at {directory_location}")
        self._dir = directory_location

    def retrieve(self, transcript:str):
        if not transcript.startswith("NM_"):
            raise ValueError(f"ERROR: genophenocorr only accepts RefGene transcripts (that begin with \"NM_\"). Could not process transcript {transcript}")
        fname = os.path.join(self._dir, f"{transcript}.json")
        if not os.path.isfile(fname):
            raise FileNotFoundError(f"Could not find transcript JSON file at {fname}")
        fh = open(fname)
        data = json.load(fh)
        fh.close()
        return self.json_to_model(data)



