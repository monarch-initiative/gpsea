import os
import json
import typing
import unittest

import hpotk
import requests

from genophenocorr.model import VariantCoordinates, TranscriptInfoAware, TranscriptCoordinates
from genophenocorr.model.genome import GenomeBuild, GenomicRegion, Strand, Contig, transpose_coordinate
#from ._api import VariantCoordinateFinder, TranscriptCoordinateService


RESPONSE_FILENAME = os.path.join(os.path.dirname(__file__), 'test_data', 'test_vv_response.json')


class TestVvPreprocessing(unittest.TestCase):
    """Test class to understand and document the structure of the JSON returned by VariantValidator queries for transcript coordinates.
    """
    @classmethod
    def setUpClass(cls) -> None:
        f = open(RESPONSE_FILENAME)
        cls._response_d = json.load(f,)

    def test_response_is_dictionary(self):
        self.assertIsNotNone(self._response_d)
        self.assertTrue(isinstance(self._response_d, dict))


    def test_keys(self):
        expected_keys = {"current_name", "current_symbol", "hgnc", "previous_symbol", "requested_symbol", "transcripts" }
        self.assertEqual(len(expected_keys), len(self._response_d.keys()))
        for key in expected_keys:
            self.assertTrue(key in self._response_d)
        
    def test_current_name(self):
        expected = "sulfite oxidase"
        # variant_identifier = list(self._response_d.keys())[0]
        self.assertEqual(expected, self._response_d.get("current_name"))

    def test_current_symbol(self):
        expected = "SUOX"
        self.assertEqual(expected, self._response_d.get("current_symbol"))

    def test_hgnc(self):
        expected = "HGNC:11460"
        self.assertEqual(expected, self._response_d.get("hgnc"))

    def test_previous_symbol(self):
        expected = ""
        self.assertEqual(expected, self._response_d.get("previous_symbol"))

    def test_requested_symbol(self):
        expected = "NM_001032386.2"  ## this is the transcript ID we requested from variantvalidator
        self.assertEqual(expected, self._response_d.get("requested_symbol"))

    def test_transcripts(self):
        transcripts = self._response_d.get("transcripts")
        self.assertTrue(isinstance(transcripts, list))
        self.assertEqual(6, len(transcripts))
        for t in transcripts:
            print(t)

    def test_transcript_zero(self):
        transcripts = self._response_d.get("transcripts")
        t0 = transcripts[0]
        self.assertTrue(isinstance(t0, dict))
        expected_keys = {'annotations', 'coding_end', 'coding_start', 'description', 'genomic_spans', 'length', 'reference', 'translation'}
        self.assertEqual(len(t0.keys()), len(expected_keys))
        for key in t0.keys():
            self.assertTrue(key in expected_keys)
        coding_start = t0.get("coding_start")
        self.assertEqual(75, coding_start)
        coding_end = t0.get("coding_end")
        self.assertEqual(1712, coding_end)
        length = t0.get("length")
        self.assertEqual(2210, length)
        reference = t0.get("reference")
        self.assertEqual("NM_001032387.2", reference)
        translation = t0.get("translation")
        self.assertEqual("NP_001027559.1", translation)
        genomic_spans = t0.get("genomic_spans")
        self.assertTrue(isinstance(genomic_spans, dict))
        #Homo sapiens chromosome 12, GRCh37.p13
        self.assertTrue("NC_000012.11" in genomic_spans)
        # Homo sapiens chromosome 12, GRCh38.p14 
        self.assertTrue("NC_000012.12" in genomic_spans)
        hg38span = genomic_spans.get("NC_000012.12")
        end_position = hg38span.get("end_position")
        self.assertEqual(56005525, end_position)
        exon_structure = hg38span.get("exon_structure")
        self.assertTrue(isinstance(exon_structure, list))
        self.assertEqual(4, len(exon_structure))
        total_exons = hg38span.get("total_exons")
        self.assertEqual(4, total_exons)