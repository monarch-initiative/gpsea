import pytest

from stairval.notepad import create_notepad
from phenopackets.schema.v2.core.base_pb2 import OntologyClass

from ._phenopacket import PhenopacketOntologyTermOnsetParser

class TestPhenopacketOntologyTermOnsetParser:

    @pytest.fixture(scope='class')
    def parser(self) -> PhenopacketOntologyTermOnsetParser:
        return PhenopacketOntologyTermOnsetParser.default_parser()

    @pytest.mark.parametrize(
        'curie, days, is_postnatal',
        [
            ('HP:0030674', 140., False),   # Antenatal onset
            ('HP:0011460', 38.5, False),   # Embryonal onset
            
            ('HP:0003577', 0., True),      # Congenital onset
            ('HP:0003623', 14.5, True),    # Neonatal onset
            ('HP:0410280', 2936.5, True),  # Pediatric onset
            ('HP:0011462', 10227.0, True), # Young adult onset

            ('HP:0003584', 25567.5, True), # Late onset
        ]
    )
    def test_process(
        self,
        curie: str,
        days: float,
        is_postnatal: bool,
        parser: PhenopacketOntologyTermOnsetParser,
    ):
        notepad = create_notepad(label='Top')
        oc = OntologyClass(id=curie, label='Whatever')
        
        age = parser.process(oc, notepad)
        
        assert age is not None
        assert age.is_postnatal == is_postnatal
        assert age.days == pytest.approx(days)

        assert not notepad.has_errors_or_warnings(include_subsections=True)
