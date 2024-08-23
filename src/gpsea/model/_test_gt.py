from ._gt import Genotypes, Genotype


class TestGenotypes:

    def test_genotypes(self):
        gts = Genotypes.from_mapping({
            'A': Genotype.HETEROZYGOUS,
            'C': Genotype.HEMIZYGOUS,
            'B': Genotype.HOMOZYGOUS_ALTERNATE,
            'E': Genotype.NO_CALL,
            'D': Genotype.HOMOZYGOUS_REFERENCE,
        })
        assert gts.for_sample('A') == Genotype.HETEROZYGOUS
        assert gts.for_sample('B') == Genotype.HOMOZYGOUS_ALTERNATE
        assert gts.for_sample('C') == Genotype.HEMIZYGOUS
        assert gts.for_sample('D') == Genotype.HOMOZYGOUS_REFERENCE
        assert gts.for_sample('E') == Genotype.NO_CALL
        assert gts.for_sample('X') is None

    def test_iteration(self):
        gts = Genotypes.from_mapping({
            'A': Genotype.HETEROZYGOUS,
            'C': Genotype.HEMIZYGOUS,
            'D': Genotype.HOMOZYGOUS_REFERENCE,
        })
        items = [(sample_id, gt) for sample_id, gt in gts]
        assert items[0] == ('A', Genotype.HETEROZYGOUS)
        assert items[1] == ('C', Genotype.HEMIZYGOUS)
        assert items[2] == ('D', Genotype.HOMOZYGOUS_REFERENCE)
