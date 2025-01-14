from ._base import SampleLabels
from ._gt import Genotypes, Genotype


class TestGenotypes:

    def test_genotypes(self):
        gts = Genotypes.from_mapping({
            SampleLabels('A'): Genotype.HETEROZYGOUS,
            SampleLabels('C'): Genotype.HEMIZYGOUS,
            SampleLabels('B'): Genotype.HOMOZYGOUS_ALTERNATE,
            SampleLabels('E'): Genotype.NO_CALL,
            SampleLabels('D'): Genotype.HOMOZYGOUS_REFERENCE,
        })
        assert gts.for_sample(SampleLabels('A')) == Genotype.HETEROZYGOUS
        assert gts.for_sample(SampleLabels('B')) == Genotype.HOMOZYGOUS_ALTERNATE
        assert gts.for_sample(SampleLabels('C')) == Genotype.HEMIZYGOUS
        assert gts.for_sample(SampleLabels('D')) == Genotype.HOMOZYGOUS_REFERENCE
        assert gts.for_sample(SampleLabels('E')) == Genotype.NO_CALL
        assert gts.for_sample(SampleLabels('X')) is None

    def test_iteration(self):
        genotypes = Genotypes.from_mapping({
            SampleLabels('A'): Genotype.HETEROZYGOUS,
            SampleLabels('C'): Genotype.HEMIZYGOUS,
            SampleLabels('D'): Genotype.HOMOZYGOUS_REFERENCE,
        })
        labels = tuple(label for label, _ in genotypes)
        gts = tuple(gt for _, gt in genotypes)

        assert len(labels) == len(gts) == len(genotypes)
        
        assert all(sample_labels.label in ('A', 'C', 'D') for sample_labels in labels)
        assert all(
            gt in (Genotype.HETEROZYGOUS, Genotype.HEMIZYGOUS, Genotype.HOMOZYGOUS_REFERENCE)
            for gt in gts
        )
