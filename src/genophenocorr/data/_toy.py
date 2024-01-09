from hpotk import TermId

from genophenocorr.model import *
from genophenocorr.model.genome import Contig, GenomicRegion, Region, Strand

CONTIG = Contig('1', 'GB_ACC', 'REFSEQ_NAME', 'UCSC_NAME', 1_000)


def make_region(start: int, end: int) -> GenomicRegion:
    return GenomicRegion(CONTIG, start, end, Strand.POSITIVE)


def get_toy_cohort() -> Cohort:
    # TODO - Lauren - implement sample patients from the terms below
    # - Arachnodactyly HP:0001166
    # - Focal clonic seizure HP:0002266
    # - Perimembranous ventricular septal defect HP:0011682
    # - Hepatosplenomegaly HP:0001433
    # - Tubularization of Bowman capsule HP:0032648
    # - Intercostal muscle weakness HP:0004878
    # - Enuresis nocturna HP:0010677
    # - Spasticity HP:0001257
    # - Chronic pancreatitis HP:0006280

    arachnodactyly_T = Phenotype(TermId.from_curie('HP:0001166'), 'Arachnodactyly', True)
    focal_clonic_seizure_T = Phenotype(TermId.from_curie('HP:0002266'), 'Focal clonic seizure', True)
    seizure_T = Phenotype(TermId.from_curie('HP:0001250'), 'Seizure', True)
    spasticity_T = Phenotype(TermId.from_curie('HP:0001257'), 'Spasticity', True)
    arachnodactyly_F = Phenotype(TermId.from_curie('HP:0001166'), 'Arachnodactyly', False)
    focal_clonic_seizure_F = Phenotype(TermId.from_curie('HP:0002266'), 'Focal clonic seizure', False)
    seizure_F = Phenotype(TermId.from_curie('HP:0001250'), 'Seizure', False)
    spasticity_F = Phenotype(TermId.from_curie('HP:0001257'), 'Spasticity', False)


    prot_feat_1 = ProteinFeature.create(FeatureInfo('domain', Region(1, 75)), FeatureType.DOMAIN)
    prot_feat_2 = ProteinFeature.create(FeatureInfo('region', Region(50, 100)), FeatureType.REGION)
    prot = ProteinMetadata('NP_09876.5', 'FakeProtein', [prot_feat_1, prot_feat_2])

    snv = Variant.create_variant_from_scratch(VariantCoordinates(make_region(280, 281), 'A', 'G', 0), 'FakeGene',
                                                  'NM_1234.5', 'NM_1234.5:c.180A>G', False, [VariantEffect.MISSENSE_VARIANT], [1],
                                                  [prot], 60, 60,
                                                  Genotypes.from_mapping({
                                                      'A': Genotype.HETEROZYGOUS, 'B': Genotype.HETEROZYGOUS,
                                                      'C': Genotype.HOMOZYGOUS_ALTERNATE,
                                                      'D': Genotype.HETEROZYGOUS, 'E': Genotype.HETEROZYGOUS,
                                                      'G': Genotype.HETEROZYGOUS, 'J': Genotype.HETEROZYGOUS,
                                                      'K': Genotype.HETEROZYGOUS, 'M': Genotype.HETEROZYGOUS,
                                                      'N': Genotype.HOMOZYGOUS_ALTERNATE,
                                                      'P': Genotype.HETEROZYGOUS,
                                                      'Q': Genotype.HOMOZYGOUS_ALTERNATE,
                                                      'R': Genotype.HETEROZYGOUS, 'T': Genotype.HETEROZYGOUS,
                                                      'V': Genotype.HETEROZYGOUS, 'Y': Genotype.HETEROZYGOUS,
                                                  })
                                              )
    deletion = Variant.create_variant_from_scratch(VariantCoordinates(make_region(360, 363), 'TTC', 'T', -2),
                                                  'FakeGene', 'NM_1234.5', 'NM_1234.5:c.261_263del',
                                                  False, [VariantEffect.FRAMESHIFT_VARIANT],
                                                  [2], [prot], 86, 87,
                                                  Genotypes.from_mapping({
                                                      'D': Genotype.HETEROZYGOUS, 'F': Genotype.HETEROZYGOUS,
                                                      'G': Genotype.HETEROZYGOUS, 'H': Genotype.HETEROZYGOUS,
                                                      'I': Genotype.HOMOZYGOUS_ALTERNATE,
                                                      'L': Genotype.HETEROZYGOUS, 'O': Genotype.HETEROZYGOUS,
                                                      'R': Genotype.HETEROZYGOUS,
                                                      'S': Genotype.HOMOZYGOUS_ALTERNATE,
                                                      'U': Genotype.HETEROZYGOUS, 'W': Genotype.HETEROZYGOUS,
                                                      'X': Genotype.HETEROZYGOUS, 'Z': Genotype.HETEROZYGOUS,
                                                  })
                                                   )
    het_dup = Variant.create_variant_from_scratch(
        VariantCoordinates(make_region(175, 176), 'T', 'TG', 1), 'FakeGene', 'NM_1234.5',
        'NM_1234.5:c.75A>G', False, [VariantEffect.FRAMESHIFT_VARIANT], [1], [prot], 25, 25,
        Genotypes.empty())  # Not used in the patients below, hence `empty()`.
    hom_dup = Variant.create_variant_from_scratch(
        VariantCoordinates(make_region(175, 176), 'T', 'TG', 1),'FakeGene', 'NM_1234.5',
        'NM_1234.5:c.75A>G', False, [VariantEffect.FRAMESHIFT_VARIANT], [1], [prot], 25, 25,
        Genotypes.empty())  # Not used in the patients below, hence `empty()`.

    patients = (
        Patient(SampleLabels('A'),
                phenotypes=(arachnodactyly_T, spasticity_F, seizure_T),
                variants=[snv],
                proteins=[prot]
                ),
        Patient(SampleLabels('B'),
                phenotypes=(arachnodactyly_T, seizure_T, spasticity_T),
                variants=[snv],
                proteins=[prot]
                ),
        Patient(SampleLabels('C'),
                phenotypes=(arachnodactyly_F, spasticity_T, seizure_T),
                variants=[snv],
                proteins=[prot]
                ),
        Patient(SampleLabels('D'),
                phenotypes=(arachnodactyly_T, spasticity_T, seizure_T),
                variants=[snv, deletion],
                proteins=[prot]
                ),
        Patient(SampleLabels('E'),
                phenotypes=(arachnodactyly_T, spasticity_T, seizure_F),
                variants=[snv],
                proteins=[prot]
                ),
        Patient(SampleLabels('F'),
                phenotypes=(arachnodactyly_F, spasticity_F, seizure_T),
                variants=[deletion],
                proteins=[prot]
                ),
        Patient(SampleLabels('G'),
                phenotypes=(arachnodactyly_T, seizure_T, spasticity_T),
                variants=[snv, deletion],
                proteins=[prot]
                ),
        Patient(SampleLabels('H'),
                phenotypes=(arachnodactyly_T, seizure_T, spasticity_F),
                variants=[deletion],
                proteins=[prot]
                ),
        Patient(SampleLabels('I'),
                phenotypes=(arachnodactyly_F, spasticity_F, seizure_T),
                variants=[deletion],
                proteins=[prot]
                ),
        Patient(SampleLabels('J'),
                phenotypes=(arachnodactyly_T, seizure_T, spasticity_T),
                variants=[snv],
                proteins=[prot]
                ),
        Patient(SampleLabels('K'),
                phenotypes=(arachnodactyly_F, spasticity_T, seizure_T),
                variants=[snv],
                proteins=[prot]
                ),
        Patient(SampleLabels('L'),
                phenotypes=(arachnodactyly_F, seizure_F, spasticity_F),
                variants=[deletion],
                proteins=[prot]
                ),
        Patient(SampleLabels('M'),
                phenotypes=(arachnodactyly_T, seizure_F, spasticity_T),
                variants=[snv],
                proteins=[prot]
                ),
        Patient(SampleLabels('N'),
                phenotypes=(arachnodactyly_F, seizure_T, spasticity_F),
                variants=[snv],
                proteins=[prot]
                ),
        Patient(SampleLabels('O'),
                phenotypes=(arachnodactyly_F, seizure_F, spasticity_T),
                variants=[deletion],
                proteins=[prot]
                ),
        Patient(SampleLabels('P'),
                phenotypes=(arachnodactyly_T, seizure_T, spasticity_F),
                variants=[snv],
                proteins=[prot]
                ),
        Patient(SampleLabels('Q'),
                phenotypes=(arachnodactyly_T, seizure_F, spasticity_F),
                variants=[snv],
                proteins=[prot]
                ),
        Patient(SampleLabels('R'),
                phenotypes=(arachnodactyly_T, seizure_T, spasticity_F),
                variants=[snv, deletion],
                proteins=[prot]
                ),
        Patient(SampleLabels('S'),
                phenotypes=(arachnodactyly_F, seizure_T, spasticity_T),
                variants=[deletion],
                proteins=[prot]
                ),
        Patient(SampleLabels('T'),
                phenotypes=(arachnodactyly_T, seizure_F, spasticity_T),
                variants=[snv],
                proteins=[prot]
                ),
        Patient(SampleLabels('U'),
                phenotypes=(arachnodactyly_F, seizure_T, spasticity_T),
                variants=[deletion],
                proteins=[prot]
                ),
        Patient(SampleLabels('V'),
                phenotypes=(arachnodactyly_T, seizure_T, spasticity_T),
                variants=[snv],
                proteins=[prot]
                ),
        Patient(SampleLabels('W'),
                phenotypes=(arachnodactyly_F, seizure_T, spasticity_T),
                variants=[deletion],
                proteins=[prot]
                ),
        Patient(SampleLabels('X'),
                phenotypes=(arachnodactyly_F, seizure_T, spasticity_T),
                variants=[deletion],
                proteins=[prot]
                ),
        Patient(SampleLabels('Y'),
                phenotypes=(arachnodactyly_T, seizure_T, spasticity_T),
                variants=[snv],
                proteins=[prot]
                ),
        Patient(SampleLabels('Z'),
                phenotypes=(arachnodactyly_F, seizure_T, spasticity_T),
                variants=[deletion],
                proteins=[prot]
                ),
    )

    return Cohort.from_patients(patients)
