from hpotk import TermId

from genophenocorr.variant import Variant, BasicVariantCoordinateFinder
from genophenocorr.phenotype import Phenotype
from genophenocorr.protein import ProteinMetadata, SimpleProteinFeature, FeatureType, FeatureInfo
from genophenocorr.patient import Patient
from genophenocorr.cohort import Cohort


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


    # TODO - make variants
    basic_var = BasicVariantCoordinateFinder()
    prot_feat_1 = SimpleProteinFeature(FeatureInfo('domain', 1, 75), FeatureType.DOMAIN)
    prot_feat_2 = SimpleProteinFeature(FeatureInfo('region', 50, 100), FeatureType.REGION)
    prot = ProteinMetadata('NP_09876.5', 'FakeProtein', [prot_feat_1, prot_feat_2])
    het_snv = Variant.create_variant_from_scratch('HetVar1', 'SNV', basic_var.find_coordinates('chr1', '280', '281', 'A', 'G', 'heterozygous'), 'FakeGene', 'NM_1234.5', 'NM_1234.5:c.180A>G', ['missense_variant'], [1], [prot], 60, 60)
    het_del = Variant.create_variant_from_scratch('HetVar2', 'indel', basic_var.find_coordinates('chr1', '360', '363', 'TTC', 'T', 'heterozygous'), 'FakeGene', 'NM_1234.5', 'NM_1234.5:c.261_263del', ['frameshift_variant'], [2], [prot], 86, 87)
    het_dup = Variant.create_variant_from_scratch('HetVar3', 'insertion', basic_var.find_coordinates('chr1', '175', '177', 'T', 'TG', 'heterozygous'), 'FakeGene', 'NM_1234.5', 'NM_1234.5:c.75A>G', ['frameshift_variant'], [1], [prot], 25, 25)
    hom_snv = Variant.create_variant_from_scratch('HomVar1', 'SNV', basic_var.find_coordinates('chr1', '280', '281', 'A', 'G', 'homozygous'), 'FakeGene', 'NM_1234.5', 'NM_1234.5:c.180A>G', ['missense_variant'], [1], [prot], 60, 60)
    hom_del = Variant.create_variant_from_scratch('HomVar2', 'indel', basic_var.find_coordinates('chr1', '360', '363', 'TTC', 'T', 'homozygous'), 'FakeGene', 'NM_1234.5', 'NM_1234.5:c.261_263del', ['frameshift_variant'], [2], [prot], 86, 87)
    hom_dup = Variant.create_variant_from_scratch('HomVar3', 'insertion', basic_var.find_coordinates('chr1', '175', '177', 'T', 'TG', 'homozygous'), 'FakeGene', 'NM_1234.5', 'NM_1234.5:c.75A>G', ['frameshift_variant'], [1], [prot], 25, 25)

    patients = (
        Patient('A',
                phenotypes=(arachnodactyly_T, spasticity_F, seizure_T),
                variants=[het_snv],
                proteins=[prot]
                ),
        Patient('B',
                phenotypes=(arachnodactyly_T, seizure_T, spasticity_T),
                variants=[het_snv],
                proteins=[prot]
                ),
        Patient('C',
                phenotypes=(arachnodactyly_F, spasticity_T, seizure_T),
                variants=[hom_snv],
                proteins=[prot]
                ),
        Patient('D',
                phenotypes=(arachnodactyly_T, spasticity_T, seizure_T),
                variants=[het_snv, het_del],
                proteins=[prot]
                ),
        Patient('E',
                phenotypes=(arachnodactyly_T, spasticity_T, seizure_F),
                variants=[het_snv],
                proteins=[prot]
                ),
        Patient('F',
                phenotypes=(arachnodactyly_F, spasticity_F, seizure_T),
                variants=[het_del],
                proteins=[prot]
                ),
        Patient('G',
                phenotypes=(arachnodactyly_T, seizure_T, spasticity_T),
                variants=[het_snv, het_del],
                proteins=[prot]
                ),
        Patient('H',
                phenotypes=(arachnodactyly_T, seizure_T, spasticity_F),
                variants=[het_del],
                proteins=[prot]
                ),
        Patient('I',
                phenotypes=(arachnodactyly_F, spasticity_F, seizure_T),
                variants=[hom_del],
                proteins=[prot]
                ),
        Patient('J',
                phenotypes=(arachnodactyly_T, seizure_T, spasticity_T),
                variants=[het_snv],
                proteins=[prot]
                ),
        Patient('K',
                phenotypes=(arachnodactyly_F, spasticity_T, seizure_T),
                variants=[het_snv],
                proteins=[prot]
                ),
        Patient('L',
                phenotypes=(arachnodactyly_F, seizure_F, spasticity_F),
                variants=[het_del],
                proteins=[prot]
                ),
        Patient('M',
                phenotypes=(arachnodactyly_T, seizure_F, spasticity_T),
                variants=[het_snv],
                proteins=[prot]
                ),
        Patient('N',
                phenotypes=(arachnodactyly_F, seizure_T, spasticity_F),
                variants=[hom_snv],
                proteins=[prot]
                ),
        Patient('O',
                phenotypes=(arachnodactyly_F, seizure_F, spasticity_T),
                variants=[het_del],
                proteins=[prot]
                ),
        Patient('P',
                phenotypes=(arachnodactyly_T, seizure_T, spasticity_F),
                variants=[het_snv],
                proteins=[prot]
                ),
        Patient('Q',
                phenotypes=(arachnodactyly_T, seizure_F, spasticity_F),
                variants=[hom_snv],
                proteins=[prot]
                ),
        Patient('R',
                phenotypes=(arachnodactyly_T, seizure_T, spasticity_F),
                variants=[het_snv, het_del],
                proteins=[prot]
                ),
        Patient('S',
                phenotypes=(arachnodactyly_F, seizure_T, spasticity_T),
                variants=[hom_del],
                proteins=[prot]
                ),
        Patient('T',
                phenotypes=(arachnodactyly_T, seizure_F, spasticity_T),
                variants=[het_snv],
                proteins=[prot]
                ),
        Patient('U',
                phenotypes=(arachnodactyly_F, seizure_T, spasticity_T),
                variants=[het_del],
                proteins=[prot]
                ),
        Patient('V',
                phenotypes=(arachnodactyly_T, seizure_T, spasticity_T),
                variants=[het_snv],
                proteins=[prot]
                ),        
        Patient('W',
                phenotypes=(arachnodactyly_F, seizure_T, spasticity_T),
                variants=[het_del],
                proteins=[prot]
                ),
        Patient('X',
                phenotypes=(arachnodactyly_F, seizure_T, spasticity_T),
                variants=[het_del],
                proteins=[prot]
                ),
        Patient('Y',
                phenotypes=(arachnodactyly_T, seizure_T, spasticity_T),
                variants=[het_snv],
                proteins=[prot]
                ),
        Patient('Z',
                phenotypes=(arachnodactyly_F, seizure_T, spasticity_T),
                variants=[het_del],
                proteins=[prot]
                ),
    )

    return Cohort.from_patients(patients)
