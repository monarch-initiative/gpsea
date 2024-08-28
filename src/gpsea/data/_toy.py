from hpotk import TermId

from gpsea.model import *
from gpsea.model.genome import Contig, GenomicRegion, Strand

CONTIG = Contig("1", "GB_ACC", "REFSEQ_NAME", "UCSC_NAME", 1_000)


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

    arachnodactyly_T = Phenotype(TermId.from_curie("HP:0001166"), True)
    focal_clonic_seizure_T = Phenotype(TermId.from_curie("HP:0002266"), True)
    seizure_T = Phenotype(TermId.from_curie("HP:0001250"), True)
    spasticity_T = Phenotype(TermId.from_curie("HP:0001257"), True)
    arachnodactyly_F = Phenotype(TermId.from_curie("HP:0001166"), False)
    focal_clonic_seizure_F = Phenotype(TermId.from_curie("HP:0002266"), False)
    seizure_F = Phenotype(TermId.from_curie("HP:0001250"), False)
    spasticity_F = Phenotype(TermId.from_curie("HP:0001257"), False)

    Disease_T = Disease(TermId.from_curie("OMIM:001234"), "Test present disease", True)
    Disease_F = Disease(TermId.from_curie("OMIM:001234"), "Test absent disease", False)

    snv = Variant.create_variant_from_scratch(
        VariantCoordinates(make_region(280, 281), "A", "G", 0),
        "FakeGene",
        "NM_1234.5",
        "NM_1234.5:c.180A>G",
        False,
        [VariantEffect.MISSENSE_VARIANT],
        [1],
        "NP_09876.5",
        "NP_09876.5:p.Unk200",
        60,
        60,
        Genotypes.from_mapping(
            {
                SampleLabels("A"): Genotype.HETEROZYGOUS,
                SampleLabels("B"): Genotype.HETEROZYGOUS,
                SampleLabels("C"): Genotype.HOMOZYGOUS_ALTERNATE,
                SampleLabels("D"): Genotype.HETEROZYGOUS,
                SampleLabels("E"): Genotype.HETEROZYGOUS,
                SampleLabels("G"): Genotype.HETEROZYGOUS,
                SampleLabels("J"): Genotype.HETEROZYGOUS,
                SampleLabels("K"): Genotype.HETEROZYGOUS,
                SampleLabels("M"): Genotype.HETEROZYGOUS,
                SampleLabels("N"): Genotype.HOMOZYGOUS_ALTERNATE,
                SampleLabels("P"): Genotype.HETEROZYGOUS,
                SampleLabels("Q"): Genotype.HOMOZYGOUS_ALTERNATE,
                SampleLabels("R"): Genotype.HETEROZYGOUS,
                SampleLabels("T"): Genotype.HETEROZYGOUS,
                SampleLabels("V"): Genotype.HETEROZYGOUS,
                SampleLabels("Y"): Genotype.HETEROZYGOUS,
            }
        ),
    )
    deletion = Variant.create_variant_from_scratch(
        VariantCoordinates(make_region(360, 363), "TTC", "T", -2),
        "FakeGene",
        "NM_1234.5",
        "NM_1234.5:c.261_263del",
        False,
        [VariantEffect.FRAMESHIFT_VARIANT],
        [2],
        "NP_09876.5",
        "NP_09876.5:p.Unk200",
        86,
        87,
        Genotypes.from_mapping(
            {
                SampleLabels("D"): Genotype.HETEROZYGOUS,
                SampleLabels("F"): Genotype.HETEROZYGOUS,
                SampleLabels("G"): Genotype.HETEROZYGOUS,
                SampleLabels("H"): Genotype.HETEROZYGOUS,
                SampleLabels("I"): Genotype.HOMOZYGOUS_ALTERNATE,
                SampleLabels("L"): Genotype.HETEROZYGOUS,
                SampleLabels("O"): Genotype.HETEROZYGOUS,
                SampleLabels("R"): Genotype.HETEROZYGOUS,
                SampleLabels("S"): Genotype.HOMOZYGOUS_ALTERNATE,
                SampleLabels("U"): Genotype.HETEROZYGOUS,
                SampleLabels("W"): Genotype.HETEROZYGOUS,
                SampleLabels("X"): Genotype.HETEROZYGOUS,
                SampleLabels("Z"): Genotype.HETEROZYGOUS,
            }
        ),
    )
    het_dup = Variant.create_variant_from_scratch(
        VariantCoordinates(make_region(175, 176), "T", "TG", 1),
        "FakeGene",
        "NM_1234.5",
        "NM_1234.5:c.75A>G",
        False,
        [VariantEffect.FRAMESHIFT_VARIANT],
        [1],
        "NP_09876.5",
        "NP_09876.5:p.Unk200",
        25,
        25,
        Genotypes.empty(),
    )  # Not used in the patients below, hence `empty()`.
    hom_dup = Variant.create_variant_from_scratch(
        VariantCoordinates(make_region(175, 176), "T", "TG", 1),
        "FakeGene",
        "NM_1234.5",
        "NM_1234.5:c.75A>G",
        False,
        [VariantEffect.FRAMESHIFT_VARIANT],
        [1],
        "NP_09876.5",
        "NP_09876.5:p.Unk200",
        25,
        25,
        Genotypes.empty(),
    )  # Not used in the patients below, hence `empty()`.

    patients = (
        Patient.from_raw_parts(
            SampleLabels("A"),
            sex=None,
            phenotypes=(arachnodactyly_T, spasticity_F, seizure_T),
            variants=[snv],
            diseases=[Disease_T],
        ),
        Patient.from_raw_parts(
            SampleLabels("B"),
            sex=None,
            phenotypes=(arachnodactyly_T, seizure_T, spasticity_T),
            variants=[snv],
            diseases=[Disease_T],
        ),
        Patient.from_raw_parts(
            SampleLabels("C"),
            sex=None,
            phenotypes=(arachnodactyly_F, spasticity_T, seizure_T),
            variants=[snv],
            diseases=[Disease_T],
        ),
        Patient.from_raw_parts(
            SampleLabels("D"),
            sex=None,
            phenotypes=(arachnodactyly_T, spasticity_T, seizure_T),
            variants=[snv, deletion],
            diseases=[Disease_T],
        ),
        Patient.from_raw_parts(
            SampleLabels("E"),
            sex=None,
            phenotypes=(arachnodactyly_T, spasticity_T, seizure_F),
            variants=[snv],
            diseases=[Disease_T],
        ),
        Patient.from_raw_parts(
            SampleLabels("F"),
            sex=None,
            phenotypes=(arachnodactyly_F, spasticity_F, seizure_T),
            variants=[deletion],
            diseases=[Disease_T],
        ),
        Patient.from_raw_parts(
            SampleLabels("G"),
            sex=None,
            phenotypes=(arachnodactyly_T, seizure_T, spasticity_T),
            variants=[snv, deletion],
            diseases=[Disease_T],
        ),
        Patient.from_raw_parts(
            SampleLabels("H"),
            sex=None,
            phenotypes=(arachnodactyly_T, seizure_T, spasticity_F),
            variants=[deletion],
            diseases=[Disease_T],
        ),
        Patient.from_raw_parts(
            SampleLabels("I"),
            sex=None,
            phenotypes=(arachnodactyly_F, spasticity_F, seizure_T),
            variants=[deletion],
            diseases=[Disease_T],
        ),
        Patient.from_raw_parts(
            SampleLabels("J"),
            sex=None,
            phenotypes=(arachnodactyly_T, seizure_T, spasticity_T),
            variants=[snv],
            diseases=[Disease_T],
        ),
        Patient.from_raw_parts(
            SampleLabels("K"),
            sex=None,
            phenotypes=(arachnodactyly_F, spasticity_T, seizure_T),
            variants=[snv],
            diseases=[Disease_T],
        ),
        Patient.from_raw_parts(
            SampleLabels("L"),
            sex=None,
            phenotypes=(arachnodactyly_F, seizure_F, spasticity_F),
            variants=[deletion],
            diseases=[Disease_T],
        ),
        Patient.from_raw_parts(
            SampleLabels("M"),
            sex=None,
            phenotypes=(arachnodactyly_T, seizure_F, spasticity_T),
            variants=[snv],
            diseases=[Disease_T],
        ),
        Patient.from_raw_parts(
            SampleLabels("N"),
            sex=None,
            phenotypes=(arachnodactyly_F, seizure_T, spasticity_F),
            variants=[snv],
            diseases=[Disease_T],
        ),
        Patient.from_raw_parts(
            SampleLabels("O"),
            sex=None,
            phenotypes=(arachnodactyly_F, seizure_F, spasticity_T),
            variants=[deletion],
            diseases=[Disease_T],
        ),
        Patient.from_raw_parts(
            SampleLabels("P"),
            sex=None,
            phenotypes=(arachnodactyly_T, seizure_T, spasticity_F),
            variants=[snv],
            diseases=[Disease_T],
        ),
        Patient.from_raw_parts(
            SampleLabels("Q"),
            sex=None,
            phenotypes=(arachnodactyly_T, seizure_F, spasticity_F),
            variants=[snv],
            diseases=[Disease_T],
        ),
        Patient.from_raw_parts(
            SampleLabels("R"),
            sex=None,
            phenotypes=(arachnodactyly_T, seizure_T, spasticity_F),
            variants=[snv, deletion],
            diseases=[Disease_T],
        ),
        Patient.from_raw_parts(
            SampleLabels("S"),
            sex=None,
            phenotypes=(arachnodactyly_F, seizure_T, spasticity_T),
            variants=[deletion],
            diseases=[Disease_T],
        ),
        Patient.from_raw_parts(
            SampleLabels("T"),
            sex=None,
            phenotypes=(arachnodactyly_T, seizure_F, spasticity_T),
            variants=[snv],
            diseases=[Disease_T],
        ),
        Patient.from_raw_parts(
            SampleLabels("U"),
            sex=None,
            phenotypes=(arachnodactyly_F, seizure_T, spasticity_T),
            variants=[deletion],
            diseases=[Disease_T],
        ),
        Patient.from_raw_parts(
            SampleLabels("V"),
            sex=None,
            phenotypes=(arachnodactyly_T, seizure_T, spasticity_T),
            variants=[snv],
            diseases=[Disease_T],
        ),
        Patient.from_raw_parts(
            SampleLabels("W"),
            sex=None,
            phenotypes=(arachnodactyly_F, seizure_T, spasticity_T),
            variants=[deletion],
            diseases=[Disease_T],
        ),
        Patient.from_raw_parts(
            SampleLabels("X"),
            sex=None,
            phenotypes=(arachnodactyly_F, seizure_T, spasticity_T),
            variants=[deletion],
            diseases=[Disease_T],
        ),
        Patient.from_raw_parts(
            SampleLabels("Y"),
            sex=None,
            phenotypes=(arachnodactyly_T, seizure_T, spasticity_T),
            variants=[snv],
            diseases=[Disease_T],
        ),
        Patient.from_raw_parts(
            SampleLabels("Z"),
            sex=None,
            phenotypes=(arachnodactyly_F, seizure_T, spasticity_T),
            variants=[deletion],
            diseases=[Disease_T],
        ),
    )

    return Cohort.from_patients(patients)
