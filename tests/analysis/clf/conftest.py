import pytest

from gpsea.model import (
    Genotype,
    Genotypes,
    Patient,
    SampleLabels,
    Sex,
    TranscriptAnnotation,
    Variant,
    VariantCoordinates,
    VariantEffect,
    VariantInfo,
)
from gpsea.model.genome import (
    GenomeBuild,
    GenomicRegion,
    Region,
    Strand,
)


"""
Genesis family - Autosomal dominant but can also be used as X dominant.

* Adam - father, unaffected
* Eve - mother, affected
* Cain - son, affected
"""


@pytest.fixture(scope="package")
def genesis_missense_mutation(
    genome_build: GenomeBuild,
    adam_label: SampleLabels,
    eve_label: SampleLabels,
    cain_label: SampleLabels,
) -> Variant:
    chr22 = genome_build.contig_by_name("chr22")
    assert chr22 is not None
    return Variant(
        variant_info=VariantInfo(
            variant_coordinates=VariantCoordinates(
                region=GenomicRegion(
                    contig=chr22,
                    start=100,
                    end=101,
                    strand=Strand.POSITIVE,
                ),
                ref="C",
                alt="G",
                change_length=0,
            )
        ),
        tx_annotations=(
            TranscriptAnnotation(
                gene_id="a_gene",
                tx_id="tx:xyz",
                hgvs_cdna=None,
                is_preferred=False,
                variant_effects=(
                    VariantEffect.MISSENSE_VARIANT,
                    VariantEffect.SPLICE_DONOR_VARIANT,
                ),
                affected_exons=(4,),
                protein_id="pt:xyz",
                hgvsp=None,
                protein_effect_coordinates=Region(40, 41),
            ),
        ),
        genotypes=Genotypes.from_mapping(
            {
                adam_label: Genotype.HOMOZYGOUS_REFERENCE,
                eve_label: Genotype.HETEROZYGOUS,
                cain_label: Genotype.HETEROZYGOUS,
            }
        ),
    )


@pytest.fixture(scope="package")
def genesis_synonymous_mutation(
    genome_build: GenomeBuild,
    adam_label: SampleLabels,
    eve_label: SampleLabels,
    cain_label: SampleLabels,
) -> Variant:
    chr22 = genome_build.contig_by_name("chr22")
    assert chr22 is not None
    return Variant(
        variant_info=VariantInfo(
            variant_coordinates=VariantCoordinates(
                region=GenomicRegion(
                    contig=chr22,
                    start=200,
                    end=201,
                    strand=Strand.POSITIVE,
                ),
                ref="T",
                alt="G",
                change_length=0,
            )
        ),
        tx_annotations=(
            TranscriptAnnotation(
                gene_id="a_gene",
                tx_id="tx:xyz",
                hgvs_cdna=None,
                is_preferred=True,
                variant_effects=(VariantEffect.SYNONYMOUS_VARIANT,),
                affected_exons=(5,),
                protein_id="pt:xyz",
                hgvsp=None,
                protein_effect_coordinates=Region(80, 81),
            ),
        ),
        genotypes=Genotypes.from_mapping(
            {
                adam_label: Genotype.HETEROZYGOUS,
                eve_label: Genotype.HOMOZYGOUS_REFERENCE,
                cain_label: Genotype.HOMOZYGOUS_REFERENCE,
            }
        ),
    )


@pytest.fixture(scope="package")
def adam_label() -> SampleLabels:
    return SampleLabels("Adam")


@pytest.fixture(scope="package")
def adam(
    adam_label: SampleLabels,
    genesis_missense_mutation: Variant,
    genesis_synonymous_mutation: Variant,
) -> Patient:
    return Patient.from_raw_parts(
        adam_label,
        sex=Sex.MALE,
        age=None,
        vital_status=None,
        phenotypes=(),
        measurements=(),
        diseases=(),
        variants=(
            genesis_missense_mutation,
            genesis_synonymous_mutation,
        ),
    )


@pytest.fixture(scope="package")
def eve_label() -> SampleLabels:
    return SampleLabels("Eve")


@pytest.fixture(scope="package")
def eve(
    eve_label: SampleLabels,
    genesis_missense_mutation: Variant,
    genesis_synonymous_mutation: Variant,
) -> Patient:
    return Patient.from_raw_parts(
        eve_label,
        sex=Sex.FEMALE,
        age=None,
        vital_status=None,
        phenotypes=(),
        measurements=(),
        diseases=(),
        variants=(
            genesis_missense_mutation,
            genesis_synonymous_mutation,
        ),
    )


@pytest.fixture(scope="package")
def cain_label() -> SampleLabels:
    return SampleLabels("Cain")


@pytest.fixture(scope="package")
def cain(
    cain_label: SampleLabels,
    genesis_missense_mutation: Variant,
    genesis_synonymous_mutation: Variant,
) -> Patient:
    return Patient.from_raw_parts(
        cain_label,
        sex=Sex.MALE,
        age=None,
        vital_status=None,
        phenotypes=(),
        measurements=(),
        diseases=(),
        variants=(
            genesis_missense_mutation,
            genesis_synonymous_mutation,
        ),
    )


"""
White family - Autosomal recessive

* Walt - father, HET
* Skyler - mother, HET
* Flynn - son, HOM_ALT
* Holly - daughter, HOM_REF
"""


@pytest.fixture(scope="package")
def white_missense_mutation(
    genome_build: GenomeBuild,
    walt_label: SampleLabels,
    skyler_label: SampleLabels,
    flynn_label: SampleLabels,
    holly_label: SampleLabels,
) -> Variant:
    chr22 = genome_build.contig_by_name("chr22")
    assert chr22 is not None
    return Variant(
        variant_info=VariantInfo(
            variant_coordinates=VariantCoordinates(
                region=GenomicRegion(
                    contig=chr22,
                    start=100,
                    end=101,
                    strand=Strand.POSITIVE,
                ),
                ref="C",
                alt="G",
                change_length=0,
            )
        ),
        tx_annotations=(
            TranscriptAnnotation(
                gene_id="a_gene",
                tx_id="tx:xyz",
                hgvs_cdna=None,
                is_preferred=True,
                variant_effects=(
                    VariantEffect.MISSENSE_VARIANT,
                    VariantEffect.SPLICE_DONOR_VARIANT,
                ),
                affected_exons=(4,),
                protein_id="pt:xyz",
                hgvsp=None,
                protein_effect_coordinates=Region(40, 41),
            ),
        ),
        genotypes=Genotypes.from_mapping(
            {
                walt_label: Genotype.HETEROZYGOUS,
                skyler_label: Genotype.HETEROZYGOUS,
                flynn_label: Genotype.HOMOZYGOUS_ALTERNATE,
                holly_label: Genotype.HOMOZYGOUS_REFERENCE,
            }
        ),
    )


@pytest.fixture(scope="package")
def white_synonymous_mutation(
    genome_build: GenomeBuild,
    walt_label: SampleLabels,
    skyler_label: SampleLabels,
    flynn_label: SampleLabels,
    holly_label: SampleLabels,
) -> Variant:
    chr22 = genome_build.contig_by_name("chr22")
    assert chr22 is not None
    return Variant(
        variant_info=VariantInfo(
            variant_coordinates=VariantCoordinates(
                region=GenomicRegion(
                    contig=chr22,
                    start=200,
                    end=201,
                    strand=Strand.POSITIVE,
                ),
                ref="T",
                alt="G",
                change_length=0,
            )
        ),
        tx_annotations=(
            TranscriptAnnotation(
                gene_id="a_gene",
                tx_id="tx:xyz",
                hgvs_cdna=None,
                is_preferred=True,
                variant_effects=(VariantEffect.SYNONYMOUS_VARIANT,),
                affected_exons=(5,),
                protein_id="pt:xyz",
                hgvsp=None,
                protein_effect_coordinates=Region(80, 81),
            ),
        ),
        genotypes=Genotypes.from_mapping(
            {
                walt_label: Genotype.HETEROZYGOUS,
                skyler_label: Genotype.HETEROZYGOUS,
                flynn_label: Genotype.HOMOZYGOUS_REFERENCE,
                holly_label: Genotype.HOMOZYGOUS_ALTERNATE,
            }
        ),
    )


@pytest.fixture(scope="package")
def walt_label() -> SampleLabels:
    return SampleLabels("Walt")


@pytest.fixture(scope="package")
def walt(
    walt_label: SampleLabels,
    white_missense_mutation: Variant,
    white_synonymous_mutation: Variant,
) -> Patient:
    return Patient.from_raw_parts(
        walt_label,
        sex=Sex.MALE,
        age=None,
        vital_status=None,
        phenotypes=(),
        measurements=(),
        diseases=(),
        variants=(
            white_missense_mutation,
            white_synonymous_mutation,
        ),
    )


@pytest.fixture(scope="package")
def skyler_label() -> SampleLabels:
    return SampleLabels("Skyler")


@pytest.fixture(scope="package")
def skyler(
    skyler_label: SampleLabels,
    white_missense_mutation: Variant,
    white_synonymous_mutation: Variant,
) -> Patient:
    return Patient.from_raw_parts(
        skyler_label,
        sex=Sex.FEMALE,
        age=None,
        vital_status=None,
        phenotypes=(),
        measurements=(),
        diseases=(),
        variants=(
            white_missense_mutation,
            white_synonymous_mutation,
        ),
    )


@pytest.fixture(scope="package")
def flynn_label() -> SampleLabels:
    return SampleLabels("Flynn")


@pytest.fixture(scope="package")
def flynn(
    flynn_label: SampleLabels,
    white_missense_mutation: Variant,
    white_synonymous_mutation: Variant,
) -> Patient:
    return Patient.from_raw_parts(
        flynn_label,
        sex=Sex.MALE,
        age=None,
        vital_status=None,
        phenotypes=(),
        measurements=(),
        diseases=(),
        variants=(
            white_missense_mutation,
            white_synonymous_mutation,
        ),
    )


@pytest.fixture(scope="package")
def holly_label() -> SampleLabels:
    return SampleLabels("Holly")


@pytest.fixture(scope="package")
def holly(
    holly_label: SampleLabels,
    white_missense_mutation: Variant,
    white_synonymous_mutation: Variant,
) -> Patient:
    return Patient.from_raw_parts(
        holly_label,
        sex=Sex.FEMALE,
        age=None,
        vital_status=None,
        phenotypes=(),
        measurements=(),
        diseases=(),
        variants=(
            white_missense_mutation,
            white_synonymous_mutation,
        ),
    )


"""
Skywalker family - X-linked recessive

* Anakin - father, homozygous reference  (possibly hemizygous reference?)
* Padme - mother, heterozygous
* Luke - son, hemizygous
* Leia - daughter, heterozygous
"""


@pytest.fixture(scope="package")
def skywalker_mutation(
    genome_build: GenomeBuild,
    anakin_label: SampleLabels,
    padme_label: SampleLabels,
    luke_label: SampleLabels,
    leia_label: SampleLabels,
) -> Variant:
    chrX = genome_build.contig_by_name("chrX")
    assert chrX is not None
    return Variant(
        variant_info=VariantInfo(
            variant_coordinates=VariantCoordinates(
                region=GenomicRegion(
                    contig=chrX,
                    start=100,
                    end=101,
                    strand=Strand.POSITIVE,
                ),
                ref="C",
                alt="G",
                change_length=0,
            )
        ),
        tx_annotations=(
            TranscriptAnnotation(
                gene_id="a_gene",
                tx_id="tx:xyz",
                hgvs_cdna=None,
                is_preferred=False,
                variant_effects=(
                    VariantEffect.MISSENSE_VARIANT,
                    VariantEffect.SPLICE_DONOR_VARIANT,
                ),
                affected_exons=(4,),
                protein_id="pt:xyz",
                hgvsp=None,
                protein_effect_coordinates=Region(40, 41),
            ),
        ),
        genotypes=Genotypes.from_mapping(
            {
                anakin_label: Genotype.HOMOZYGOUS_REFERENCE,
                padme_label: Genotype.HETEROZYGOUS,
                luke_label: Genotype.HEMIZYGOUS,
                leia_label: Genotype.HETEROZYGOUS,
            }
        ),
    )


@pytest.fixture(scope="package")
def anakin_label() -> SampleLabels:
    return SampleLabels("Anakin")


@pytest.fixture(scope="package")
def anakin(
    anakin_label: SampleLabels,
    skywalker_mutation: Variant,
) -> Patient:
    return Patient.from_raw_parts(
        anakin_label,
        sex=Sex.MALE,
        age=None,
        vital_status=None,
        phenotypes=(),
        measurements=(),
        diseases=(),
        variants=(skywalker_mutation,),
    )


@pytest.fixture(scope="package")
def padme_label() -> SampleLabels:
    return SampleLabels("Padme")


@pytest.fixture(scope="package")
def padme(
    padme_label: SampleLabels,
    skywalker_mutation: Variant,
) -> Patient:
    return Patient.from_raw_parts(
        padme_label,
        sex=Sex.FEMALE,
        age=None,
        vital_status=None,
        phenotypes=(),
        measurements=(),
        diseases=(),
        variants=(skywalker_mutation,),
    )


@pytest.fixture(scope="package")
def luke_label() -> SampleLabels:
    return SampleLabels("Luke")


@pytest.fixture(scope="package")
def luke(
    luke_label: SampleLabels,
    skywalker_mutation: Variant,
) -> Patient:
    return Patient.from_raw_parts(
        luke_label,
        sex=Sex.MALE,
        age=None,
        vital_status=None,
        phenotypes=(),
        measurements=(),
        diseases=(),
        variants=(skywalker_mutation,),
    )


@pytest.fixture(scope="package")
def leia_label() -> SampleLabels:
    return SampleLabels("Leia")


@pytest.fixture(scope="package")
def leia(
    leia_label: SampleLabels,
    skywalker_mutation: Variant,
) -> Patient:
    return Patient.from_raw_parts(
        leia_label,
        sex=Sex.FEMALE,
        age=None,
        vital_status=None,
        phenotypes=(),
        measurements=(),
        diseases=(),
        variants=(skywalker_mutation,),
    )
