import pytest

import hpotk

from gpsea.model import (
    Cohort,
    SampleLabels,
    Patient,
    Sex,
    Phenotype,
    Variant,
    Genotype,
    Genotypes,
    VariantCoordinates,
    VariantEffect,
)
from gpsea.model.genome import GenomeBuild
from gpsea.analysis.predicate import VariantPredicate, true


@pytest.fixture(scope="package")
def always_false_variant_predicate() -> VariantPredicate:
    return ~true()


@pytest.fixture(scope="package")
def degenerated_cohort(
    genome_build: GenomeBuild,
) -> Cohort:
    """
    This cohort is "degenerated" because the members are annotated
    with unspecific HPO terms - just *Phenotypic abnormality*.
    """
    chr_22 = genome_build.contig_by_name("22")
    assert chr_22 is not None
    labels_a = SampleLabels("A")
    labels_b = SampleLabels("B")

    return Cohort.from_patients(
        members=(
            Patient.from_raw_parts(
                labels=labels_a,
                sex=Sex.UNKNOWN_SEX,
                age=None,
                vital_status=None,
                phenotypes=(
                    Phenotype.from_raw_parts(
                        term_id=hpotk.TermId.from_curie("HP:0000118"),
                        is_observed=True,
                    ),
                ),
                measurements=(),
                diseases=(),
                variants=(
                    Variant.create_variant_from_scratch(
                        VariantCoordinates.from_vcf_literal(
                            chr_22,
                            1001,
                            "C",
                            "G",
                        ),
                        "FAKE_GENE",
                        "NM_123.4",
                        hgvs_cdna="c.123C>G",
                        is_preferred=True,
                        consequences=(VariantEffect.MISSENSE_VARIANT,),
                        exons_effected=(2,),
                        protein_id=None,
                        hgvsp=None,
                        protein_effect_start=None,
                        protein_effect_end=None,
                        genotypes=Genotypes.single(
                            sample_id=labels_a,
                            genotype=Genotype.HETEROZYGOUS,
                        ),
                    ),
                ),
            ),
            Patient.from_raw_parts(
                labels=labels_b,
                sex=Sex.UNKNOWN_SEX,
                age=None,
                vital_status=None,
                phenotypes=(
                    Phenotype.from_raw_parts(
                        term_id=hpotk.TermId.from_curie("HP:0000118"),
                        is_observed=True,
                    ),
                ),
                measurements=(),
                diseases=(),
                variants=(
                    Variant.create_variant_from_scratch(
                        VariantCoordinates.from_vcf_literal(
                            chr_22,
                            1001,
                            "C",
                            "G",
                        ),
                        "FAKE_GENE",
                        "NM_123.4",
                        hgvs_cdna="c.123C>G",
                        is_preferred=True,
                        consequences=(VariantEffect.MISSENSE_VARIANT,),
                        exons_effected=(2,),
                        protein_id=None,
                        hgvsp=None,
                        protein_effect_start=None,
                        protein_effect_end=None,
                        genotypes=Genotypes.single(
                            sample_id=labels_b,
                            genotype=Genotype.HOMOZYGOUS_ALTERNATE,
                        ),
                    ),
                ),
            ),
        ),
    )
