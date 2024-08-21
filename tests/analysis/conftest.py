import pytest

import hpotk

from genophenocorr.model import *
from genophenocorr.model.genome import *
from genophenocorr.analysis.predicate.genotype import VariantPredicate


class AlwaysFalseVariantPredicate(VariantPredicate):
    def get_question(self) -> str:
        return "No question asked, just always returns False"

    def test(self, variant: Variant) -> bool:
        return False


@pytest.fixture(scope="package")
def always_false_variant_predicate() -> VariantPredicate:
    return AlwaysFalseVariantPredicate()


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
            Patient(
                labels=labels_a,
                phenotypes=(
                    Phenotype(
                        term_id=hpotk.TermId.from_curie("HP:0000118"),
                        is_observed=True,
                    ),
                ),
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
            Patient(
                labels=labels_b,
                phenotypes=(
                    Phenotype(
                        term_id=hpotk.TermId.from_curie("HP:0000118"),
                        is_observed=True,
                    ),
                ),
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
