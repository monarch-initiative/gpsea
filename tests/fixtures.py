import hpotk
import pytest

from genophenocorr.model import *
from genophenocorr.model.genome import GRCh38, GenomicRegion, Region, Strand


def make_region(contig: str, start: int, end: int) -> GenomicRegion:
    return GenomicRegion(GRCh38.contig_by_name(contig), start, end, Strand.POSITIVE)


@pytest.fixture
def toy_cohort() -> Cohort:

    prot = ProteinMetadata(protein_id='NP_037407.4', label='Ankyrin repeat domain-containing protein 11',
            protein_features=(ProteinFeature.create(feature_type=FeatureType.REPEAT, info=FeatureInfo('ANK 1', Region(start=167, end=196))),
            ProteinFeature.create(feature_type=FeatureType.REPEAT, info=FeatureInfo('ANK 2', Region(start=200, end=229))),
            ProteinFeature.create(feature_type=FeatureType.REPEAT, info=FeatureInfo('ANK 3', Region(start=233, end=262))),
            ProteinFeature.create(feature_type=FeatureType.REPEAT, info=FeatureInfo('ANK 4', Region(start=266, end=292))),
            ProteinFeature.create(feature_type=FeatureType.REGION, info=FeatureInfo('Disordered', Region(start=1, end=90))),
            ProteinFeature.create(feature_type=FeatureType.REGION, info=FeatureInfo('Disordered', Region(start=128, end=169))),
            ProteinFeature.create(feature_type=FeatureType.REGION, info=FeatureInfo('Disordered', Region(start=289, end=380))),
            ProteinFeature.create(feature_type=FeatureType.REGION, info=FeatureInfo('Disordered', Region(start=398, end=647))),
            ProteinFeature.create(feature_type=FeatureType.REGION, info=FeatureInfo('Disordered', Region(start=723, end=783))),
            ProteinFeature.create(feature_type=FeatureType.REGION, info=FeatureInfo('Disordered', Region(start=881, end=1043))),
            ProteinFeature.create(feature_type=FeatureType.REGION, info=FeatureInfo('Disordered', Region(start=1059, end=1393))),
            ProteinFeature.create(feature_type=FeatureType.REGION, info=FeatureInfo('Disordered', Region(start=1424, end=1710))),
            ProteinFeature.create(feature_type=FeatureType.REGION, info=FeatureInfo('Disordered', Region(start=1814, end=1836))),
            ProteinFeature.create(feature_type=FeatureType.REGION, info=FeatureInfo('Disordered', Region(start=1988, end=2019))),
            ProteinFeature.create(feature_type=FeatureType.REGION, info=FeatureInfo('Disordered', Region(start=2131, end=2406))),
            ProteinFeature.create(feature_type=FeatureType.REGION, info=FeatureInfo('Important for protein degradation', Region(start=2369, end=2663)))))

    phenos = get_test_phenotypes()

    dup = Variant(VariantCoordinates(make_region("16", 89279849, 89279850), ref='G', alt='GC', change_length=1),
                  [
                      TranscriptAnnotation('ANKRD11', 'NM_013275.6', 'NM_013275.6:c.6691dup', False, [VariantEffect.FRAMESHIFT_VARIANT], [9],
                                           [prot], Region(2230, 2231))
                  ],
                  Genotypes.from_mapping({SampleLabels('HetSingleVar'): Genotype.HETEROZYGOUS}))
    indel = Variant(VariantCoordinates(make_region("16", 89284600, 89284602), ref='GG', alt='A', change_length=-1),
                    [
                        TranscriptAnnotation('ANKRD11', 'NM_013275.6', 'NM_013275.6:c.1940_1941delinsT', False, [VariantEffect.FRAMESHIFT_VARIANT],
                                             [9], [prot], Region(646, 647))
                    ],
                    Genotypes.from_mapping({SampleLabels('HetDoubleVar1'): Genotype.HETEROZYGOUS}))
    snv_stop_gain = Variant(VariantCoordinates(make_region("16", 89280751, 89280752), ref='G', alt='T', change_length=0),
                            [
                                TranscriptAnnotation('ANKRD11', 'NM_013275.6', 'NM_013275.6:c.5790C>A', False, [VariantEffect.STOP_GAINED], [9], [prot],
                             Region(1929, 1930))],
                            Genotypes.from_mapping({SampleLabels('HetDoubleVar1'): Genotype.HETEROZYGOUS}))
    snv_missense = Variant(VariantCoordinates(make_region("16", 89275127, 89275128), ref='G', alt='A', change_length=0),
                           [
                               TranscriptAnnotation('ANKRD11', 'NM_013275.6', 'NM_013275.6:c.7534C>T', False, [VariantEffect.MISSENSE_VARIANT], [10],
                             [prot], Region(2511, 2512))
                           ],
                           Genotypes.from_mapping({SampleLabels('HetDoubleVar2'): Genotype.HETEROZYGOUS}))
    del_frameshift = Variant(VariantCoordinates(make_region("16", 89279707, 89279725), ref='AGTGTTCGGGGCGGGGCC', alt='A', change_length=-17),
                             [
                                 TranscriptAnnotation('ANKRD11', 'NM_013275.6', 'NM_013275.6:c.6817_6833del', False, [VariantEffect.FRAMESHIFT_VARIANT],
                              [9], [prot], Region(2272, 2278))
                             ],
                             Genotypes.from_mapping({SampleLabels('HetDoubleVar2'): Genotype.HETEROZYGOUS}))
    del_small = Variant(VariantCoordinates(make_region("16", 89279457, 89279459), ref='TG', alt='T', change_length=-1),
                        [
                            TranscriptAnnotation('ANKRD11', 'NM_013275.6', 'NM_013275.6:c.7083del', False, [VariantEffect.FRAMESHIFT_VARIANT], [9],
                             [prot], Region(2360, 2362))
                        ],
                        Genotypes.from_mapping({SampleLabels('HomoVar'): Genotype.HOMOZYGOUS_ALTERNATE}))
    del_large = Variant(VariantCoordinates(make_region("16", 89_190_070, 89_439_815), ref='N', alt='<DEL>', change_length=-249_745),
                        [
                            TranscriptAnnotation('ANKRD11', 'NM_013275.6', None, False,
                                 [VariantEffect.STOP_LOST, VariantEffect.FEATURE_TRUNCATION, VariantEffect.CODING_SEQUENCE_VARIANT, VariantEffect.FIVE_PRIME_UTR_VARIANT,
                                  VariantEffect.THREE_PRIME_UTR_VARIANT, VariantEffect.INTRON_VARIANT], [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13],
                                 [prot], None)
                        ],
                        Genotypes.from_mapping({SampleLabels('LargeCNV'): Genotype.HETEROZYGOUS}))

    patients = (
        Patient(SampleLabels('HetSingleVar'),
                phenotypes=(phenos['arachnodactyly_T'], phenos['spasticity_F'], phenos['focal_clonic_seizure_T']),
                variants=(dup,),
                proteins=[prot]
                ),
        Patient(SampleLabels('HetDoubleVar1'),
                phenotypes=(phenos['arachnodactyly_T'], phenos['seizure_T'], phenos['spasticity_T']),
                variants=(indel, snv_stop_gain),
                proteins=[prot]
                ),
        Patient(SampleLabels('HetDoubleVar2'),
                phenotypes=(phenos['arachnodactyly_F'], phenos['spasticity_T'], phenos['seizure_T']),
                variants=(snv_missense, del_frameshift),
                proteins=[prot]
                ),
        Patient(SampleLabels('HomoVar'),
                phenotypes=(phenos['arachnodactyly_T'], phenos['spasticity_T'], phenos['seizure_T']),
                variants=(del_small,),
                proteins=[prot]
                ),
        Patient(SampleLabels('LargeCNV'),
                phenotypes=(phenos['arachnodactyly_T'], phenos['spasticity_T'], phenos['seizure_F']),
                variants=(del_large,),
                proteins=[prot]
                ),
    )

    return Cohort.from_patients(patients)


def get_test_phenotypes():
    phenotypes = {}

    phenotypes['arachnodactyly_T'] = Phenotype(hpotk.TermId.from_curie('HP:0001166'), "Arachnodactyly", True)
    phenotypes['seizure_T'] = Phenotype(hpotk.TermId.from_curie('HP:0001250'), "Seizure", True)
    phenotypes['focal_clonic_seizure_T'] = Phenotype(hpotk.TermId.from_curie('HP:0002266'), "Focal clonic seizure", True)
    phenotypes['spasticity_T'] = Phenotype(hpotk.TermId.from_curie('HP:0001257'), "Spasticity", True)
    phenotypes['arachnodactyly_F'] = Phenotype(hpotk.TermId.from_curie('HP:0001166'), "Arachnodactyly", False)
    phenotypes['seizure_F'] = Phenotype(hpotk.TermId.from_curie('HP:0001250'), "Seizure", False)
    phenotypes['spasticity_F'] = Phenotype(hpotk.TermId.from_curie('HP:0001257'), "Spasticity", False)
    phenotypes['focal_clonic_seizure_F'] = Phenotype(hpotk.TermId.from_curie('HP:0002266'), "Focal clonic seizure", False)

    return phenotypes

