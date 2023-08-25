from hpotk import TermId

from genophenocorr.variant import Variant
from genophenocorr.phenotype import Phenotype
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

    arachnodactyly = Phenotype(TermId.from_curie('HP:0001166'), 'Arachnodactyly', True)
    focal_clonic_seizure = Phenotype(TermId.from_curie('HP:0002266'), 'Focal clonic seizure', True)
    seizure = Phenotype(TermId.from_curie('HP:0001250'), 'Seizure', True)

    # TODO - make variants
    hom = Variant('hom', '', None, (), 'homozygous')
    het = Variant('het', '', None, (), 'heterozygous')

    patients = (
        Patient('A',
                phenotypes=(arachnodactyly, focal_clonic_seizure),
                variants=(hom,),
                proteins=()
                ),
        Patient('B',
                phenotypes=(arachnodactyly, seizure),
                variants=(het,),
                proteins=()
                )
    )

    return Cohort.from_patients(patients)
