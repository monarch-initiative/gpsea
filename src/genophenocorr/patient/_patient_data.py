import typing

from genophenocorr.phenotype import Phenotype
from genophenocorr.variant import Variant
from genophenocorr.protein import ProteinMetadata


class Patient:
    def __init__(self, patient_id: str,
                 phenotypes: typing.Sequence[Phenotype],
                 variants: typing.Sequence[Variant],
                 proteins: typing.Sequence[ProteinMetadata]):
        self._id = patient_id
        # TODO(lnrekerle) - We should wrap the sequences below in a tuple, otherwise the `__hash__` will explode if the
        #  caller provides a list.
        self._phenotypes = phenotypes
        self._variants = variants
        self._proteins = proteins

    @property
    def patient_id(self) -> str:
        return self._id

    @property
    def phenotypes(self) -> typing.Sequence[Phenotype]:
        return self._phenotypes
    
    @property
    def variants(self) -> typing.Sequence[Variant]:
        return self._variants

    @property
    def proteins(self) -> typing.Sequence[ProteinMetadata]:
        return self._proteins


    # def has_hpo(self, hpo, all_hpo):
    #     if not isinstance(all_hpo, defaultdict):
    #         for h in self.phenotype_ids:
    #             if h == hpo:
    #                 return True
    #         return False
    #     else:
    #         for h in self.phenotype_ids:
    #             if h in all_hpo.get(hpo):
    #                 return True
    #         return False

    def __str__(self) -> str:
        return f"Patient(patient_id:{self.patient_id}, " \
            f"variants:{[var.variant_string for var in self.variants]}, " \
            f"phenotypes:{[pheno.identifier for pheno in self.phenotypes]}, " \
            f"proteins:{[prot.protein_id for prot in self.proteins]})"

    def __repr__(self) -> str:
        return str(self)

    def __eq__(self, other) -> bool:
        return isinstance(other, Patient) \
            and self.patient_id == other.patient_id \
            and self.variants == other.variants \
            and self.phenotypes == other.phenotypes \
            and self.proteins == other.proteins
    
    def __hash__(self) -> int:
        return hash((self.patient_id, self.variants, self.phenotypes, self.proteins))
