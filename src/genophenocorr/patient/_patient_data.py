import typing

from genophenocorr.phenotype import Phenotype
from genophenocorr.variant import Variant
from genophenocorr.protein import ProteinMetadata


class Patient:
    """A class that represents an individual patient

    Attributes:
        patient_id (string): A string unique to this Patient object
        phenotypes (Sequence[Phenotype]): A list of Phenotype objects
        variants (Sequence[Variant]): A list of Variant objects
        proteins (Sequence[ProteinMetadata]): A list of ProteinMetadata objects
    """
    def __init__(self, patient_id: str,
                 phenotypes: typing.Sequence[Phenotype],
                 variants: typing.Sequence[Variant],
                 proteins: typing.Sequence[ProteinMetadata]):
        """Constructs all necessary attributes for a Patient object

        Args:
            patient_id (string): A string unique to this Patient object
            phenotypes (Sequence[Phenotype]): A list of Phenotype objects
            variants (Sequence[Variant]): A list of Variant objects
            proteins (Sequence[ProteinMetadata]): A list of ProteinMetadata objects
        """
        self._id = patient_id
        self._phenotypes = tuple(phenotypes)
        self._variants = tuple(variants)
        self._proteins = tuple(proteins)

    @property
    def patient_id(self) -> str:
        """
        Returns:
            string: Patient ID unique to this Patient object
        """
        return self._id

    @property
    def phenotypes(self) -> typing.Sequence[Phenotype]:
        """
        Returns:
            Sequence[Phenotype]: A list of Phenotype objects associated with this Patient object
        """
        return self._phenotypes
    
    @property
    def variants(self) -> typing.Sequence[Variant]:
        """
        Returns:
            Sequence[Variant]: A list of Variant objects associated with this Patient object
        """
        return self._variants

    @property
    def proteins(self) -> typing.Sequence[ProteinMetadata]:
        """
        Returns:
            Sequence[ProteinMetadata]: A list of ProteinMetadata objects associated with this Patient object
        """
        return self._proteins

    def __str__(self) -> str:
        return f"Patient(patient_id:{self.patient_id}, " \
            f"variants:{[var.variant_string for var in self.variants]}, " \
            f"phenotypes:{[pheno.identifier for pheno in self.phenotypes]}, " 
            #f"proteins:{[prot.protein_id for prot in self.proteins]})"

    def __repr__(self) -> str:
        return str(self)

    def __eq__(self, other) -> bool:
        return isinstance(other, Patient) \
            and self.patient_id == other.patient_id \
            and self.variants == other.variants \
            and self.phenotypes == other.phenotypes 
            #and self.proteins == other.proteins
    
    def __hash__(self) -> int:
        return hash((self.patient_id, self.variants, self.phenotypes)) #self.proteins))
