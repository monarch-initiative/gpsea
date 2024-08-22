import abc
import enum
import typing

import hpotk

from .genome import Region


class FeatureInfo:
    """
    `FeatureInfo` represents a protein feature
    (e.g. a repeated sequence given the name "ANK 1" in protein "Ankyrin repeat domain-containing protein 11")
"""

    def __init__(self, name: str, region: Region):
        self._name = hpotk.util.validate_instance(name, str, 'name')
        self._region = hpotk.util.validate_instance(region, Region, 'region')

    @property
    def name(self) -> str:
        """
        Returns:
            string: the name of the protein feature
        """
        return self._name

    @property
    def start(self) -> int:
        """
        Returns:
            integer: A 0-based (excluded) start coordinate of the protein feature.
        """
        return self._region.start

    @property
    def end(self) -> int:
        """
        Returns:
            integer: A 0-based (included) end coordinate of the protein feature.
        """
        return self._region.end

    @property
    def region(self) -> Region:
        """
        Returns:
            Region: a protein region spanned by the feature.
        """
        return self._region

    def __len__(self):
        return len(self._region)

    def __eq__(self, other) -> bool:
        return isinstance(other, FeatureInfo) \
            and self.name == other.name \
            and self.region == other.region

    def __hash__(self):
        return hash((self._name, self._region))

    def __str__(self) -> str:
        return f"FeatureInfo(name={self.name}, start={self.start}, end={self.end})"

    def __repr__(self) -> str:
        return str(self)


class FeatureType(enum.Enum):
    """
    An enum representing the protein feature types supported in GPSEA.
    """
    
    REPEAT = enum.auto()
    """
    A repeated sequence motif or repeated domain within the protein.
    """
    
    MOTIF = enum.auto()
    """
    A short (usually not more than 20 amino acids) conserved sequence motif of biological significance.
    """

    DOMAIN = enum.auto()
    """
    A specific combination of secondary structures organized into a characteristic three-dimensional structure or fold.
    """

    REGION = enum.auto()
    """
    A region of interest that cannot be described in other subsections.
    """


class ProteinFeature(metaclass=abc.ABCMeta):

    @staticmethod
    def create(info: FeatureInfo, feature_type: FeatureType):
        return SimpleProteinFeature(info, feature_type)

    @property
    @abc.abstractmethod
    def info(self) -> FeatureInfo:
        pass

    @property
    @abc.abstractmethod
    def feature_type(self) -> FeatureType:
        pass

    def to_string(self) -> str:
        return f"{self.feature_type.name}-{self.info.name}-{self.info.region}"


class SimpleProteinFeature(ProteinFeature):
    """
    An implementation of a `ProteinFeature`.
    """
    # Not part of the public API.

    def __init__(self, info: FeatureInfo, feature_type: FeatureType):
        """Constructs all necessary attributes for a SimpleProteinFeature

        Args:
            info (FeatureInfo): A FeatureInfo object, describing name and location of the feature
            feature_type (FeatureType): A FeatureType object, limited to specific feature types
        """
        if not isinstance(info, FeatureInfo):
            raise ValueError(f"info must be type FeatureInfo but was type {type(info)}")
        self._info = info
        if not isinstance(feature_type, FeatureType):
            raise ValueError(f"feature_type must be type FeatureType but was type {type(feature_type)}")
        self._type = feature_type

    @property
    def info(self) -> FeatureInfo:
        """
        Returns:
            FeatureInfo: A FeatureInfo object, describing name and location of the feature
        """
        return self._info

    @property
    def feature_type(self) -> FeatureType:
        """
        Returns:
            FeatureType: A FeatureType object, limited to specific feature types (e.g. REGION, REPEAT, MOTIF, DOMAIN)
        """
        return self._type

    def __str__(self) -> str:
        return f"SimpleProteinFeature(type={self._type}, " \
               f"info={self._info})"

    def __repr__(self) -> str:
        return str(self)

    def __eq__(self, other) -> bool:
        return isinstance(other, SimpleProteinFeature) \
            and self._type == other._type \
            and self._info == other._info

    def __hash__(self) -> int:
        return hash((self._type, self._info))


class ProteinMetadata:
    """
    An info regarding a protein sequence, including an ID, a label,
    and location of protein features, such as motifs, domains, or other regions.
    """

    def __init__(
        self,
        protein_id: str,
        label: str,
        protein_features: typing.Sequence[ProteinFeature],
        protein_length: int = 0,
    ):
        if not isinstance(protein_id, str):
            raise ValueError(f"Protein ID must be type string but is type {type(protein_id)}")
        self._id = protein_id
        if not isinstance(label, str):
            raise ValueError(f"Protein label must be type string but is type {type(label)}")
        self._label = label
        if not all(isinstance(x, ProteinFeature) for x in protein_features):
            raise ValueError(
                f"Protein Features must be a list of type ProteinFeature but is type {type(protein_features)}")
        self._features = tuple(protein_features)
        self._protein_length = protein_length

    @property
    def protein_id(self) -> str:
        """
        Returns:
            string: A string unique to this protein
        """
        return self._id

    @property
    def label(self) -> str:
        """
        Returns:
            string: The full name of the protein
        """
        return self._label

    @property
    def protein_features(self) -> typing.Sequence[ProteinFeature]:
        """
        Returns:
            Sequence[ProteinFeature]: A sequence of ProteinFeatures objects
        """
        return self._features
    
    @property
    def protein_length(self) -> int:
        """
        Returns:
            int: length of protein
        """
        return self._protein_length

    def domains(self) -> typing.Iterable[ProteinFeature]:
        """
        Returns:
            Iterable[ProteinFeature]: A subgroup of protein_features, where the ProteinFeature object has a FeatureType equal to "DOMAIN"
        """
        return filter(lambda f: f.feature_type == FeatureType.DOMAIN, self.protein_features)

    def repeats(self) -> typing.Iterable[ProteinFeature]:
        """
        Returns:
            Iterable[ProteinFeature]: A subgroup of protein_features, where the ProteinFeature object has a FeatureType equal to "REPEAT"
        """
        return filter(lambda f: f.feature_type == FeatureType.REPEAT, self.protein_features)

    def regions(self) -> typing.Iterable[ProteinFeature]:
        """
        Returns:
            Iterable[ProteinFeature]: A subgroup of protein_features, where the ProteinFeature object has a FeatureType equal to "REGIONS"
        """
        return filter(lambda f: f.feature_type == FeatureType.REGION, self.protein_features)

    def motifs(self) -> typing.Iterable[ProteinFeature]:
        """
        Returns:
            Iterable[ProteinFeature]: A subgroup of protein_features, where the ProteinFeature object has a FeatureType equal to "MOTIF"
        """
        return filter(lambda f: f.feature_type == FeatureType.MOTIF, self.protein_features)

    def get_features_variant_overlaps(self, region: Region) -> typing.Collection[ProteinFeature]:
        """
        Get a collection of protein features that overlap with the `region`.
        Args:
            region: the query region.

        Returns:
            Collection[ProteinFeature]: a collection of overlapping protein features.
        """
        affected_features = set()
        for feat in self.protein_features:
            if feat.info.region.overlaps_with(region):
                affected_features.add(feat)

        return affected_features

    def __str__(self) -> str:
        return f"ProteinMetadata(id={self.protein_id}, " \
               f"label={self.label}, " \
               f"features={str(self.protein_features)})"

    def __eq__(self, other) -> bool:
        return isinstance(other, ProteinMetadata) \
            and self.label == other.label \
            and self.protein_features == other.protein_features \
            and self.protein_id == other.protein_id

    def __hash__(self) -> int:
        return hash((self.protein_id, self.label, self._features))

    def __repr__(self) -> str:
        return str(self)
