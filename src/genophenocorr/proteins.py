import enum
import abc
import logging
import typing

import requests


class FeatureInfo:

    def __init__(self, name: str, start: int, end: int):
        if not isinstance(name, str):
            raise ValueError(f"name must be type string but was type {type(name)}")
        self._name = name
        if not isinstance(start, int):
            raise ValueError(f"start must be an integer but was type {type(start)}")
        self._start = start
        if not isinstance(end, int):
            raise ValueError(f"end must be an integer but was type {type(end)}")
        self._end = end

        if self._start > self._end:
            raise ValueError(f"The start value must come before end but {self._start} is greater than {self._end}")

    @property
    def name(self) -> str:
        return self._name

    @property
    def start(self) -> int:
        """
        Get a 0-based (excluded) start coordinate of the protein region.
        """
        return self._start

    @property
    def end(self) -> int:
        """
        Get a 0-based (included) end coordinate of the protein region.
        """
        return self._end

    def __len__(self):
        return self._end - self._start

    def __eq__(self, other) -> bool:
        return isinstance(other, FeatureInfo) \
            and self.name == other.name \
            and self.start == other.start \
            and self.end == other.end

    def __str__(self) -> str:
        return f"FeatureInfo(name={self.name}, start={self.start}, end={self.end})"
    
    def __repr__(self) -> str:
        return str(self)


class FeatureType(enum.Enum):
    """
    
    """
    REPEAT = enum.auto()
    MOTIF = enum.auto()
    DOMAIN = enum.auto()
    REGION = enum.auto()


class ProteinFeature(metaclass=abc.ABCMeta):

    @property
    @abc.abstractmethod
    def info(self) -> FeatureInfo:
        pass

    @property
    @abc.abstractmethod
    def feature_type(self) -> FeatureType:
        pass


class SimpleProteinFeature(ProteinFeature):

    def __init__(self, info: FeatureInfo, feature_type: FeatureType):
        if not isinstance(info, FeatureInfo):
            raise ValueError(f"info must be type FeatureInfo but was type {type(info)}")
        self._info = info
        if not isinstance(feature_type, FeatureType):
            raise ValueError(f"feature_type must be type FeatureType but was type {type(feature_type)}")
        self._type = feature_type

    @property
    def info(self) -> FeatureInfo:
        return self._info

    @property
    def feature_type(self) -> FeatureType:
        return self._type

    def __str__(self) -> str:
        return f"SimpleProteinFeature(type={self.feature_type}, " \
            f"info={str(self.info)})"

    def __repr__(self) -> str:
        return str(self)

    def __eq__(self, other) -> bool:
        return isinstance(other, ProteinFeature) \
            and self.feature_type == other.feature_type \
            and self.info == other.info

    def __hash__(self) -> int:
        return hash((self.feature_type, self.info))


class ProteinMetadata:

    def __init__(self, protein_id: str, label: str, protein_features: typing.Sequence[ProteinFeature]):
        if not isinstance(protein_id, str):
            raise ValueError(f"Protein ID must be type string but is type {type(protein_id)}")
        self._id = protein_id
        if not isinstance(label, str):
            raise ValueError(f"Protein label must be type string but is type {type(label)}")
        self._label = label
        if not all(isinstance(x, ProteinFeature) for x in protein_features):
            raise ValueError(f"Protein Features must be a list of type ProteinFeature but is type {type(protein_features)}")
        self._features = protein_features

    @property
    def protein_id(self) -> str:
        return self._id

    @property
    def label(self) -> str:
        return self._label

    @property
    def protein_features(self) -> typing.Sequence[ProteinFeature]:
        return self._features

    def domains(self) -> typing.Sequence[ProteinFeature]:
        return list(filter(lambda f: f.feature_type == FeatureType.DOMAIN, self.protein_features))

    def repeats(self) -> typing.Sequence[ProteinFeature]:
        return list(filter(lambda f: f.feature_type == FeatureType.REPEAT, self.protein_features))

    def regions(self) -> typing.Sequence[ProteinFeature]:
        return list(filter(lambda f: f.feature_type == FeatureType.REGION, self.protein_features))

    def motifs(self) -> typing.Sequence[ProteinFeature]:
        return list(filter(lambda f: f.feature_type == FeatureType.MOTIF, self.protein_features))

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
        return hash((self.protein_id, self.label, self.protein_features))

    def __repr__(self) -> str:
        return str(self)


class ProteinMetadataService(metaclass=abc.ABCMeta):

    @abc.abstractmethod
    def annotate(self, protein_id: str) -> typing.Sequence[ProteinMetadata]:
        """
        Get metadata for given protein ID (e.g. NP_037407.4).

        The function returns an empty sequence if there is no info for the protein ID.
        """
        pass


class UniprotProteinMetadataService(ProteinMetadataService):

    def __init__(self):
        self._logger = logging.getLogger(__name__)
        self._url = 'https://rest.uniprot.org/uniprotkb/search?query=%s&fields=accession,id,gene_names,gene_primary,protein_name,ft_domain,ft_motif,ft_region,ft_repeat,xref_refseq'

    def annotate(self, protein_id: str) -> typing.Sequence[ProteinMetadata]:
        if not isinstance(protein_id, str):
            # TODO - log a warning?
            raise ValueError(f'Protein ID must be a str but it was {type(protein_id)}')
        api_url = self._url % (protein_id)
        r = requests.get(api_url).json()
        results = r['results']
        if len(results) == 0:
            self._logger.warning(f"No proteins found for ID {protein_id}. Please verify refseq ID.")
            return []
        protein_list = []
        for protein in results:
            protein_name = protein['proteinDescription']['recommendedName']['fullName']['value']
            all_features_list = []
            for feature in protein['features']:
                feat_type = feature['type']
                feat_name = feature['description']
                feat_start = int(feature['location']['start']['value'])
                feat_end = int(feature['location']['end']['value'])
                feat = SimpleProteinFeature(FeatureInfo(feat_name, feat_start, feat_end), feat_type.upper())
                all_features_list.append(feat)
            protein_list.append(ProteinMetadata(protein_id, protein_name, all_features_list))

        return protein_list

# TODO - caching protein annotations?

            