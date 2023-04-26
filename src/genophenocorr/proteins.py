import enum
import abc
import logging
import typing

import requests


class ProteinRegion:

    def __init__(self, start: int, end: int):
        # TODO(lnrekerle) - check instance and that it is non-negative and start <= end
        self._start = start
        self._end = end

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

    # TODO - eq, hash, str, repr, __gt__


class ProteinFeatureType(enum.Enum):
    # TODO - add python doc strings with a little description, if possible.
    REPEAT = enum.auto()
    MOTIF = enum.auto()
    DOMAIN = enum.auto()
    REGION = enum.auto()


class ProteinFeature(metaclass=abc.ABCMeta):

    @property
    @abc.abstractmethod
    def region(self) -> ProteinRegion:
        pass

    @property
    @abc.abstractmethod
    def feature_type(self) -> ProteinFeatureType:
        pass


class SimpleProteinFeature(ProteinFeature):

    def __init__(self, region: ProteinRegion, feature_type: ProteinFeatureType):
        # TODO - add instance checks
        self._region = region
        self._type = feature_type

    @property
    def region(self) -> ProteinRegion:
        return self._region

    @property
    def feature_type(self) -> ProteinFeatureType:
        return self._type

    # TODO - eq, hash, str, repr


class ProteinMetadata:

    def __init__(self, protein_id: str, label: str, protein_features: typing.Sequence[ProteinFeature]):
        # TODO - add instance checks
        self._id = protein_id
        self._label = label
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
        return list(filter(lambda f: f.feature_type == ProteinFeatureType.DOMAIN, self.protein_features))

    # TODO(lnrekerle) - implement for other enum types if you think it is appropriate

    # TODO - eq, hash, str, repr


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

    def annotate(self, protein_id: str) -> typing.Sequence[ProteinMetadata]:
        if not isinstance(protein_id, str):
            # TODO - log a warning?
            self._logger.debug(f'Protein ID must be a str but it was {type(protein_id)}')
            return []

        # TODO - implement using the code below

        return []

# TODO - caching protein annotations?


class Protein:
    def __init__(self, protein_id, variantss_in_protein):
        self._protein_id = protein_id
        self._variants_in_protein = variantss_in_protein
        self._protein_json = Protein.findProtein(protein_id)
        if self._protein_json is not None:
            if not self.__verify_isoforms():
                print("WARNING - Current Protein ID or Transcript ID does not match with Described Isoform")
            if self._protein_json.get('proteinDescription') is None or len(self._protein_json.get('proteinDescription')) == 0:
                self._protein_name = self._protein_id
            else:
                if self._protein_json.get('proteinDescription').get('recommendedName') is None: 
                    if self._protein_json.get('proteinDescription').get('submissionNames') is None:
                        self._protein_name = self._protein_id
                    elif len(self._protein_json.get('proteinDescription').get('submissionNames')) > 0:
                        self._protein_name = self._protein_json.get('proteinDescription').get('submissionNames')[0].get('fullName').get('value')
                elif len(self._protein_json.get('proteinDescription').get('recommendedName')) > 0:
                        self._protein_name = self._protein_json.get('proteinDescription').get('recommendedName').get('fullName').get('value') 
                else:
                    self._protein_name = self._protein_id
            self._features = self.__seperate_by_feature()
        else:
            self._protein_name = None
            # self._features = pd.DataFrame()

    @staticmethod
    def findProtein(protein_id):
        if protein_id is not None:
            fields = "accession,id,gene_names,gene_primary,protein_name,feature_count,ft_domain,ft_motif,ft_region,ft_repeat,xref_ensembl,cc_alternative_products" 
            # Additional features we could add - ,ft_zn_fing,ft_chain,ft_coiled,ft_compbias,ft_mod_res,ft_crosslnk,ft_var_seq,ft_variant,ft_conflict,ft_helix
            url = 'https://rest.uniprot.org/uniprotkb/search?query='+ protein_id +'&fields=' + fields
            try:
                json = requests.get(url, timeout=10).json()
            except requests.exceptions.Timeout:
                print(f'Timeout Occurred with Protein ID {protein_id}')
                protein = None
            if len(json['results']) > 0:
                protein = json['results'][0]
            else:
                print(f'No Protein found with ID {protein_id}')
                protein = None
        else:
            protein = None
        return protein

    def __verify_isoforms(self):
        result = False
        if self._protein_json.get('comments') is not None:
            if len(self._protein_json.get('comments')) > 0:
                Isoforms = self._protein_json.get('comments')[0].get('isoforms')
            else:
                return True
        else:
            return False
        defaultIso = None
        for iso in Isoforms:
            if iso.get('isoformSequenceStatus') == 'Displayed':
                defaultIso = iso.get('isoformIds')[0]
                break
        ref = self._protein_json.get('uniProtKBCrossReferences')
        for r in ref:
            if r.get('database') == "Ensembl" and r.get('isoformId') == defaultIso:
                for k in r.get('properties'):
                    if k.get('key') == 'ProteinId' and k.get('value').split('.')[0] == self.id:
                        result = True
        return result

    # def __seperate_by_feature(self):
    #     if self._protein_json is not None:
    #         outputDict = {}
    #         try:
    #             for feat in self._protein_json.get('features'):
    #                 if feat.get('location').get('start').get('modifier') == 'EXACT':
    #                     key = feat.get('type') + ': ' + feat.get('description')
    #                     i = 1
    #                     startKey = key
    #                     while key in outputDict.keys():
    #                         key = startKey + " " + str(i)
    #                         i+=1
    #                     outputDict[key] = {'type': feat.get('type'), 'start':feat.get('location').get('start').get('value'),'end':feat.get('location').get('end').get('value')}
    #             finalOut = pd.DataFrame.from_dict(outputDict)
    #             return finalOut.T
    #         except(TypeError):
    #             print('Protein has no features')
    #             return pd.DataFrame()
    #     else:
    #         return pd.DataFrame()
    #
    # def add_vars_to_features(self, varSeries):
    #     self._features = pd.concat([self._features,varSeries], axis=1)
    #     #return self._features


    @property
    def id(self):
        return self._protein_id

    @property
    def label(self):
        return self._protein_name

    @property
    def variants_in_protein(self):
        return self._variants_in_protein

    @property
    def features(self):
        return self._features

    @property
    def domains(self):
        if 'Domain' in self.features.keys():
            return self.features.get('Domain')
        else:
            return None

    @property
    def regions(self):
        if 'Region' in self.features.keys():
            return self.features.get('Region')
        else:
            return None

    @property
    def motifs(self):
        if 'Motif' in self.features.keys():
            return self.features.get('Motif')
        else:
            return None

    @property
    def repeats(self):
        if 'Repeat' in self.features.keys():
            return self.features.get('Repeat')
        else:
            return None



            