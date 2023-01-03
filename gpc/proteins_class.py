import requests
import pandas as pd
import re

class Protein:
    def __init__(self, protein_id):
        self._protein_id = protein_id
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
            self._features = pd.DataFrame()

    @staticmethod
    def findProtein(protein_id):
        if protein_id is not None:
            fields = "accession,id,gene_names,gene_primary,protein_name,feature_count,ft_domain,ft_motif,ft_region,ft_repeat,xref_ensembl,cc_alternative_products" 
            # Additional features we could add - ,ft_zn_fing,ft_chain,ft_coiled,ft_compbias,ft_mod_res,ft_crosslnk,ft_var_seq,ft_variant,ft_conflict,ft_helix
            url = 'https://rest.uniprot.org/uniprotkb/search?query=xref:ensembl-'+ protein_id +'&fields=' + fields
            json = requests.get(url).json()
            if len(json['results']) > 0:
                protein = requests.get(url).json()['results'][0]
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

    def __seperate_by_feature(self):
        if self._protein_json is not None:
            outputDict = {}
            try:
                for feat in self._protein_json.get('features'):
                    if feat.get('location').get('start').get('modifier') == 'EXACT':
                        key = feat.get('type') + ': ' + feat.get('description')
                        i = 1
                        startKey = key
                        while key in outputDict.keys():
                            key = startKey + " " + str(i)
                            i+=1
                        outputDict[key] = {'type': feat.get('type'), 'start':feat.get('location').get('start').get('value'),'end':feat.get('location').get('end').get('value')}
                finalOut = pd.DataFrame.from_dict(outputDict)
                return finalOut.T
            except(TypeError):
                print('Protein has no features')
                return pd.DataFrame()
        else:
            return pd.DataFrame()

    def add_vars_to_features(self, varSeries):
        self._features = pd.concat([self._features,varSeries], axis=1)
        #return self._features


    @property
    def id(self):
        return self._protein_id

    @property
    def label(self):
        return self._protein_name

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



            