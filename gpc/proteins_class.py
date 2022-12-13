import requests

class Protein:
    def __init__(self, transcript):
        self._transcript = transcript

        self._protein_id = transcript.protein_id
        self._protein_json = Protein.findProtein(transcript.protein_id)
        if not self.__verify_isoforms():
            print("WARNING - Current Protein ID or Transcript ID does not match with Described Isoform")
        self._protein_name = self._protein_json.get('proteinDescription').get('recommendedName').get('fullName').get('value')
        self._features = self.__seperate_by_feature()

    @staticmethod
    def findProtein(protein_id):
        fields = "accession,id,gene_names,gene_primary,protein_name,feature_count,ft_domain,ft_motif,ft_region,ft_repeat,ft_zn_fing,ft_chain,ft_coiled,ft_compbias,ft_mod_res,ft_crosslnk,ft_var_seq,ft_variant,ft_conflict,ft_helix,xref_ensembl,cc_alternative_products"
        url = 'https://rest.uniprot.org/uniprotkb/search?query=xref:ensembl-'+ protein_id +'&fields=' + fields
        protJson = requests.get(url).json()['results'][0]
        if protJson is not None:
            return protJson
        else:
            print('WARNING - Protein Not Found')
            return None

    def __verify_isoforms(self):
        result = False
        try:
            Isoforms = self._protein_json.get('comments')[0].get('isoforms')
        except (IndexError):
            print('No Isoforms found, no need to verify')
            return True
        defaultIso = None
        for iso in Isoforms:
            if iso.get('isoformSequenceStatus') == 'Displayed':
                defaultIso = iso.get('isoformIds')[0]
                break
        ref = self._protein_json.get('uniProtKBCrossReferences')
        for r in ref:
            if r.get('database') == "Ensembl" and r.get('id').split('.')[0] == self._transcript.id and r.get('isoformId') == defaultIso:
                for k in r.get('properties'):
                    if k.get('key') == 'ProteinId' and k.get('value').split('.')[0] == self.protein_id:
                        result = True
        return result

    def __seperate_by_feature(self):
        outputDict = {}
        for feat in self._protein_json.get('features'):
            if feat.get('location').get('start').get('modifier') == 'EXACT':
                finalFeat = {'start':feat.get('location').get('start').get('value'),'end':feat.get('location').get('end').get('value')}
                if feat.get('type') in outputDict.keys():
                    outputDict[feat.get('type')].append(finalFeat)
                else:
                    outputDict[feat.get('type')] = []
                    outputDict[feat.get('type')].append(finalFeat)
        return outputDict


    @property
    def protein_id(self):
        return self._protein_id

    @property
    def protein_name(self):
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


