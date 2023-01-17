import pyensembl
import varcode as vc
from .patient_class import Patient
import re
import warnings

    
def is_var_type(pat, varType):
    """"" Determines if the given Patient has a variant of the given variant type

    Args:
        pat (Patient) : Patient Class
        varType (str) : Options include - 
                        missense
                        nonsense
                        duplication
                        deletion
                        insertion
                        transition
                        transversion
                        indel

    Returns:
        boolean :   True if the Patient does have a variant with given variant type
    """""

    if pat.variant.variant_type == varType:
        return True
    else:
        return False

def is_not_var_type(pat, varType):
    """"" Determines if the given Patient does not have a variant of the given variant type

    Args:
        pat (Patient) : Patient Class
        varType (str) : Options include - 
                        missense
                        nonsense
                        duplication
                        deletion
                        insertion
                        transition
                        transversion
                        indel

    Returns:
        boolean :   True if the Patient does NOT have a variant with given variant type
    """""

    if pat.variant.variant_type == varType:
        return False
    else:
        return True


def is_var_match(pat, variant):
    """ Determines if the given Patient has the given Variant

    Args:
        pat (Patient) : Patient Class
        variant (str OR Variant) : Either a string formatted 
                                'chr:start:reference:alternative'
                                OR a class Variant

    Returns:
        boolean : True if the Patient does have the given Variant
    
    """

    if pat.variant.variant is not None:
        if isinstance(variant, str):
            test_var = verify_var(variant)
        else:
            test_var = variant.variant
        if pat.variant.variant == test_var:
            return True
        else:
            return False
    else:
        return False

def is_not_var_match(pat, variant):
    """ Determines if the given Patient does NOT have the given Variant

    Args:
        pat (Patient) : Patient Class
        variant (str OR Variant) : Either a string formatted 
                                'chr:start:reference:alternative'
                                OR a class Variant

    Returns:
        boolean : True if the Patient does NOT have the given Variant
    
    """

    if pat.variant.variant is not None:
        if isinstance(variant, str):
            test_var = verify_var(variant)
        else:
            test_var = variant.variant
        if pat.variant.variant == test_var:
            return False
        else:
            return True
    else:
        return False

def verify_var(variant):
    """
    Args:
        variant (str) : 'chr:start:reference:alternative'
            i.e.    - '3:12345:A:G'
                    - '3:15432:AG:A'
                    - '3:98765:A:AG'

    Returns:
        Variant : a variant of class Variant 
    """

    contig, start, ref, alt = variant.split(':')
    var = vc.Variant(contig, start, ref, alt, ensembl = pyensembl.ensembl_grch37)
    return var

def in_feature(pat, feature):
    """Given a specific patient and feature, determine True
    or False that the variant effect is within that feature
    
    Args:
        pat (Patient) : Patient Class 
        feature (str) : Either a feature ID or a feature type
                        feature types include:
                                Domain
                                Region
                                Motif
                                Repeat
    Returns:
        boolean : True if the variant location is in the 
                    specified feature or feature type
    """

    featureDF = pat.protein.features
    loc = pat.variant.protein_effect_location
    isIn = False
    if loc is not None and not featureDF.empty:
        for row in featureDF.iterrows():
            if row[0] == feature or row[1]['type'] == feature:
                if row[1]['start'] <= loc <= row[1]['end']:
                    isIn = True
    return isIn

def not_in_feature(pat, feature):
    """ Given a specific patient and feature, determine True
    or False that the variant effect is NOT within that feature

    Args:
        pat (Patient) : Patient Class 
        feature (str) : Either a feature ID or a feature type
                        feature types include:
                                Domain
                                Region
                                Motif
                                Repeat
    Returns:
        boolean : True if the variant location is NOT in the 
                    specified feature or feature type
    """

    featureDF = pat.protein.features
    loc = pat.variant.protein_effect_location
    isIn = False
    if loc is not None and not featureDF.empty:
        for row in featureDF.iterrows():
            if row[0] == feature or row[1]['type'] == feature:
                if not row[1]['start'] <= loc <= row[1]['end']:
                    isIn = True
    return isIn
