import pyensembl
import varcode as vc
from .patient._patient_data import Patient
import re
import warnings


ALLOWED_VARIANT_TYPES = {'stop_gained',
                        'frameshift_variant',
                        'missense_variant',
                        'SNV',
                        'deletion',
                        'insertion',
                        'copy_number_change',
                        'copy_number_decrease',
                        'copy_number_increase'}


def is_var_type(var, pat, varType):
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
    if varType not in ALLOWED_VARIANT_TYPES:
        raise ValueError(f"Did not recognize variant type {varType}")
    if varType in var.variant_types:
        return True
    else:
        return False

def is_not_var_type(var, pat, varType):
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
    if varType not in ALLOWED_VARIANT_TYPES:
        raise ValueError(f"Did not recognize variant type {varType}")
    if varType in var.variant_types:
        return False
    else:
        return True


def is_var_match(var, pat, variant):
    """ Determines if the given Patient has the given Variant

    Args:
        pat (Patient) : Patient Class
        variant (str OR Variant) : Either a string formatted 
                                'chr:start:reference:alternative'
                                OR a class Variant

    Returns:
        boolean : True if the Patient does have the given Variant
    
    """
    if variant == var.variant_string:
        return True
    else:
        return False

def is_not_var_match(var, pat, variant):
    """ Determines if the given Patient does NOT have the given Variant

    Args:
        pat (Patient) : Patient Class
        variant (str OR Variant) : Either a string formatted 
                                'chr:start:reference:alternative'
                                OR a class Variant

    Returns:
        boolean : True if the Patient does NOT have the given Variant
    
    """
    if variant == var.variant_string:
        return False
    else:
        return True

def in_exon(var, pat, exon):
    if exon in var.exons_effected:
        return True
    else:
        return False

def not_in_exon(var, pat, exon):
    if exon in var.exons_effected:
        return False
    else:
        return True


def in_feature(var, pat, feature):
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
    for prot in pat.proteins:
        featureDF = prot.features
        loc = var.protein_effect_location
        if loc is not None and not featureDF.empty:
            for row in featureDF.iterrows():
                if row[0] == feature or row[1]['type'] == feature:
                    if row[1]['start'] <= loc <= row[1]['end']:
                        return True
    return False

def not_in_feature(var, pat, feature):
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
    for prot in pat.proteins:
        featureDF = prot.features
        loc = var.protein_effect_location
        if loc is not None and not featureDF.empty:
            for row in featureDF.iterrows():
                if row[0] == feature or row[1]['type'] == feature:
                    if row[1]['start'] <= loc <= row[1]['end']:
                        return False
    return True
