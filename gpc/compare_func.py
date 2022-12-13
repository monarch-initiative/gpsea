import pyensembl
import varcode as vc
from .patient_class import Patient
import re

    
def is_missense(pat, holder = None):
    if pat.variant.variant.effects().top_priority_effect().short_description.endswith("*") or not pat._variant.variant.is_snv:
        return False
    else:
        return True
    
def is_nonsense(pat, holder = None):
    if pat.variant.variant.effects().top_priority_effect().short_description.endswith("*") and pat._variant.variant.is_snv:
        return True
    else:
        return False
        
def is_deletion(pat, holder = None):
    if pat.variant.variant is not None:
        return pat.variant.variant.is_deletion
    elif pat.disease_label is not None:
        if 'deletion' in pat.disease_label:
            return True
    else:
        return False

def is_insertion(pat, holder = None):
    return pat.variant.variant.is_insertion

def is_transition(pat, holder = None):
    return pat.variant.variant.is_transition

def is_transversion(pat, holder = None):
    return pat.variant.variant.is_transversion

def is_indel(pat, holder = None):
    return pat.variant.variant.is_indel

def is_duplication(pat, holder = None):
    if pat.variant.variant is not None:
        if 'dup' in pat.variant.variant.effects().top_priority_effect().short_description:
            return True
        else:
            return False
    elif pat.disease_label is not None:
        if 'duplication' in pat.disease_label:
            return True
    else:
        return False

def is_var_match(pat, variant):
    if pat.variant.variant is not None:
        if not type(variant) == vc.Variant:
            test_var = verify_var(variant)
        else:
            test_var = variant
        if pat.variant.variant == test_var:
            return True
        else:
            return False
    else:
        return False

def verify_var(variant):
    """
    Variant is a String formatted:
    'chr:start:reference:alternative"
    i.e.    - '3:12345:A:G'
            - '3:15432:AG:A'
            - '3:98765:A:AG'
    """
    contig, start, ref, alt = variant.split(':')
    var = vc.Variant(contig, start, ref, alt, ensembl = pyensembl.ensembl_grch37)
    return var

def in_domain(pat, domain):
    loc = int(re.sub(r'[^0-9]', '', pat.variant.top_effect.short_description))
    if domain[0] <= loc <= domain[1]:
        return True
    return False

def is_motif_match(pat, motif):
    loc = int(re.sub(r'[^0-9]', '', pat.variant.top_effect.short_description))
    if motif[0] <= loc <= motif[1]:
            return True
    return False

