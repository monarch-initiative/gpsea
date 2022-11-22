import pyensembl
import varcode as vc
from .patient_class import Patient

    
def is_missense(pat):
    if pat.variant.variant.effects().top_priority_effect().short_description.endswith("*") or not pat._variant.variant.is_snv:
        return False
    else:
        return True
    
def is_nonsense(pat):
    if pat.variant.variant.effects().top_priority_effect().short_description.endswith("*") and pat._variant.variant.is_snv:
        return True
    else:
        return False
        
def is_deletion(pat):
    if pat.variant.variant is not None:
        return pat.variant.variant.is_deletion
    elif pat.disease_label is not None:
        if 'deletion' in pat.disease_label:
            return True
    else:
        return False

def is_insertion(pat):
    return pat.variant.variant.is_insertion

def is_transition(pat):
    return pat.variant.variant.is_transition

def is_transversion(pat):
    return pat.variant.variant.is_transversion

def is_indel(pat):
    return pat.variant.variant.is_indel

def is_duplication(pat):
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

