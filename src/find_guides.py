"""
Functions to find CRISPR guide RNAs in DNA sequences
"""

from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd


def find_pam_sites(sequence, pam_sequence="GG"):
    """
    Find all PAM sites (NGG) in a DNA sequence.
    
    Args:
        sequence (str): DNA sequence to search
        pam_sequence (str): PAM pattern (default: "GG")
    
    Returns:
        list: Positions of PAM sites
    """
    sequence = str(sequence).upper() 
    pam_sites = []
    
    # Search for PAM (NGG = any base + GG)
    for i in range(len(sequence) - 2):
        if sequence[i+1:i+3] == pam_sequence:
            pam_sites.append(i)
    
    return pam_sites


def extract_guide_sequence(sequence, pam_position, guide_length=20):
    """
    Extract 20bp guide RNA upstream of PAM.
    
    Args:
        sequence (str): DNA sequence
        pam_position (int): Position of PAM
        guide_length (int): Length of guide (default: 20)
    
    Returns:
        str or None: Guide sequence
    """
    sequence = str(sequence).upper()
    
    guide_start = pam_position - guide_length
    guide_end = pam_position
    
    if guide_start < 0:
        return None
    
    guide = sequence[guide_start:guide_end]
    
    if len(guide) != guide_length:
        return None
    
    return guide


def find_all_guides(sequence):
    """
    Find all possible guide RNAs in a sequence.
    
    Args:
        sequence (str): DNA sequence
    
    Returns:
        pandas.DataFrame: All guides with positions
    """
    sequence = str(sequence).upper()
    guides = []
    
    # Forward strand
    pam_sites = find_pam_sites(sequence)
    
    for pam_pos in pam_sites:
        guide = extract_guide_sequence(sequence, pam_pos)
        
        if guide:
            pam = sequence[pam_pos:pam_pos+3]
            
            guides.append({
                'guide_sequence': guide,
                'pam_site': pam_pos,
                'pam_sequence': pam,
                'strand': '+',
                'full_target': guide + pam
            })
    
    # Reverse complement
    rev_seq = str(Seq(sequence).reverse_complement())
    rev_pam_sites = find_pam_sites(rev_seq)
    
    for pam_pos in rev_pam_sites:
        guide = extract_guide_sequence(rev_seq, pam_pos)
        
        if guide:
            pam = rev_seq[pam_pos:pam_pos+3]
            original_pos = len(sequence) - pam_pos - 3
            
            guides.append({
                'guide_sequence': guide,
                'pam_site': original_pos,
                'pam_sequence': pam,
                'strand': '-',
                'full_target': guide + pam
            })
    
    return pd.DataFrame(guides)


# Test it
if __name__ == "__main__":
    test_seq = "ATGCATGCATGCATGCATGCAGGCTAGCTAGCTAGCTAGC"
    
    print("Test Sequence:")
    print(test_seq)
    print(f"\nLength: {len(test_seq)} bp\n")
    
    pam_sites = find_pam_sites(test_seq)
    print(f"Found {len(pam_sites)} PAM sites: {pam_sites}\n")
    
    guides_df = find_all_guides(test_seq)
    print(f"Found {len(guides_df)} guide RNAs:\n")
    print(guides_df)