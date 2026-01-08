"""
Functions to predict off-target effects for CRISPR guide RNAs
"""

import pandas as pd
from Bio.Seq import Seq


def count_mismatches(seq1, seq2):
    """
    Count mismatches between two sequences.
    
    Args:
        seq1 (str): First sequence
        seq2 (str): Second sequence
    
    Returns:
        int: Number of mismatches
    
    Example:
        >>> count_mismatches("ATGC", "ATCC")
        1
    """
    if len(seq1) != len(seq2):
        return float('inf')  # Different lengths = not comparable
    
    mismatches = sum(1 for a, b in zip(seq1, seq2) if a != b)
    return mismatches


def find_similar_sequences(guide_sequence, target_sequence, max_mismatches=4):
    """
    Find all sequences in target that are similar to guide (potential off-targets).
    
    Off-target effects can occur when guide RNA is similar but not identical
    to other genomic locations. Typically, sites with ≤4 mismatches are concerning.
    
    Args:
        guide_sequence (str): 20bp guide RNA sequence
        target_sequence (str): Full sequence to search
        max_mismatches (int): Maximum mismatches to consider (default: 4)
    
    Returns:
        list: List of dicts with off-target information
    """
    guide_sequence = guide_sequence.upper()
    target_sequence = target_sequence.upper()
    guide_length = len(guide_sequence)
    
    off_targets = []
    
    # Search forward strand
    for i in range(len(target_sequence) - guide_length + 1):
        substring = target_sequence[i:i + guide_length]
        mismatches = count_mismatches(guide_sequence, substring)
        
        if mismatches <= max_mismatches and mismatches > 0:  # >0 to exclude perfect match
            off_targets.append({
                'position': i,
                'sequence': substring,
                'mismatches': mismatches,
                'strand': '+'
            })
    
    # Search reverse complement
    rev_target = str(Seq(target_sequence).reverse_complement())
    for i in range(len(rev_target) - guide_length + 1):
        substring = rev_target[i:i + guide_length]
        mismatches = count_mismatches(guide_sequence, substring)
        
        if mismatches <= max_mismatches and mismatches > 0:
            # Convert position back to forward strand coordinates
            original_pos = len(target_sequence) - i - guide_length
            off_targets.append({
                'position': original_pos,
                'sequence': substring,
                'mismatches': mismatches,
                'strand': '-'
            })
    
    return off_targets


def score_offtarget_risk(off_targets):
    """
    Calculate off-target risk score based on number and quality of off-targets.
    
    Scoring logic:
    - 0 mismatches: On-target (not counted)
    - 1 mismatch: Very high risk (50 points per site)
    - 2 mismatches: High risk (25 points per site)
    - 3 mismatches: Medium risk (10 points per site)
    - 4 mismatches: Low risk (5 points per site)
    
    Lower score = lower risk (better)
    
    Args:
        off_targets (list): List of off-target sites
    
    Returns:
        float: Risk score (0 = no off-targets, higher = more risk)
    """
    if not off_targets:
        return 0.0
    
    risk_score = 0.0
    
    # Weight off-targets by mismatch count
    mismatch_weights = {
        1: 50,   # 1 mismatch = very high risk
        2: 25,   # 2 mismatches = high risk
        3: 10,   # 3 mismatches = medium risk
        4: 5     # 4 mismatches = low risk
    }
    
    for ot in off_targets:
        weight = mismatch_weights.get(ot['mismatches'], 0)
        risk_score += weight
    
    return risk_score


def assess_offtarget_risk(guide_sequence, target_sequence, max_mismatches=4):
    """
    Complete off-target assessment for a guide RNA.
    
    Args:
        guide_sequence (str): Guide RNA sequence
        target_sequence (str): Target sequence to search
        max_mismatches (int): Maximum mismatches to consider
    
    Returns:
        dict: Off-target assessment with:
            - off_targets: List of potential off-target sites
            - risk_score: Numerical risk score
            - risk_level: 'Low', 'Medium', 'High', or 'Very High'
            - num_offtargets: Total number of off-targets found
    """
    # Find off-targets
    off_targets = find_similar_sequences(guide_sequence, target_sequence, max_mismatches)
    
    # Calculate risk score
    risk_score = score_offtarget_risk(off_targets)
    
    # Categorize risk level
    if risk_score == 0:
        risk_level = 'None'
    elif risk_score < 25:
        risk_level = 'Low'
    elif risk_score < 100:
        risk_level = 'Medium'
    elif risk_score < 200:
        risk_level = 'High'
    else:
        risk_level = 'Very High'
    
    return {
        'off_targets': off_targets,
        'risk_score': risk_score,
        'risk_level': risk_level,
        'num_offtargets': len(off_targets)
    }


def add_offtarget_scores(guides_df, target_sequence, max_mismatches=4):
    """
    Add off-target information to all guides in DataFrame.
    
    WARNING: This can be slow for large sequences!
    Consider using a subset of guides or smaller target sequence.
    
    Args:
        guides_df (pd.DataFrame): DataFrame with guide sequences
        target_sequence (str): Target sequence to search
        max_mismatches (int): Maximum mismatches to consider
    
    Returns:
        pd.DataFrame: DataFrame with added off-target columns
    """
    print(f"Analyzing off-targets for {len(guides_df)} guides...")
    print(f"Target sequence: {len(target_sequence):,} bp")
    print("This may take a few minutes for large sequences...\n")
    
    scored_df = guides_df.copy()
    
    # Initialize columns
    scored_df['num_offtargets'] = 0
    scored_df['offtarget_risk_score'] = 0.0
    scored_df['offtarget_risk_level'] = 'None'
    
    # Analyze each guide
    for idx, row in scored_df.iterrows():
        if idx % 100 == 0:
            print(f"  Processed {idx}/{len(scored_df)} guides...")
        
        assessment = assess_offtarget_risk(
            row['guide_sequence'],
            target_sequence,
            max_mismatches
        )
        
        scored_df.at[idx, 'num_offtargets'] = assessment['num_offtargets']
        scored_df.at[idx, 'offtarget_risk_score'] = assessment['risk_score']
        scored_df.at[idx, 'offtarget_risk_level'] = assessment['risk_level']
    
    print(f"✅ Completed off-target analysis!\n")
    
    return scored_df


def filter_by_offtarget_risk(guides_df, max_risk_level='Medium'):
    """
    Filter guides to keep only those with acceptable off-target risk.
    
    Args:
        guides_df (pd.DataFrame): Guides with off-target scores
        max_risk_level (str): Maximum acceptable risk ('Low', 'Medium', 'High')
    
    Returns:
        pd.DataFrame: Filtered guides
    """
    risk_order = {'None': 0, 'Low': 1, 'Medium': 2, 'High': 3, 'Very High': 4}
    max_risk_value = risk_order.get(max_risk_level, 2)
    
    filtered = guides_df[
        guides_df['offtarget_risk_level'].map(risk_order) <= max_risk_value
    ]
    
    print(f"Filtered: {len(filtered)}/{len(guides_df)} guides "
          f"with risk ≤ {max_risk_level}")
    
    return filtered


# Test functions
if __name__ == "__main__":
    print("Testing off-target detection...\n")
    
    # Test mismatch counting
    print("Test 1: Mismatch counting")
    print(f"ATGC vs ATGC: {count_mismatches('ATGC', 'ATGC')} mismatches (expected: 0)")
    print(f"ATGC vs ATCC: {count_mismatches('ATGC', 'ATCC')} mismatches (expected: 1)")
    print(f"ATGC vs TTCA: {count_mismatches('ATGC', 'TTCA')} mismatches (expected: 3)")
    
    # Test off-target finding
    print("\nTest 2: Off-target detection")
    guide = "ATGCATGCATGCATGCATGC"
    target = "ATGCATGCATGCATGCATGC" + "N"*100 + "ATGCATGCATGCATGCATCC"  # 1 mismatch
    
    off_targets = find_similar_sequences(guide, target, max_mismatches=2)
    print(f"Found {len(off_targets)} potential off-targets")
    
    for ot in off_targets:
        print(f"  Position {ot['position']}: {ot['mismatches']} mismatches")
    
    # Test risk scoring
    print("\nTest 3: Risk scoring")
    risk_score = score_offtarget_risk(off_targets)
    print(f"Risk score: {risk_score}")
    
    assessment = assess_offtarget_risk(guide, target, max_mismatches=3)
    print(f"Risk level: {assessment['risk_level']}")