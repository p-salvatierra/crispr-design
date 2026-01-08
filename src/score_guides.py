"""
Functions to score CRISPR guide RNA efficiency
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def calculate_gc_content(sequence):
    """
    Calculate GC content (percentage of G and C bases).
    
    Optimal GC content for guide RNAs is typically 40-60%.
    Too low or too high can affect cutting efficiency.
    
    Args:
        sequence (str): DNA sequence
    
    Returns:
        float: GC content as percentage (0-100)
    
    Example:
        >>> calculate_gc_content("ATGCATGC")
        50.0
    """
    sequence = sequence.upper()
    
    if len(sequence) == 0:
        return 0.0
    
    gc_count = sequence.count('G') + sequence.count('C')
    gc_content = (gc_count / len(sequence)) * 100
    
    return round(gc_content, 2)


def has_poly_t(sequence, threshold=4):
    """
    Check if sequence contains poly-T stretch (4+ T's in a row).
    
    Poly-T sequences (TTTT) can cause premature transcription termination,
    making them unfavorable for guide RNAs.
    
    Args:
        sequence (str): DNA sequence
        threshold (int): Number of consecutive T's to flag (default: 4)
    
    Returns:
        bool: True if poly-T found, False otherwise
    
    Example:
        >>> has_poly_t("ATGCTTTTGC")
        True
        >>> has_poly_t("ATGCTTTGC")
        False
    """
    sequence = sequence.upper()
    poly_t = 'T' * threshold
    
    return poly_t in sequence


def score_gc_content(gc_content):
    """
    Score GC content on scale of 0-100.
    
    Optimal range: 40-60% (score = 100)
    Acceptable range: 30-70% (score = 50-100)
    Outside range: Lower scores
    
    Args:
        gc_content (float): GC percentage
    
    Returns:
        float: Score from 0-100
    """
    if 40 <= gc_content <= 60:
        # Optimal range
        return 100.0
    elif 30 <= gc_content < 40:
        # Below optimal - linear penalty
        return 50 + (gc_content - 30) * 5
    elif 60 < gc_content <= 70:
        # Above optimal - linear penalty
        return 50 + (70 - gc_content) * 5
    else:
        # Outside acceptable range - harsh penalty
        if gc_content < 30:
            return max(0, gc_content * 1.67)  # 0-50 score
        else:  # gc_content > 70
            return max(0, (100 - gc_content) * 1.67)  # 0-50 score


def calculate_position_score(pam_position, sequence_length):
    """
    Score based on guide position in the gene.
    
    Guides closer to the 5' end (start) of a gene are often preferred
    for knockout experiments.
    
    Args:
        pam_position (int): Position of PAM site
        sequence_length (int): Total length of sequence
    
    Returns:
        float: Score from 0-100
    """
    # Calculate relative position (0-1)
    relative_position = pam_position / sequence_length
    
    # Prefer guides in first 50% of gene
    if relative_position <= 0.5:
        return 100.0
    else:
        # Linear decrease from 100 to 50 as you go from 50% to 100%
        return 100 - (relative_position - 0.5) * 100


def calculate_efficiency_score(guide_sequence, pam_position=None, sequence_length=None):
    """
    Calculate overall efficiency score for a guide RNA.
    
    Combines multiple factors:
    - GC content (40% weight)
    - Poly-T presence (30% penalty if present)
    - Position in gene (30% weight, if provided)
    
    Args:
        guide_sequence (str): 20bp guide sequence
        pam_position (int, optional): Position of PAM in gene
        sequence_length (int, optional): Total length of target sequence
    
    Returns:
        float: Overall efficiency score (0-100)
    """
    # Calculate GC content score
    gc_content = calculate_gc_content(guide_sequence)
    gc_score = score_gc_content(gc_content)
    
    # Check for poly-T
    has_polyt = has_poly_t(guide_sequence)
    polyt_penalty = 30 if has_polyt else 0
    
    # Calculate position score if position provided
    if pam_position is not None and sequence_length is not None:
        position_score = calculate_position_score(pam_position, sequence_length)
        # Weighted average: 40% GC, 30% position
        overall_score = (gc_score * 0.4) + (position_score * 0.3) + 30
    else:
        # Just use GC score if no position info
        overall_score = gc_score
    
    # Apply poly-T penalty
    overall_score = max(0, overall_score - polyt_penalty)
    
    return round(overall_score, 2)


def score_all_guides(guides_df, sequence_length=None):
    """
    Score all guides in a DataFrame.
    
    Args:
        guides_df (pd.DataFrame): DataFrame from find_all_guides()
        sequence_length (int, optional): Length of target sequence
    
    Returns:
        pd.DataFrame: Original DataFrame with added scoring columns
    """
    # Make a copy to avoid modifying original
    scored_df = guides_df.copy()
    
    # Calculate GC content for each guide
    scored_df['gc_content'] = scored_df['guide_sequence'].apply(calculate_gc_content)
    
    # Check for poly-T
    scored_df['has_poly_t'] = scored_df['guide_sequence'].apply(has_poly_t)
    
    # Calculate efficiency score
    if sequence_length:
        scored_df['efficiency_score'] = scored_df.apply(
            lambda row: calculate_efficiency_score(
                row['guide_sequence'],
                row['pam_site'],
                sequence_length
            ),
            axis=1
        )
    else:
        scored_df['efficiency_score'] = scored_df['guide_sequence'].apply(
            lambda seq: calculate_efficiency_score(seq)
        )
    
    # Sort by efficiency score (best first)
    scored_df = scored_df.sort_values('efficiency_score', ascending=False)
    
    # Add rank
    scored_df['rank'] = range(1, len(scored_df) + 1)
    
    return scored_df


def get_top_guides(scored_df, n=10, min_score=50):
    """
    Get top N guides above minimum score threshold.
    
    Args:
        scored_df (pd.DataFrame): Scored guides DataFrame
        n (int): Number of top guides to return (default: 10)
        min_score (float): Minimum efficiency score (default: 50)
    
    Returns:
        pd.DataFrame: Top guides meeting criteria
    """
    # Filter by minimum score
    filtered = scored_df[scored_df['efficiency_score'] >= min_score]
    
    # Return top N
    return filtered.head(n)


def visualize_guide_scores(scored_df, top_n=20):
    """
    Create visualizations of guide RNA scores.
    
    Args:
        scored_df (pd.DataFrame): Scored guides
        top_n (int): Number of top guides to show
    """
    # Get top guides
    top_guides = scored_df.head(top_n)
    
    # Create figure with subplots
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    # 1. Efficiency scores bar plot
    ax1 = axes[0, 0]
    ax1.barh(range(len(top_guides)), top_guides['efficiency_score'])
    ax1.set_yticks(range(len(top_guides)))
    ax1.set_yticklabels([f"Guide {i+1}" for i in range(len(top_guides))])
    ax1.set_xlabel('Efficiency Score')
    ax1.set_title(f'Top {top_n} Guide RNA Efficiency Scores')
    ax1.axvline(x=50, color='red', linestyle='--', alpha=0.5, label='Min threshold')
    ax1.legend()
    ax1.invert_yaxis()
    
    # 2. GC content distribution
    ax2 = axes[0, 1]
    ax2.hist(scored_df['gc_content'], bins=20, edgecolor='black', alpha=0.7)
    ax2.axvline(x=40, color='green', linestyle='--', alpha=0.5, label='Optimal range')
    ax2.axvline(x=60, color='green', linestyle='--', alpha=0.5)
    ax2.set_xlabel('GC Content (%)')
    ax2.set_ylabel('Number of Guides')
    ax2.set_title('GC Content Distribution')
    ax2.legend()
    
    # 3. Score vs GC content scatter
    ax3 = axes[1, 0]
    colors = ['red' if pt else 'blue' for pt in scored_df['has_poly_t']]
    ax3.scatter(scored_df['gc_content'], scored_df['efficiency_score'], 
                c=colors, alpha=0.6, s=50)
    ax3.set_xlabel('GC Content (%)')
    ax3.set_ylabel('Efficiency Score')
    ax3.set_title('Efficiency Score vs GC Content')
    ax3.axhline(y=50, color='gray', linestyle='--', alpha=0.3)
    ax3.axvline(x=40, color='gray', linestyle='--', alpha=0.3)
    ax3.axvline(x=60, color='gray', linestyle='--', alpha=0.3)
    
    # Legend for colors
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor='blue', label='No Poly-T'),
                      Patch(facecolor='red', label='Has Poly-T')]
    ax3.legend(handles=legend_elements)
    
    # 4. Position distribution
    ax4 = axes[1, 1]
    forward = scored_df[scored_df['strand'] == '+']
    reverse = scored_df[scored_df['strand'] == '-']
    
    ax4.scatter(forward['pam_site'], forward['efficiency_score'], 
                label='Forward (+)', alpha=0.6, s=50)
    ax4.scatter(reverse['pam_site'], reverse['efficiency_score'], 
                label='Reverse (-)', alpha=0.6, s=50)
    ax4.set_xlabel('Position in Sequence')
    ax4.set_ylabel('Efficiency Score')
    ax4.set_title('Guide Score vs Position')
    ax4.legend()
    
    plt.tight_layout()
    return fig

# Test the functions
if __name__ == "__main__":
    # Test individual functions
    print("Testing GC Content:")
    print(f"ATGC: {calculate_gc_content('ATGC')}% (expected: 50%)")
    print(f"GGGG: {calculate_gc_content('GGGG')}% (expected: 100%)")
    print(f"AAAA: {calculate_gc_content('AAAA')}% (expected: 0%)")
    
    print("\nTesting Poly-T Detection:")
    print(f"ATGCTTTTGC: {has_poly_t('ATGCTTTTGC')} (expected: True)")
    print(f"ATGCTTTGC: {has_poly_t('ATGCTTTGC')} (expected: False)")
    
    print("\nTesting GC Scoring:")
    print(f"50% GC: {score_gc_content(50)} (expected: 100)")
    print(f"35% GC: {score_gc_content(35)} (expected: 75)")
    print(f"65% GC: {score_gc_content(65)} (expected: 75)")
    
    print("\nTesting Efficiency Score:")
    good_guide = "ATGCATGCATGCATGCATGC"  # 50% GC, no poly-T
    bad_guide = "AAAATTTTAAAATTTTAAAA"   # 0% GC, has poly-T
    
    print(f"Good guide: {calculate_efficiency_score(good_guide)}")
    print(f"Bad guide: {calculate_efficiency_score(bad_guide)}")