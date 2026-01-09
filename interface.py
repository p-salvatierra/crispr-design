"""
CRISPR Guide RNA Design Tool - Web Interface
Interactive tool for designing CRISPR-Cas9 guide RNAs
"""

import streamlit as st
import pandas as pd
import sys
from io import StringIO
from Bio import SeqIO
from Bio.Seq import Seq

# Add src to path
sys.path.append('src')

from find_guides import find_all_guides
from score_guides import score_all_guides, get_top_guides, visualize_guide_scores
from offtarget_prediction import add_offtarget_scores, filter_by_offtarget_risk

# Page config
st.set_page_config(
    page_title="CRISPR Guide Designer",
    page_icon="🧬",
    layout="wide"
)

# Title and description
st.title("🧬 CRISPR Guide RNA Design Tool")
st.markdown("""
Design optimal CRISPR-Cas9 guide RNAs with efficiency scoring and off-target prediction.

**Features:**
- Find all PAM sites (NGG) in your sequence
- Score guides by GC content, poly-T, and position
- Predict off-target effects
- Export results to CSV
""")

st.divider()

# Sidebar - Input options
st.sidebar.header("Input Options")

input_method = st.sidebar.radio(
    "Choose input method:",
    ["Paste Sequence", "Upload FASTA File", "Use Example"]
)

sequence = None

# Handle different input methods
if input_method == "Paste Sequence":
    st.sidebar.info("Paste your DNA sequence (FASTA format or plain text)")
    sequence_input = st.sidebar.text_area(
        "DNA Sequence:",
        height=200,
        placeholder="Paste sequence here..."
    )
    
    if sequence_input:
        # Remove FASTA header if present
        if sequence_input.startswith(">"):
            lines = sequence_input.split("\n")
            sequence = "".join(lines[1:]).replace(" ", "").replace("\n", "")
        else:
            sequence = sequence_input.replace(" ", "").replace("\n", "")

elif input_method == "Upload FASTA File":
    uploaded_file = st.sidebar.file_uploader("Choose a FASTA file", type=['fasta', 'fa', 'txt'])
    
    if uploaded_file is not None:
        stringio = StringIO(uploaded_file.getvalue().decode("utf-8"))
        record = SeqIO.read(stringio, "fasta")
        sequence = str(record.seq)
        st.sidebar.success(f"✅ Loaded: {record.id} ({len(sequence):,} bp)")

else:  # Use Example
    st.sidebar.info("Using example BRCA1 sequence (first 2000 bp)")
    # Example sequence - first part of BRCA1
    sequence = "ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAAATCTTAGAGTGTCCCATCTGTCTGGAGTTGATCAAGGAACCTGTCTCCACAAAGTGTGACCACATATTTTGCAAATTTTGCATGCTGAAACTTCTCAACCAGAAGAAAGGGCCTTCACAGTGTCCTTTATGTAAGAATGATATAACCAAAAGGAGCCTACAAGAAAGTACGAGATTTAGTCAACTTGTTGAAGAGCTATTGAAAATCATTTGTGCTTTTCAGCTTGACACAGGTTTGGAGTATGCAAACAGCTATAATTTTGCAAAAAAGGAAAATAACTCTCCTGAACATCTAAAAGATGAAGTTTCTATCATCCAAAGTATGGGCTACAGAAACCGTGCCAAAAGACTTCTACAGAGTGAACCCGAAAATCCTTCCTTGCAGGAAACCAGTCTCAGTGTCCAACTCTCTAACCTTGGAACTGTGAGAACTCTGAGGACAAAGCAGCGGATACAACCTCAAAAGACGTCTGTCTACATTGAATTGGGATCTGATTCTTCTGAAGATACCGTTAATAAGGCAACTTATTGCAGTGTGGGAGATCAAGAATTGTTACAAATCACCCCTCAAGGAACCAGGGATGAAATCAGTTTGGATTCTGCAAAAAAGGCTGCTTGTGAATTTTCTGAGACGGATGTAACAAATACTGAACATCATCAACCCAGTAATAATGATTTGAACACCACTGAGAAGCGTGCAGCTGAGAGGCATCCAGAAAAGTATCAGGGTAGTTCTGTTTCAAACTTGCATGTGGAGCCATGTGGCACAAATACTCATGCCAGCTCATTACAGCATGAGAACAGCAGTTTATTACTCACTAAAGACAGAATGAATGTAGAAAAGGCTGAATTCTGTAATAAAAGCAAACAGCCTGGCTTAGCAAGGAGCCAACATAACAGATGGGCTGGAAGTAAGGAAACATGTAATGATAGGCGGACTCCCAGCACAGAAAAAAAGGTAGATCTGAATGCTGATCCCCTGTGTGAGAGAAAAGAATGGAATA"
    
    st.sidebar.success(f"✅ Loaded example ({len(sequence):,} bp)")

# Analysis parameters
st.sidebar.header("Analysis Parameters")

min_efficiency_score = st.sidebar.slider(
    "Minimum Efficiency Score:",
    min_value=0,
    max_value=100,
    value=50,
    help="Filter guides below this efficiency score"
)

include_offtarget = st.sidebar.checkbox(
    "Include Off-Target Analysis",
    value=False,
    help="⚠️ Slower for large sequences"
)

if include_offtarget:
    max_mismatches = st.sidebar.slider(
        "Max Mismatches:",
        min_value=1,
        max_value=4,
        value=3,
        help="Maximum mismatches to consider for off-targets"
    )
    
    max_risk_level = st.sidebar.selectbox(
        "Maximum Risk Level:",
        ["Low", "Medium", "High"],
        index=1
    )

num_guides_to_show = st.sidebar.slider(
    "Number of Top Guides:",
    min_value=5,
    max_value=50,
    value=10
)

# Run analysis button
run_analysis = st.sidebar.button("🚀 Design Guide RNAs", type="primary", use_container_width=True)

# Main content
if sequence and run_analysis:
    # Validate sequence
    sequence = sequence.upper()
    valid_bases = set('ATGCN')
    if not set(sequence).issubset(valid_bases):
        st.error("❌ Invalid sequence! Only A, T, G, C, N characters allowed.")
        st.stop()
    
    if len(sequence) < 23:
        st.error("❌ Sequence too short! Minimum 23 bp required.")
        st.stop()
    
    # Show sequence info
    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("Sequence Length", f"{len(sequence):,} bp")
    with col2:
        gc = (sequence.count('G') + sequence.count('C')) / len(sequence) * 100
        st.metric("GC Content", f"{gc:.1f}%")
    with col3:
        st.metric("AT Content", f"{100-gc:.1f}%")
    
    st.divider()
    
    # Step 1: Find guides
    with st.spinner("🔍 Finding PAM sites and extracting guide RNAs..."):
        guides_df = find_all_guides(sequence)
    
    st.success(f"✅ Found {len(guides_df)} potential guide RNAs")
    
    # Step 2: Score guides
    with st.spinner("📊 Scoring guide efficiency..."):
        scored_guides = score_all_guides(guides_df, sequence_length=len(sequence))
    
    # Step 3: Off-target analysis (optional)
    if include_offtarget:
        # Only analyze top guides to save time
        top_for_offtarget = scored_guides.head(50)
        
        with st.spinner(f"🎯 Analyzing off-targets (this may take 1-2 minutes)..."):
            scored_guides_subset = add_offtarget_scores(
                top_for_offtarget,
                sequence,
                max_mismatches=max_mismatches
            )
            
            # Filter by risk
            scored_guides_subset = filter_by_offtarget_risk(
                scored_guides_subset,
                max_risk_level=max_risk_level
            )
        
        # Use the filtered subset
        final_guides = scored_guides_subset
    else:
        final_guides = scored_guides
    
    # Filter by minimum score
    final_guides = final_guides[final_guides['efficiency_score'] >= min_efficiency_score]
    
    # Get top N
    top_guides = final_guides.head(num_guides_to_show)
    
    if len(top_guides) == 0:
        st.warning("⚠️ No guides meet the criteria. Try lowering the minimum score or risk level.")
        st.stop()
    
    st.success(f"✅ {len(top_guides)} high-quality guides selected")
    
    st.divider()
    
    # Display results
    st.header("📋 Top Guide RNAs")
    
    # Format dataframe for display
    display_df = top_guides.copy()
    display_df['rank'] = range(1, len(display_df) + 1)
    
    columns_to_show = ['rank', 'guide_sequence', 'pam_sequence', 'efficiency_score', 
                       'gc_content', 'has_poly_t', 'strand', 'pam_site']
    
    if include_offtarget:
        columns_to_show.extend(['num_offtargets', 'offtarget_risk_level', 'offtarget_risk_score'])
    
    st.dataframe(
        display_df[columns_to_show],
        use_container_width=True,
        hide_index=True
    )
    
    # Detailed view of top 5
    st.subheader("🔬 Detailed View - Top 5 Guides")
    
    for idx, row in top_guides.head(5).iterrows():
        with st.expander(f"**Guide #{row.get('rank', idx+1)}** - Score: {row['efficiency_score']:.1f}"):
            col1, col2 = st.columns(2)
            
            with col1:
                st.markdown("**Sequence Information:**")
                st.code(f"Guide:  5'-{row['guide_sequence']}-3'\nPAM:       {row['pam_sequence']}")
                st.markdown(f"**Position:** {row['pam_site']:,} bp ({row['strand']} strand)")
                st.markdown(f"**Full Target:** `{row['full_target']}`")
            
            with col2:
                st.markdown("**Quality Metrics:**")
                st.markdown(f"**Efficiency Score:** {row['efficiency_score']:.1f}/100")
                st.markdown(f"**GC Content:** {row['gc_content']:.1f}%")
                st.markdown(f"**Poly-T Present:** {'⚠️ Yes' if row['has_poly_t'] else '✅ No'}")
                
                if include_offtarget:
                    risk_color = {
                        'None': '🟢',
                        'Low': '🟢', 
                        'Medium': '🟡',
                        'High': '🟠',
                        'Very High': '🔴'
                    }
                    st.markdown(f"**Off-Target Risk:** {risk_color.get(row['offtarget_risk_level'], '⚪')} {row['offtarget_risk_level']}")
                    st.markdown(f"**Off-Target Sites:** {row['num_offtargets']}")
    
    # Visualization
    st.divider()
    st.header("📊 Analysis Visualizations")
    
    fig = visualize_guide_scores(final_guides, top_n=min(50, len(final_guides)))
    st.pyplot(fig)
    
    # Export options
    st.divider()
    st.header("💾 Export Results")
    
    col1, col2 = st.columns(2)
    
    with col1:
        # CSV export
        csv = top_guides.to_csv(index=False)
        st.download_button(
            label="📥 Download Top Guides (CSV)",
            data=csv,
            file_name="crispr_guides.csv",
            mime="text/csv"
        )
    
    with col2:
        # Summary report
        report = f"""CRISPR Guide RNA Design Report
================================

Sequence Information:
- Length: {len(sequence):,} bp
- GC Content: {gc:.1f}%

Analysis Parameters:
- Minimum Efficiency Score: {min_efficiency_score}
- Off-Target Analysis: {'Yes' if include_offtarget else 'No'}

Results:
- Total Guides Found: {len(guides_df)}
- High-Quality Guides: {len(top_guides)}

Top 5 Recommended Guides:
"""
        for idx, row in top_guides.head(5).iterrows():
            report += f"\n{idx+1}. {row['guide_sequence']} (Score: {row['efficiency_score']:.1f})"
        
        st.download_button(
            label="📄 Download Summary Report",
            data=report,
            file_name="crispr_design_report.txt",
            mime="text/plain"
        )

elif not sequence:
    # Welcome message
    st.info("👈 Choose an input method in the sidebar to get started!")
    
    st.markdown("""
    ### How to Use:
    
    1. **Input your sequence** using one of three methods:
       - Paste DNA sequence directly
       - Upload a FASTA file
       - Use the example sequence
    
    2. **Adjust parameters** (optional):
       - Set minimum efficiency score
       - Enable off-target analysis
       - Choose number of guides to display
    
    3. **Click "Design Guide RNAs"** to run the analysis
    
    4. **Review results** and download your guides!
    
    ### About This Tool:
    
    This tool helps design CRISPR-Cas9 guide RNAs by:
    - Finding all NGG PAM sites in your sequence
    - Scoring guides based on GC content, poly-T sequences, and position
    - Optionally predicting off-target effects
    - Providing downloadable results
    
    **Developed by:** Paula Salvatierra  
    **GitHub:** https://github.com/p-salvatierra/crispr-design
    """)

else:
    st.info("👆 Click the button in the sidebar to start analysis!")