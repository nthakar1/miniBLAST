import streamlit as st
import pandas as pd
import time
import glob
import os
import io
from Bio import SeqIO
from main import miniBLASTn
from database import fetch_mixed_database
from datatypes import BlastConfig, BLASTN_PARAMS, BLASTP_PARAMS, BLOSUM
from statistics import calculate_bit_score, calculate_e_value
import html

st.title("miniBLAST")
st.write("Local alignment of nucleotide or protein sequences against a database.")
#TODO: create separate tabs for miniBLASTn and miniBLASTp
# Configurable parameters that users can change. Help text will appear as a question mark tooltip.
st.sidebar.header("Parameters")
s1 = st.sidebar.number_input("Ungapped extension threshold (s1)", value=20, min_value=1,
                                help= "Minimum ungapped extension score needed to proceed to gapped alignment. Default value= 20.")
A  = st.sidebar.number_input("Two-hit diagonal window (A)",        value=40, min_value=1,
                                help="Maximum diagonal distance between two hits to trigger extension. Default value= 40.")
step = st.sidebar.number_input("Sample every Nth reference",        value=1, min_value=1,
                                help="1 = align against every sequence (slow), 100 = sample every 100th")

# Upload a FASTA query sequence with file type= fna, fasta, or fa
st.header("Query Sequence")
uploaded = st.file_uploader("Upload a FASTA file (.fna / .fasta):", type=["fna", "fasta", "fa"])
raw_text=st.text_input(label= "Or type a custom sequence:")

# Once a user uploads a FASTA file or types a custom sequence, check how many sequences are in the fasta file. If there is more than 1 sequence in the FASTA file, allow them to choose which sequence will be their query sequence.
query = None
if uploaded:
    segments = [str(r.seq) for r in SeqIO.parse(io.TextIOWrapper(uploaded, encoding="utf-8"), "fasta")]
    if len(segments) > 1:
        seg_index = st.selectbox(
            "Select query sequence",
            range(len(segments)),
            format_func=lambda i: f"Sequence {i+1} ({len(segments[i])} bp)"
        )
    else:
        seg_index = 0
    query = segments[seg_index]
    st.success(f"Loaded {len(segments)} sequence(s). Using sequence {seg_index+1}.")
    st.code(query[:80] + ("..." if len(query) > 80 else ""), language=None)
elif raw_text.strip():
    query = raw_text.strip()
    st.success(f"Using typed sequence ({len(query)} bp).")
    st.code(query[:80] + ("..." if len(query) > 80 else ""), language=None)


# Configuring database with a dropdown menu that pulls FASTA/fna files from our repository with database in the filename. Users can select which database they want to align against.
st.header("Database") 
# FOR LATER TODO: User builds custom database using database.py 
db_files = [f for f in glob.glob("*.fasta") + glob.glob("*.fna")
            if "database" in f.lower()]


if not db_files:
    st.error("No .fasta or .fna files found in the project directory.")
    st.stop()

selected_db = st.selectbox("Select database", sorted(db_files))
st.caption(f"Selected: `{selected_db}`  ({os.path.getsize(selected_db) // 1_000_000} MB)")

# Streamlit app can run only once query sequence(s) are uploaded. SeqIO will read the fasta files from the selected database.
if st.button("Run miniBLAST", disabled=(query is None)): 
    with st.spinner(f"Loading {selected_db}..."):
        db = list(SeqIO.parse(selected_db, "fasta"))
    st.success(f"Database loaded: {len(db)} sequences")

    
    references = list(range(0, len(db), int(step)))
    results = []
    progress = st.progress(0, text="Aligning...")
    status   = st.empty()

    # begin alignment of query sequence to reference sequence from database
    t0 = time.perf_counter()
    for idx, i in enumerate(references):
        ref = str(db[i].seq)
        ref_id = db[i].id
        status.text(f"Aligning to reference {i} ({ref_id})")

        alignment = miniBLASTn(ref, query, s1=int(s1), A=int(A))

        if alignment:
            alignment["bit_score"] = round(calculate_bit_score(alignment["score"], BLASTN_PARAMS), 4)
            alignment["e_value"]   = round(calculate_e_value(alignment["score"], len(query), len(db), BLASTN_PARAMS), 4)
            results.append({
                "ref_index":      i,
                "ref_id":         ref_id,
                "score":          round(alignment["score"], 2),
                "bit_score":      round(alignment["bit_score"],2),
                "e_value": float(f"{calculate_e_value(alignment['score'], len(query), len(db), BLASTN_PARAMS):.2e}"),
                "query_coverage": round(alignment["query_coverage"], 2),
                "pct_identity":   round(alignment["pct_identity"], 2),
                "position":       alignment["position"],
                "aligned_query":  alignment["alignment"][0],
                "aligned_ref":    alignment["alignment"][1],
            })
        progress.progress((idx + 1) / len(references),
                          text=f"Aligned {idx+1}/{len(references)}")

    elapsed = round(time.perf_counter() - t0, 2)
    status.empty()
    st.success(f"Done in {elapsed}s — {len(results)} alignments returned.")

    if results:
        df = pd.DataFrame(results)

        st.header("Results")
        st.dataframe(
            df[["ref_index", "ref_id", "score", "bit_score","e_value", "query_coverage", "pct_identity", "position"]],
            use_container_width=True
        )

        # Top hit information
        best = df.loc[df["e_value"].idxmin()]
        st.subheader(f"Top hit: {best['ref_id']}  (E_value: {best['e_value']})")
        st.text("Query: " + best["aligned_query"])
        st.text("Ref:   " + best["aligned_ref"])

        # Download results file
        csv = df.drop(columns=["aligned_query", "aligned_ref"]).to_csv(index=False)
        st.download_button("Download results CSV", csv, "blast_results.csv", "text/csv")
    else:
        st.warning("No alignments passed the threshold.")

