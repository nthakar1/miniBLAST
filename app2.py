import streamlit as st
import pandas as pd
import time
import glob
import os
import io
import html
import dataclasses

from Bio import SeqIO

from main2 import miniBLASTn, miniBLASTp
from datatypes import BLASTN_PARAMS, BLASTP_PARAMS, scoring_presets
from blastStats import calculate_bit_score, calculate_e_value, compute_s1_threshold

# Initialize Session State
if "results_df" not in st.session_state:
    st.session_state.results_df = None

if "run_metadata" not in st.session_state:
    st.session_state.run_metadata = None

# -----------------------------
# ALIGNMENT VIEWER HELPERS
# -----------------------------
def alignment_events(aligned_query: str, aligned_ref: str):
    events = []
    for q, r in zip(aligned_query, aligned_ref):
        if q == r and q != "-":
            events.append("match")
        elif q == "-" and r != "-":
            events.append("gap_query")
        elif r == "-" and q != "-":
            events.append("gap_ref")
        else:
            events.append("mismatch")
    return events


def render_overview_bar(aligned_query: str, aligned_ref: str, width_px: int = 1250, height_px: int = 30):
    events = alignment_events(aligned_query, aligned_ref)
    n = len(events)

    if n == 0:
        st.info("No alignment available.")
        return

    spans = []
    for i, ev in enumerate(events):
        if ev == "match":
            color = "#d8d8d8"
        elif ev == "mismatch":
            color = "#e53935"
        elif ev == "gap_query":
            color = "#fb8c00"
        else:
            color = "#1e88e5"

        left = (i / n) * 100
        width = max(100 / n, 0.08)

        spans.append(
            f"<div style='position:absolute;left:{left}%;top:0;height:100%;width:{width}%;background:{color};'></div>"
        )

    html_block = f"""
    <div style="margin-top:6px;margin-bottom:14px;">
      <div style="font-size:14px;font-weight:600;margin-bottom:6px;">Alignment overview</div>
      <div style="
          position:relative;
          width:min(100%, {width_px}px);
          height:{height_px}px;
          border:1px solid #999;
          background:#f3f3f3;
          overflow:hidden;
          border-radius:4px;
      ">
        {''.join(spans)}
      </div>
      <div style="font-size:12px;color:#555;margin-top:5px;">
        Gray = match &nbsp;&nbsp; Red = mismatch &nbsp;&nbsp; Orange = gap in query &nbsp;&nbsp; Blue = gap in reference
      </div>
    </div>
    """
    st.markdown(html_block, unsafe_allow_html=True)


def render_alignment_window(aligned_query: str, aligned_ref: str, start: int = 0, window: int = 120):
    end = min(len(aligned_query), start + window)

    q_seg = aligned_query[start:end]
    r_seg = aligned_ref[start:end]

    q_html = []
    mid_html = []
    r_html = []

    for q, r in zip(q_seg, r_seg):
        q_disp = html.escape(q)
        r_disp = html.escape(r)

        if q == r and q != "-":
            q_html.append(f"<span>{q_disp}</span>")
            r_html.append(f"<span>{r_disp}</span>")
            mid_html.append("<span>|</span>")
        elif q == "-" and r != "-":
            q_html.append("<span style='background:#1e88e5;color:white;'>-</span>")
            r_html.append(f"<span style='background:#1e88e5;color:white;'>{r_disp}</span>")
            mid_html.append("<span style='color:#1e88e5;'> </span>")
        elif r == "-" and q != "-":
            q_html.append(f"<span style='background:#fb8c00;color:white;'>{q_disp}</span>")
            r_html.append("<span style='background:#fb8c00;color:white;'>-</span>")
            mid_html.append("<span style='color:#fb8c00;'> </span>")
        else:
            q_html.append(f"<span style='background:#e53935;color:white;'>{q_disp}</span>")
            r_html.append(f"<span style='background:#e53935;color:white;'>{r_disp}</span>")
            mid_html.append("<span style='color:#e53935;'>•</span>")

    ruler_numbers = "".join(str((start + i + 1) % 10) for i in range(end - start))

    # Use a fixed-width CSS label so every row's content starts at the same
    # horizontal offset.  Plain HTML collapses multiple spaces inside <b> tags,
    # which shifts the Ref / match rows relative to Query / Ruler – this avoids
    # that by giving each label an explicit inline-block width.
    lbl = "display:inline-block;width:7ch;font-weight:bold;flex-shrink:0;"

    html_block = f"""
    <div style="
        font-family:Menlo,Consolas,monospace;
        font-size:15px;
        line-height:1.7;
        background:#fafafa;
        border:1px solid #ddd;
        padding:12px;
        border-radius:8px;
        overflow-x:auto;
        white-space:nowrap;
    ">
      <div><span style="{lbl}">Pos</span>{start + 1}&nbsp;&ndash;&nbsp;{end}</div>
      <div><span style="{lbl}">Ruler</span>{ruler_numbers}</div>
      <div><span style="{lbl}">Query</span>{''.join(q_html)}</div>
      <div><span style="{lbl}">&nbsp;</span>{''.join(mid_html)}</div>
      <div><span style="{lbl}">Ref</span>{''.join(r_html)}</div>
    </div>
    """
    st.markdown(html_block, unsafe_allow_html=True)


def render_alignment_viewer(hit_row: pd.Series, key_prefix=""):
    aligned_query = hit_row["aligned_query"]
    aligned_ref = hit_row["aligned_ref"]

    st.subheader("Alignment Viewer")

    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Score", hit_row["score"])
    c2.metric("Bit score", hit_row["bit_score"])
    c3.metric("E-value", hit_row["e_value"])
    c4.metric("% identity", hit_row["pct_identity"])

    st.caption(
        f"Reference: {hit_row['ref_id']} | Query coverage: {hit_row['query_coverage']} | Position: {hit_row['position']}"
    )

    render_overview_bar(aligned_query, aligned_ref)

    window = st.slider("Zoom window size", min_value=40, max_value=300, value=120, step=10, key=f"{key_prefix}_zoom")
    max_start = max(0, len(aligned_query) - window)
    if max_start > 0:
        start = st.slider("Alignment start", min_value=0, max_value=max_start, value=0, step=1, key=f"{key_prefix}_start")
    else:
        start = 0

    render_alignment_window(aligned_query, aligned_ref, start=start, window=window)

# PAGE SETUP

tab1, tab2 = st.tabs(["miniBLASTn", "miniBLASTp"])

with tab1:
    st.title("miniBLASTn")
    st.write("Local alignment of nucleotide query sequence(s) against a database of DNA sequences.")

    # because sidebars exist globally, VS is removing sidebar method to ensure that each parameter specific to DNA or proteins stays in that tab
    k= st.number_input(
        "Word length",
        value=7,
        min_value=1,
        help="Initial seed length prior to extension",
        key="n_k"
    )

    scoring_parameters= st.selectbox(
        "Select a match reward and mismatch penalty scoring system for nucleotide sequence alignment. ",
        options= list(scoring_presets.keys()),
        key="n_scoring"
    )

    e_threshold= st.number_input(
        "Expected number of chance matches in a random model.",
        value=0.001,
        key="n_E",
        format="%.5g"
    )

    A = st.number_input(
        "Two-hit diagonal window (A)",
        value=40,
        min_value=1,
        help="Maximum diagonal distance between two hits to trigger extension. Default = 40.",
        key="n_A"
    )

    step = st.number_input(
        "Sample every Nth reference",
        value=1,
        min_value=1,
        help="1 = align against every sequence (slow), 100 = sample every 100th reference.",
        key="n_step"
    )

    # set the other parameters associated with the specific match reward and mismatch penalty the user chose
    selected_preset= scoring_presets[scoring_parameters]
    run_params=dataclasses.replace(BLASTN_PARAMS,**selected_preset)

    # -----------------------------
    # QUERY INPUT
    # -----------------------------
    st.header("Query Sequence")
    uploaded = st.file_uploader("Upload a FASTA file (.fna / .fasta / .fa):", type=["fna", "fasta", "fa"],key="n_qupload")
    raw_text = st.text_input("Or type a custom sequence:",key="n_qcustom")

    query = None

    if uploaded:
        records = list(SeqIO.parse(io.TextIOWrapper(uploaded, encoding="utf-8"), "fasta"))
        segments = [str(r.seq) for r in records]

        if len(segments) > 1:
            seg_index = st.selectbox(
                "Select query sequence",
                range(len(segments)),
                format_func=lambda i: f"Sequence {i + 1} ({len(segments[i])} bp)",
                key="n_seg"
            )
        else:
            seg_index = 0

        query = segments[seg_index]
        st.success(f"Loaded {len(segments)} sequence(s). Using sequence {seg_index + 1}.")
        st.code(query[:120] + ("..." if len(query) > 120 else ""))

    elif raw_text.strip():
        query = raw_text.strip().upper()
        st.success(f"Using typed sequence ({len(query)} bp).")
        st.code(query[:120] + ("..." if len(query) > 120 else ""))


    # -----------------------------
    # DATABASE SELECTION
    # -----------------------------
    st.header("Database")
    db_files = [
        f for f in glob.glob("*.fasta") + glob.glob("*.fna")
        if "database" in f.lower()
    ]

    if not db_files:
        st.error("No .fasta or .fna files containing 'database' were found in this directory.")
        st.stop()

    selected_db = st.selectbox("Select database", sorted(db_files),key="n_db")
    st.caption(f"Selected: `{selected_db}` ({os.path.getsize(selected_db) // 1_000_000} MB)")


    # -----------------------------
    # MAIN RUN
    # -----------------------------
    if st.button("Run miniBLAST", disabled=(query is None),key="n_run"):
        with st.spinner(f"Loading {selected_db}..."):
            db = list(SeqIO.parse(selected_db, "fasta"))
        st.success(f"Database loaded: {len(db)} sequences")

        effective_s1 = round(compute_s1_threshold(query, db, run_params, e_threshold), 2)
        st.info(f"s1 threshold (auto-computed): {effective_s1}  |  Mode: {"blastn"}")

        references = list(range(0, len(db), int(step)))
        results = []

        progress = st.progress(0, text="Aligning...")
        status = st.empty()

        t0 = time.perf_counter()

        for idx, i in enumerate(references):
            ref = str(db[i].seq)
            ref_id = db[i].id

            status.text(f"Aligning to reference {i} ({ref_id})")

            alignment = miniBLASTn(ref, query, s1=effective_s1, A=int(A))

            if alignment:
                raw_score = alignment["score"]
                bit_score = round(calculate_bit_score(raw_score, run_params), 4)
                e_value = float(f"{calculate_e_value(raw_score, len(query), len(db), run_params):.2e}")

                results.append({
                    "ref_index": i,
                    "ref_id": ref_id,
                    "score": round(raw_score, 2),
                    "bit_score": round(bit_score, 2),
                    "e_value": e_value,
                    "query_coverage": round(alignment["query_coverage"], 2),
                    "pct_identity": round(alignment["pct_identity"], 2),
                    "position": alignment["position"],
                    "aligned_query": alignment["alignment"][0],
                    "aligned_ref": alignment["alignment"][1],
                })

            progress.progress((idx + 1) / len(references), text=f"Aligned {idx + 1}/{len(references)}")

        elapsed = round(time.perf_counter() - t0, 2)
        status.empty()

        st.success(f"Done in {elapsed}s — {len(results)} alignments returned.")

        if results:
            df = pd.DataFrame(results)
            df = df.sort_values(by=["e_value", "bit_score", "score"], ascending=[True, False, False]).reset_index(drop=True)

            st.session_state.results_df = df
            st.session_state.run_metadata = {
                "elapsed": elapsed,
                "db_name": selected_db,
                "query_length": len(query),
                "num_results": len(results),
                "s1": effective_s1,
                "A": int(A),
                "step": int(step),
            }
        else:
            st.session_state.results_df_n = None
            st.session_state.run_metadata_n = None
            st.warning("No alignments passed the threshold.")


    # -----------------------------
    # DISPLAY SAVED RESULTS
    # -----------------------------
    if st.session_state.results_df is not None:
        df = st.session_state.results_df
        meta = st.session_state.run_metadata

        if meta is not None:
            st.success(
                f"Showing saved results from last run — DB: {meta['db_name']} | "
                f"Query length: {meta['query_length']} bp | "
                f"Results: {meta['num_results']} | "
                f"s1={meta['s1']}, A={meta['A']}, step={meta['step']} | "
                f"Run time: {meta['elapsed']}s"
            )

        st.header("Results")
        st.dataframe(
            df[[
                "ref_index",
                "ref_id",
                "score",
                "bit_score",
                "e_value",
                "query_coverage",
                "pct_identity",
                "position"
            ]],
            width="stretch",
            hide_index=True
        )

        best = df.iloc[0]
        st.subheader(f"Top hit: {best['ref_id']} (E-value: {best['e_value']})")

        chosen_ref = st.selectbox(
            "Choose hit to visualize",
            df["ref_id"].tolist(),
            index=0,
            key="n_chosen"
        )
        chosen_row = df[df["ref_id"] == chosen_ref].iloc[0]
        render_alignment_viewer(chosen_row, key_prefix="n")

        csv = df.drop(columns=["aligned_query", "aligned_ref"]).to_csv(index=False)
        st.download_button(
            "Download results CSV",
            csv,
            "blast_results.csv",
            "text/csv",
            key="n_csv"
        )
    

    
with tab2:
    st.title("miniBLASTp")
    st.write("Local alignment of protein query sequence(s) against a database of protein sequences.")

    # because sidebars exist globally, VS is removing sidebar method to ensure that each parameter specific to DNA or proteins stays in that tab
    k= st.number_input(
        "Word length",
        value=7,
        min_value=1,
        help="Initial seed length prior to extension",
        key="p_k"
    )
    
    threshold_T=st.number_input(
        "Add a numerical threshold value T that each protein alignment seed must reach prior to extension",
        value=11,
        min_value=0,
        key="p_T"
    )

    e_threshold= st.number_input(
        "Expected number of chance matches in a random model.",
        value=0.001,
        key="p_E",
        format="%.5g"
    )
    
    A = st.number_input(
        "Two-hit diagonal window (A)",
        value=40,
        min_value=1,
        help="Maximum diagonal distance between two hits to trigger extension. Default = 40.",
        key="p_A"
    )

    step = st.number_input(
        "Sample every Nth reference",
        value=1,
        min_value=1,
        help="1 = align against every sequence (slow), 100 = sample every 100th reference.",
        key="p_step"
    )

    # set the BLASTP T value based on user input
    run_params=dataclasses.replace(BLASTP_PARAMS,threshold_T=threshold_T)

    # -----------------------------
    # QUERY INPUT
    # -----------------------------
    st.header("Query Sequence")
    uploaded = st.file_uploader("Upload a FASTA file (.fna / .fasta / .fa):", type=["fna", "fasta", "fa"],key="p_qupload")
    raw_text = st.text_input("Or type a custom sequence:", key="p_qcustom")

    query = None

    if uploaded:
        records = list(SeqIO.parse(io.TextIOWrapper(uploaded, encoding="utf-8"), "fasta"))
        segments = [str(r.seq) for r in records]

        if len(segments) > 1:
            seg_index = st.selectbox(
                "Select query sequence",
                range(len(segments)),
                format_func=lambda i: f"Sequence {i + 1} ({len(segments[i])} bp)",
                key="p_seg"
            )
        else:
            seg_index = 0

        query = segments[seg_index]
        st.success(f"Loaded {len(segments)} sequence(s). Using sequence {seg_index + 1}.")
        st.code(query[:120] + ("..." if len(query) > 120 else ""))

    elif raw_text.strip():
        query = raw_text.strip().upper()
        st.success(f"Using typed sequence ({len(query)} bp).")
        st.code(query[:120] + ("..." if len(query) > 120 else ""))


    # -----------------------------
    # DATABASE SELECTION
    # -----------------------------
    st.header("Database")
    db_files = [
        f for f in glob.glob("*.fasta") + glob.glob("*.fna")
        if "database" in f.lower()
    ]

    if not db_files:
        st.error("No .fasta or .fna files containing 'database' were found in this directory.")
        st.stop()

    selected_db = st.selectbox("Select database", sorted(db_files), key="p_db")
    st.caption(f"Selected: `{selected_db}` ({os.path.getsize(selected_db) // 1_000_000} MB)")


    # -----------------------------
    # MAIN RUN
    # -----------------------------
    if st.button("Run miniBLAST", disabled=(query is None),key="p_run"):
        with st.spinner(f"Loading {selected_db}..."):
            db = list(SeqIO.parse(selected_db, "fasta"))
        st.success(f"Database loaded: {len(db)} sequences")

        effective_s1 = round(compute_s1_threshold(query, db, run_params), 2)
        st.info(f"s1 threshold (auto-computed): {effective_s1}  |  Mode: {"blastp"}")

        references = list(range(0, len(db), int(step)))
        results = []

        progress = st.progress(0, text="Aligning...")
        status = st.empty()

        t0 = time.perf_counter()

        for idx, i in enumerate(references):
            ref = str(db[i].seq)
            ref_id = db[i].id

            status.text(f"Aligning to reference {i} ({ref_id})")

            alignment = miniBLASTp(ref, query, s1=effective_s1, A=int(A),threshT=threshold_T)

            if alignment:
                raw_score = alignment["score"]
                bit_score = round(calculate_bit_score(raw_score, run_params), 4)
                e_value = float(f"{calculate_e_value(raw_score, len(query), len(db), run_params):.2e}")

                results.append({
                    "ref_index": i,
                    "ref_id": ref_id,
                    "score": round(raw_score, 2),
                    "bit_score": round(bit_score, 2),
                    "e_value": e_value,
                    "query_coverage": round(alignment["query_coverage"], 2),
                    "pct_identity": round(alignment["pct_identity"], 2),
                    "position": alignment["position"],
                    "aligned_query": alignment["alignment"][0],
                    "aligned_ref": alignment["alignment"][1],
                })

            progress.progress((idx + 1) / len(references), text=f"Aligned {idx + 1}/{len(references)}")

        elapsed = round(time.perf_counter() - t0, 2)
        status.empty()

        st.success(f"Done in {elapsed}s — {len(results)} alignments returned.")

        if results:
            df = pd.DataFrame(results)
            df = df.sort_values(by=["e_value", "bit_score", "score"], ascending=[True, False, False]).reset_index(drop=True)

            st.session_state.results_df = df
            st.session_state.run_metadata = {
                "elapsed": elapsed,
                "db_name": selected_db,
                "query_length": len(query),
                "num_results": len(results),
                "s1": effective_s1,
                "A": int(A),
                "step": int(step),
            }
        else:
            st.session_state.results_df = None
            st.session_state.run_metadata = None
            st.warning("No alignments passed the threshold.")


    # -----------------------------
    # DISPLAY SAVED RESULTS
    # -----------------------------
    if st.session_state.results_df is not None:
        df = st.session_state.results_df
        meta = st.session_state.run_metadata

        if meta is not None:
            st.success(
                f"Showing saved results from last run — DB: {meta['db_name']} | "
                f"Query length: {meta['query_length']} bp | "
                f"Results: {meta['num_results']} | "
                f"s1={meta['s1']}, A={meta['A']}, step={meta['step']} | "
                f"Run time: {meta['elapsed']}s"
            )

        st.header("Results")
        st.dataframe(
            df[[
                "ref_index",
                "ref_id",
                "score",
                "bit_score",
                "e_value",
                "query_coverage",
                "pct_identity",
                "position"
            ]],
            width="stretch",
            hide_index=True
        )

        best = df.iloc[0]
        st.subheader(f"Top hit: {best['ref_id']} (E-value: {best['e_value']})")

        chosen_ref = st.selectbox(
            "Choose hit to visualize",
            df["ref_id"].tolist(),
            index=0,
            key="p_chosen"
        )
        chosen_row = df[df["ref_id"] == chosen_ref].iloc[0]
        render_alignment_viewer(chosen_row, key_prefix="p")

        csv = df.drop(columns=["aligned_query", "aligned_ref"]).to_csv(index=False)
        st.download_button(
            "Download results CSV",
            csv,
            "blast_results.csv",
            "text/csv",
            key="p_csv"
        )
