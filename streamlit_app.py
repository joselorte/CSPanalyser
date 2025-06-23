import streamlit as st
import pandas as pd
from csp_utils import (
    run_csp_analysis,
    create_plotly_chart,
    create_matplotlib_plot
)

st.set_page_config(layout="wide")
st.title("🧪 Chemical Shift Perturbation (CSP) Analyzer")

with st.sidebar:
    st.header("📂 Upload Inputs")
    initial_file = st.file_uploader("Initial .xpk File", type="xpk")
    final_file = st.file_uploader("Final .xpk File", type="xpk")
    nucleus = st.selectbox("Select Nucleus Type", ["N", "C"])

if initial_file and final_file and nucleus:
    x_scale = 0.17 if nucleus == "N" else 0.251
    df = run_csp_analysis(initial_file, final_file, nucleus, x_scale)
    st.success("✅ CSP Analysis Complete")

    # --- FILTERING ---
    st.sidebar.header("🧪 Filter Options")
    max_res = int(df["Residue"].max())
    selected_range = st.slider("Residue Number Range", 1, max_res, (1, max_res))
    atom_list = sorted(df["Atom"].unique())
    selected_atoms = st.multiselect("Select Atoms", atom_list, default=atom_list)

    filtered_df = df[
        (df["Residue"] >= selected_range[0]) &
        (df["Residue"] <= selected_range[1]) &
        (df["Atom"].isin(selected_atoms))
    ]

    # --- THRESHOLDS BASED ON FULL DATA ---
    mean_csp = df["CSP"].mean()
    std_csp = df["CSP"].std()
    threshold1 = mean_csp + std_csp
    threshold2 = mean_csp + 2 * std_csp

    # --- ANNOTATION RESIDUE HIGHLIGHTING ---
    level1_df = filtered_df[filtered_df["CSP"] > threshold1]
    level2_df = filtered_df[filtered_df["CSP"] > threshold2]

    st.subheader("🔎 Selection Criteria Applied")
    st.markdown(f"""
    - **Residue range:** {selected_range[0]} to {selected_range[1]}
    - **Atoms included:** {', '.join(selected_atoms)}
    - **CSP thresholds (calculated using all CSP values):**
        - Level 1: > {threshold1:.3f}
        - Level 2: > {threshold2:.3f}
    - **Grey shaded regions:** Residues not found in either or both files
    """)

    st.markdown("🔬 **Highly Perturbed Residues:**")
    if not level1_df.empty:
        st.markdown(f"✅ *{len(level1_df)}* residues exceed CSP > mean + 1σ")
        st.dataframe(level1_df[["Residue", "Atom", "CSP"]], use_container_width=True)

        if not level2_df.empty:
            st.markdown(f"🔥 *{len(level2_df)}* residues exceed CSP > mean + 2σ")
            st.dataframe(level2_df[["Residue", "Atom", "CSP"]], use_container_width=True)
        else:
            st.info("No residues exceed CSP > mean + 2σ")
    else:
        st.info("No residues exceed CSP > mean + 1σ")

    # --- INTERACTIVE PLOT ---
    st.subheader("📊 CSP Plot (Interactive)")
    fig = create_plotly_chart(
        filtered_df,
        nucleus,
        thresholds=(threshold1, threshold2),
        res_range=selected_range
    )
    st.plotly_chart(fig, use_container_width=True)

    with st.expander("🧠 What Does This Plot Show?", expanded=False):
        st.markdown("""
        This bar chart visualizes **chemical shift perturbation (CSP)** values across selected residues and atoms.

        - **X-axis:** Residue number (based on your slider).
        - **Y-axis:** CSP = √((ΔH)² + α×(ΔX)²), where α = 0.17 (¹⁵N) or 0.251 (¹³C)
        - **Bar colors:** Atom names (¹³C) or uniform (¹⁵N)
        - **Red dashed line:** mean + 1σ
        - **Red dotted line:** mean + 2σ
        - **Grey shading:** Residues not assigned or missing in input files
        """)

    # --- HIGH-RES DOWNLOADABLE PLOTS ---
    st.subheader("🖼️ Download High-Resolution Plots")
    img_buf_clean = create_matplotlib_plot(
        filtered_df,
        nucleus,
        annotate=False,
        res_range=selected_range,
        thresholds=(threshold1, threshold2)
    )
    img_buf_annotated = create_matplotlib_plot(
        filtered_df,
        nucleus,
        annotate=True,
        res_range=selected_range,
        thresholds=(threshold1, threshold2)
    )

    col1, col2 = st.columns(2)
    with col1:
        st.download_button(
            "📥 Download Clean Plot (PNG)",
            data=img_buf_clean,
            file_name="csp_plot_clean.png",
            mime="image/png"
        )
    with col2:
        st.download_button(
            "📥 Download Annotated Plot (PNG)",
            data=img_buf_annotated,
            file_name="csp_plot_annotated.png",
            mime="image/png"
        )

    # --- DATA TABLE + CSV EXPORT ---
    st.subheader("📋 CSP Table")
    st.dataframe(filtered_df, use_container_width=True)
    csv = filtered_df.to_csv(index=False).encode("utf-8")
    st.download_button(
        "📥 Download Filtered CSV",
        data=csv,
        file_name="filtered_csp.csv",
        mime="text/csv"
    )
