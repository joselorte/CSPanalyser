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

    # --- ANNOTATION ---
    threshold1 = df["CSP"].mean() + df["CSP"].std()
    threshold2 = df["CSP"].mean() + 2 * df["CSP"].std()
    level1_df = filtered_df[filtered_df["CSP"] > threshold1]
    level2_df = filtered_df[filtered_df["CSP"] > threshold2]

    st.subheader("🔎 Selection Criteria Applied")
    st.markdown(f"""
    - **Residue range:** {selected_range[0]} to {selected_range[1]}
    - **Atoms included:** {', '.join(selected_atoms)}
    - **CSP thresholds:**
        - Level 1: CSP > {threshold1:.3f}
        - Level 2: CSP > {threshold2:.3f}
    - **Grey shaded regions:** Unassigned residues and residue numbers without CSP data
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
    fig = create_plotly_chart(filtered_df, nucleus, thresholds=(threshold1, threshold2), res_range=selected_range)
    st.plotly_chart(fig, use_container_width=True)

    with st.expander("🧠 What Does This Plot Show?", expanded=False):
        st.markdown("""
        This bar chart visualizes **chemical shift perturbation (CSP)** values across selected residues and atoms, comparing two NMR conditions.

        - **X-axis:** Residue number, filtered based on your selection.
        - **Y-axis:** CSP magnitude, computed as  
          *CSP = √((ΔH)² + (α × ΔX)²)*  
          where α is a nucleus-specific scaling factor (¹⁵N: 0.17, ¹³C: 0.251).
        - **Bar colors:** Distinguish atom names (for ¹³C); uniform color is used for ¹⁵N.
        - **Threshold lines:**  
          - Red dashed line = mean + 1σ  
          - Red dotted line = mean + 2σ
        - **Grey shaded areas:** Mark missing or unassigned residues—either because data wasn't present in your uploaded `.xpk` files or the residue was skipped during acquisition.
        """)

    # --- HIGH-RES PLOT DOWNLOADS ---
    st.subheader("🖼️ Download High-Resolution Plots")

    img_buf_clean = create_matplotlib_plot(filtered_df, nucleus, annotate=False, res_range=selected_range)
    img_buf_annotated = create_matplotlib_plot(filtered_df, nucleus, annotate=True, res_range=selected_range)

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

    # --- RESULTS TABLE + DOWNLOAD ---
    st.subheader("📋 CSP Table")
    st.dataframe(filtered_df, use_container_width=True)

    csv = filtered_df.to_csv(index=False).encode("utf-8")
    st.download_button("📥 Download Filtered CSV", data=csv, file_name="filtered_csp.csv", mime="text/csv")
