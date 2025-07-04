import streamlit as st
import pandas as pd
from csp_utils import (
    run_csp_analysis,
    create_plotly_chart,
    create_matplotlib_plot,
    run_intensity_analysis,
    create_combined_matplotlib_plot,
    create_intensity_plot
)

st.set_page_config(layout="wide", initial_sidebar_state="expanded")
st.title("🧪 CSP & Intensity Analyser")

with st.sidebar:
    st.header("📂 Upload Files")
    format = st.selectbox("Peaklist Format", ["xpk", "nmrfx", "sparky"])
    initial_file = st.file_uploader("Initial .xpk File", type="xpk")
    final_file = st.file_uploader("Final .xpk File", type="xpk")
    nucleus = st.selectbox("Nucleus Type", ["N", "C"])

st.markdown("""
> **Welcome to the CSP & Intensity Analyzer App**  
>
> Upload two `.xpk` NMR peak files (Initial and Final), select the nucleus type, and the app will compute:
> - **Chemical Shift Perturbation (CSP)** profiles with interactive and exportable plots
> - **Intensity Ratio (Final/Initial)** analysis with normalization
> - High-resolution PNG and LaTeX outputs for figure-ready results
>
> Use the sidebar to apply filters, set thresholds, or customize your plotting range.
""")

if initial_file and final_file and nucleus:
    x_scale = 0.17 if nucleus == "N" else 0.251
    df = run_csp_analysis(initial_file, final_file, nucleus, x_scale)
    st.success("✅ CSP Analysis Complete")

    # Filter controls
    st.sidebar.header("🔎 Filter CSP")
    max_res = int(df["Residue"].max())
    selected_range = st.slider("Residue Range", 1, max_res, (1, max_res))
    atom_list = sorted(df["Atom"].unique())
    selected_atoms = st.multiselect("Atoms", atom_list, default=atom_list)

    filtered_df = df[
        (df["Residue"] >= selected_range[0]) &
        (df["Residue"] <= selected_range[1]) &
        (df["Atom"].isin(selected_atoms))
    ]

    # CSP thresholds
    mean_csp = df["CSP"].mean()
    std_csp = df["CSP"].std()
    threshold1 = mean_csp + std_csp
    threshold2 = mean_csp + 2 * std_csp

    with st.sidebar.expander("📊 CSP Stats", expanded=True):
        st.markdown(f"**Mean CSP**: `{mean_csp:.5f}`")
        st.markdown(f"**Std Dev**: `{std_csp:.5f}`")
        st.markdown(f"**Matched rows**: `{len(df)}`")
        st.markdown(f"**Zero CSPs**: `{sum(df['CSP'] == 0)}`")
        st.download_button(
            label="📥 Download All CSPs",
            data=df.to_csv(index=False).encode("utf-8"),
            file_name="full_matched_csp.csv",
            mime="text/csv"
        )

    # CSP plot + tables
    st.header("🔬 Chemical Shift Perturbation")
    fig = create_plotly_chart(filtered_df, nucleus, thresholds=(threshold1, threshold2), res_range=selected_range)
    st.plotly_chart(fig, use_container_width=True)
    st.dataframe(filtered_df, use_container_width=True)

    st.download_button("📥 Download Filtered CSPs", data=filtered_df.to_csv(index=False).encode("utf-8"),
                       file_name="filtered_csp.csv", mime="text/csv")

    # CSP Matplotlib
    st.subheader("🖼️ CSP Plots (PNG)")
    clean_csp = create_matplotlib_plot(filtered_df, mode="CSP", annotate=False,
                                       res_range=selected_range, thresholds=(threshold1, threshold2))
    annotated_csp = create_matplotlib_plot(filtered_df, mode="CSP", annotate=True,
                                           res_range=selected_range, thresholds=(threshold1, threshold2))
    col1, col2 = st.columns(2)
    with col1:
        st.download_button("Clean CSP Plot", clean_csp, "csp_clean.png", "image/png")
    with col2:
        st.download_button("Annotated CSP Plot", annotated_csp, "csp_annotated.png", "image/png")

    # Intensity Analysis
    st.header("📈 Intensity Ratio Analysis")
    with st.expander("⚙️ Intensity Settings", expanded=True):
        norm_res = st.number_input("Residue for normalization", min_value=1, value=1)
        scale_factor = st.number_input("Scaling factor (final intensities)", min_value=0.0, value=1.0, step=0.1)

    if st.button("Run Intensity Analysis"):
        intensity_df = run_intensity_analysis(initial_file, final_file, norm_res, scale_factor)

        st.subheader("📊 Intensity Ratio Plot")
        fig_intensity = create_intensity_plot(intensity_df, res_range=selected_range)
        st.plotly_chart(fig_intensity, use_container_width=True)

        st.dataframe(intensity_df, use_container_width=True)
        st.download_button("📥 Download Intensity CSV", intensity_df.to_csv(index=False).encode("utf-8"),
                           file_name="intensity_ratio.csv", mime="text/csv")
        
        combined_buf = create_combined_matplotlib_plot(
            df_csp=filtered_df,
            df_intensity=intensity_df,
            res_range=selected_range,
            thresholds_csp=(threshold1, threshold2),
            thresholds_ratio=(0.5, 2.0)
        )

        st.subheader("🖼️ Side-by-Side CSP + Intensity Plot")


        st.download_button(
            "📥 Download Combined CSP + Intensity Plot",
            combined_buf,
            file_name="combined_csp_intensity.png",
            mime="image/png"
        )

        # Intensity Matplotlib
        st.subheader("🖼️ Intensity Ratio Plots (PNG)")
        clean_ratio = create_matplotlib_plot(intensity_df, mode="Ratio", annotate=False,
                                             res_range=selected_range, thresholds=(0.5, 2.0))
        annotated_ratio = create_matplotlib_plot(intensity_df, mode="Ratio", annotate=True,
                                                 res_range=selected_range, thresholds=(0.5, 2.0))
        col3, col4 = st.columns(2)
        with col3:
            st.download_button("Clean Intensity Plot", clean_ratio, "intensity_clean.png", "image/png")
        with col4:
            st.download_button("Annotated Intensity Plot", annotated_ratio, "intensity_annotated.png", "image/png")

