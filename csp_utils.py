import math
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import plotly.express as px
from io import BytesIO
from itertools import groupby
from operator import itemgetter

# --- Unified Parser --- NOT IN USE ANYMORE
def parse_peaklist(file_obj):
    file_obj.seek(0)  # 🔄 Reset pointer before reading
    lines = file_obj.read().decode("utf-8").splitlines()[6:]
    data = []
    for line in lines:
        parts = line.strip().split()
        if len(parts) < 17 or parts[1] == '{}':
            continue
        try:
            label = parts[1]
            residue = int(label.split('{')[1].split('.')[0])
            atom = label.split('{')[1].split('.')[1].rstrip('}')
            h = float(parts[2])
            x = float(parts[9])
            intensity = float(parts[16])
            data.append({
                "Residue": residue,
                "Atom": atom,
                "H": h,
                "X": x,
                "Intensity": intensity
            })
        except:
            continue
    return pd.DataFrame(data)

def parse_peaklist(file_obj, format="xpk"):
    import pandas as pd
    file_obj.seek(0)
    lines = file_obj.read().decode("utf-8").splitlines()
    data = []

    if format == "xpk":
        lines = lines[6:]
        for line in lines:
            parts = line.strip().split()
            if len(parts) < 17 or parts[1] == '{}':
                continue
            try:
                label = parts[1]
                residue = int(label.split('{')[1].split('.')[0])
                atom = label.split('{')[1].split('.')[1].rstrip('}')
                h = float(parts[2])
                x = float(parts[9])
                intensity = float(parts[16])
                data.append({
                    "Residue": residue,
                    "Atom": atom,
                    "H": h,
                    "X": x,
                    "Intensity": intensity
                })
            except:
                continue

    elif format == "nmrfx":
        for line in lines:
            parts = line.strip().split()
            if len(parts) < 5:
                continue
            try:
                label = parts[1]
                residue = int(label.split('.')[0])
                atom = label.split('.')[1]
                h = float(parts[2])
                x = float(parts[6])
                intensity = float(parts[10])
                data.append({
                    "Residue": residue,
                    "Atom": atom,
                    "H": h,
                    "X": x,
                    "Intensity": intensity
                })
            except:
                continue

    elif format == "sparky":
        for line in lines:
            parts = line.strip().split()
            if len(parts) < 4:
                continue
            try:
                h = float(parts[0])
                x = float(parts[1])
                label = parts[2]
                residue = int(''.join(filter(str.isdigit, label)))
                atom = ''.join(filter(str.isalpha, label))
                intensity = float(parts[3]) if len(parts) >= 4 else 0.0
                data.append({
                    "Residue": residue,
                    "Atom": atom,
                    "H": h,
                    "X": x,
                    "Intensity": intensity
                })
            except:
                continue

    return pd.DataFrame(data)


# --- CSP Analysis ---
def run_csp_analysis(initial_file, final_file, nucleus, x_scale):
    df_init = parse_peaklist(initial_file)
    df_final = parse_peaklist(final_file)

    merged = pd.merge(
        df_init, df_final,
        on=["Residue", "Atom"],
        suffixes=("_init", "_final")
    )
    merged = merged[merged["Residue"] >= 1]

    merged["CSP"] = merged.apply(
        lambda row: math.sqrt((row["H_final"] - row["H_init"])**2 + (x_scale * (row["X_final"] - row["X_init"]))**2),
        axis=1
    )
    merged["ΔH"] = merged["H_final"] - merged["H_init"]
    merged["ΔX"] = merged["X_final"] - merged["X_init"]

    return merged[["Residue", "Atom", "CSP", "ΔH", "ΔX"]].sort_values(["Residue", "Atom"])

# --- Intensity Analysis ---
def run_intensity_analysis(initial_file, final_file, normalize_residue, final_scale):
    df_init = parse_peaklist(initial_file)
    df_final = parse_peaklist(final_file)

    df_init = df_init[df_init["Residue"] >= 1]
    df_final = df_final[df_final["Residue"] >= 1]

    df_init_norm = df_init.groupby("Residue")["Intensity"].mean().reset_index()
    df_final_norm = df_final.groupby("Residue")["Intensity"].mean().reset_index()

    merged = pd.merge(df_init_norm, df_final_norm, on="Residue", suffixes=("_init", "_final"))

    init_ref = merged.loc[merged["Residue"] == normalize_residue, "Intensity_init"]
    final_ref = merged.loc[merged["Residue"] == normalize_residue, "Intensity_final"]

    init_ref_val = init_ref.values[0] if not init_ref.empty else 1.0
    final_ref_val = final_ref.values[0] if not final_ref.empty else 1.0

    merged["Initial_Norm"] = merged["Intensity_init"] / init_ref_val if init_ref_val != 0 else 0
    merged["Final_Norm"] = (merged["Intensity_final"] / final_ref_val) * final_scale if final_ref_val != 0 else 0
    merged["Ratio"] = merged.apply(
        lambda row: row["Final_Norm"] / row["Initial_Norm"] if row["Initial_Norm"] != 0 else 0,
        axis=1
    )

    return merged[["Residue", "Initial_Norm", "Final_Norm", "Ratio"]]


# --- CSP Plot ---
def create_plotly_chart(df, nucleus, thresholds=None, res_range=None):
    min_res = res_range[0] if res_range else int(df["Residue"].min())
    max_res = res_range[1] if res_range else int(df["Residue"].max())
    full_range = range(min_res, max_res + 1)
    present_residues = sorted(df["Residue"].unique())
    missing_residues = sorted(set(full_range) - set(present_residues))

    color_mode = "Atom" if nucleus == "C" else None
    fig = px.bar(
        df,
        x="Residue",
        y="CSP",
        color=color_mode,
        hover_data=["Atom", "ΔH", "ΔX"],
        color_discrete_sequence=px.colors.qualitative.Set2,
        labels={"CSP": "Chemical Shift Perturbation", "Residue": "Residue Number"},
        title="Interactive CSP Profile"
    )
    fig.update_layout(
        template="simple_white",
        xaxis=dict(tickmode="linear", dtick=5, range=[min_res - 1, max_res + 1]),
        yaxis=dict(title="CSP", rangemode="tozero")
    )

    for k, g in groupby(enumerate(missing_residues), lambda ix: ix[0] - ix[1]):
        gap = list(map(itemgetter(1), g))
        fig.add_vrect(x0=gap[0] - 0.5, x1=gap[-1] + 0.5, fillcolor="lightgrey", opacity=0.25, layer="below", line_width=0)

    if thresholds:
        fig.add_shape(type="line", x0=min_res - 1, x1=max_res + 1, y0=thresholds[0], y1=thresholds[0],
                      line=dict(color="red", width=2, dash="dash"))
        fig.add_shape(type="line", x0=min_res - 1, x1=max_res + 1, y0=thresholds[1], y1=thresholds[1],
                      line=dict(color="firebrick", width=2, dash="dot"))

    return fig

# --- Intensity Plot ---
def create_intensity_plot(df, res_range=None):
    import plotly.express as px
    import plotly.graph_objects as go

    df_plot = df.copy()

    min_res = res_range[0] if res_range else int(df_plot["Residue"].min())
    max_res = res_range[1] if res_range else int(df_plot["Residue"].max())
    full_range = set(range(min_res, max_res + 1))
    present = set(df_plot["Residue"])
    missing = sorted(full_range - present)

    fig = px.bar(
        df_plot,
        x="Residue",
        y="Ratio",
        hover_data=["Initial_Norm", "Final_Norm"],
        color="Ratio",
        color_continuous_scale="Blues",
        title="Normalized Intensity Ratio (Final / Initial)"
    )

    y_max = df_plot["Ratio"].max() * 1.1  # Or round it up slightly for headroom

    fig.update_layout(
        template="simple_white",
        xaxis=dict(tickmode="linear", dtick=5),
        yaxis=dict(range=[0, y_max], zeroline=True, zerolinewidth=1),
        xaxis_title="Residue",
        yaxis_title="Intensity Ratio"
    )



    # Add shaded regions for unassigned residues
    from itertools import groupby
    from operator import itemgetter

    for _, g in groupby(enumerate(missing), lambda ix: ix[0] - ix[1]):
        gap = list(map(itemgetter(1), g))
        fig.add_vrect(
            x0=gap[0] - 0.5,
            x1=gap[-1] + 0.5,
            fillcolor="lightgrey",
            opacity=0.25,
            layer="below",
            line_width=0
        )

    return fig


def create_matplotlib_plot(df, mode="CSP", annotate=False, res_range=None, thresholds=None):
    import matplotlib.pyplot as plt
    import numpy as np
    import math
    from io import BytesIO
    from itertools import groupby
    from operator import itemgetter

    # Filter valid residues
    df = df[df["Residue"] >= 1].copy()
    df.sort_values(["Residue", "Atom"] if "Atom" in df.columns else "Residue", inplace=True)

    # Extract x/y data
    x = df["Residue"]
    y = df["CSP"] if mode == "CSP" else df["Ratio"]
    ylabel = "Chemical Shift Perturbation (CSP)" if mode == "CSP" else "Normalized Intensity Ratio"
    title = "CSP vs Residue Number" if mode == "CSP" else "Intensity Ratio (Final / Initial)"

    # Determine range
    min_res = res_range[0] if res_range else int(x.min()) if not x.empty else 1
    max_res = res_range[1] if res_range else int(x.max()) if not x.empty else min_res
    full_range = range(min_res, max_res + 1)
    df = df[(df["Residue"] >= min_res) & (df["Residue"] <= max_res)]
    x = df["Residue"]
    y = df["CSP"] if mode == "CSP" else df["Ratio"]

    # Setup plot
    fig_width = min(max(8, len(full_range) * 0.4), 20)
    fig, ax = plt.subplots(figsize=(fig_width, 6))
    bar_color = "steelblue" if mode == "CSP" else "royalblue"
    ax.bar(x, y, color=bar_color, edgecolor="black", linewidth=0.8)

    # Shade unassigned residues
    missing = sorted(set(full_range) - set(x))
    for _, g in groupby(enumerate(missing), lambda ix: ix[0] - ix[1]):
        gap = list(map(itemgetter(1), g))
        ax.axvspan(gap[0] - 0.5, gap[-1] + 0.5, color="lightgrey", alpha=0.4)

    # X-axis styling
    ax.set_xlim(min_res - 0.5, max_res + 0.5)
    tick_start = (min_res // 5) * 5
    ax.set_xticks(np.arange(tick_start, max_res + 1, 5))
    ax.set_xlabel("Residue")
    ax.tick_params(axis='x', rotation=45)

    # Y-axis styling
    if mode == "CSP":
        ymax = max(0.1, y.max())
        ymax = round(ymax + 0.05, 2)
        ax.set_ylim(0, ymax)
        ax.set_yticks(np.arange(0, ymax + 0.01, 0.1))
    else:
        ymax = max(2.0, y.max())
        ymax = round(ymax * 1.1, 2)
        ax.set_ylim(0, ymax)
        ax.set_yticks(np.arange(0, ymax + 0.1, 0.5))

    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(axis='x', linestyle='--', alpha=0.4)

    # Threshold lines
    if annotate and thresholds:
        t1, t2 = thresholds
        if mode == "CSP":
            ax.axhline(t1, color="red", linestyle="--", linewidth=1.5)
            ax.text(min_res, t1 + 0.01, f"1σ: {t1:.2f}", color="red", fontsize=10)
            ax.axhline(t2, color="firebrick", linestyle=":", linewidth=1.5)
            ax.text(min_res, t2 + 0.01, f"2σ: {t2:.2f}", color="firebrick", fontsize=10)
        else:
            ax.axhline(t1, color="green", linestyle="--", linewidth=1.5)
            ax.text(min_res, t1 + 0.05, f"{t1}×", color="green", fontsize=10)
            ax.axhline(t2, color="darkgreen", linestyle=":", linewidth=1.5)
            ax.text(min_res, t2 + 0.05, f"{t2}×", color="darkgreen", fontsize=10)

    # Output as high-res PNG
    plt.tight_layout()
    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=300)
    plt.close()
    buf.seek(0)
    return buf

def create_combined_matplotlib_plot(df_csp, df_intensity, res_range=None, thresholds_csp=None, thresholds_ratio=None):
    import matplotlib.pyplot as plt
    import numpy as np
    from io import BytesIO

    df_csp = df_csp[df_csp["Residue"] >= 1].copy()
    df_intensity = df_intensity[df_intensity["Residue"] >= 1].copy()

    if res_range:
        df_csp = df_csp[(df_csp["Residue"] >= res_range[0]) & (df_csp["Residue"] <= res_range[1])]
        df_intensity = df_intensity[(df_intensity["Residue"] >= res_range[0]) & (df_intensity["Residue"] <= res_range[1])]

    x_csp = df_csp["Residue"]
    y_csp = df_csp["CSP"]
    x_int = df_intensity["Residue"]
    y_int = df_intensity["Ratio"]

    fig_width = min(max(10, (res_range[1] - res_range[0]) * 0.5), 22) if res_range else 14
    fig, axes = plt.subplots(2, 1, figsize=(fig_width, 10), sharex=True)

    # --- CSP Plot ---
    axes[0].bar(x_csp, y_csp, color="steelblue", edgecolor="black", linewidth=0.8)
    axes[0].set_ylabel("CSP")
    axes[0].set_title("CSP vs Residue Number")
    axes[0].grid(axis='x', linestyle='--', alpha=0.4)
    ymax_csp = np.round(y_csp.max() + 0.1, 2)
    axes[0].set_ylim(0, ymax_csp)
    if thresholds_csp:
        t1, t2 = thresholds_csp
        axes[0].axhline(t1, color="red", linestyle="--")
        axes[0].axhline(t2, color="firebrick", linestyle=":")

    # --- Intensity Plot ---
    axes[1].bar(x_int, y_int, color="royalblue", edgecolor="black", linewidth=0.8)
    axes[1].set_ylabel("Intensity Ratio")
    axes[1].set_title("Intensity Ratio (Final / Initial)")
    axes[1].grid(axis='x', linestyle='--', alpha=0.4)
    ymax_int = np.round(max(2.0, y_int.max()) * 1.1, 2)
    axes[1].set_ylim(0, ymax_int)
    if thresholds_ratio:
        t1, t2 = thresholds_ratio
        axes[1].axhline(t1, color="green", linestyle="--")
        axes[1].axhline(t2, color="darkgreen", linestyle=":")

    # X-axis config

    min_res = res_range[0] if res_range else int(min(x_csp.min(), x_int.min()))
    max_res = res_range[1] if res_range else int(max(x_csp.max(), x_int.max()))
    tick_start = (min_res // 5) * 5

    x_min = int(min(x_csp.min(), x_int.min()))
    x_max = int(max(x_csp.max(), x_int.max()))
    axes[1].set_xlabel("Residue")
    tick_start = (min_res // 5) * 5
    axes[1].set_xticks(np.arange(tick_start, max_res + 1, 5))
    axes[1].set_xlim(min_res - 0.5, max_res + 0.5)
    axes[1].tick_params(axis='x', rotation=45)

    plt.tight_layout()
    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=300)
    plt.close()
    buf.seek(0)
    return buf

