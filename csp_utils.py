import math
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import plotly.express as px
from plotly.graph_objects import Figure
from io import BytesIO
from itertools import groupby
from operator import itemgetter

def parse_xpk(file_obj):
    lines = file_obj.read().decode("utf-8").splitlines()[6:]
    data = []
    for line in lines:
        parts = line.strip().split()
        if len(parts) < 17 or parts[1] == '{}':
            continue
        try:
            label = parts[1]
            residue = int(label.split('{')[1].split('.')[0])
            atom_name = label.split('{')[1].split('.')[1].rstrip('}')
            h_shift = float(parts[2])
            x_shift = float(parts[9])
            data.append((residue, atom_name, h_shift, x_shift))
        except:
            continue
    return data

def calculate_csp(h1, h2, x1, x2, scale):
    return math.sqrt((h1 - h2) ** 2 + (scale * (x1 - x2)) ** 2)

def run_csp_analysis(initial_file, final_file, nucleus, x_scale):
    initial_data = parse_xpk(initial_file)
    final_data = parse_xpk(final_file)

    results = []
    for res1, atom1, h1, x1 in initial_data:
        for res2, atom2, h2, x2 in final_data:
            if res1 == res2 and atom1 == atom2:
                csp = calculate_csp(h1, h2, x1, x2, x_scale)
                results.append({
                    "Residue": res1,
                    "Atom": atom1,
                    "CSP": csp,
                    "ΔH": h2 - h1,
                    "ΔX": x2 - x1
                })

    df = pd.DataFrame(results)
    df.sort_values(["Residue", "Atom"], inplace=True)
    return df

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
        yaxis=dict(title="CSP", rangemode="tozero"),
        margin=dict(t=60, b=40),
        legend_title="Atom" if color_mode else None,
        showlegend=True
    )

    for k, g in groupby(enumerate(missing_residues), lambda ix: ix[0] - ix[1]):
        gap = list(map(itemgetter(1), g))
        fig.add_vrect(
            x0=gap[0] - 0.5,
            x1=gap[-1] + 0.5,
            fillcolor="lightgrey",
            opacity=0.25,
            layer="below",
            line_width=0,
        )

    if thresholds:
        fig.add_shape(
            type="line",
            x0=min_res - 1,
            x1=max_res + 1,
            y0=thresholds[0],
            y1=thresholds[0],
            line=dict(color="red", width=2, dash="dash")
        )
        fig.add_shape(
            type="line",
            x0=min_res - 1,
            x1=max_res + 1,
            y0=thresholds[1],
            y1=thresholds[1],
            line=dict(color="firebrick", width=2, dash="dot")
        )

    return fig

def create_matplotlib_plot(df, nucleus, annotate=False, res_range=None, thresholds=None):
    df = df[df["Residue"] >= 1]
    df = df.sort_values(["Residue", "Atom"])
    residues = sorted(df["Residue"].unique())
    atoms = sorted(df["Atom"].unique())
    n_atoms = len(atoms)

    min_res = res_range[0] if res_range else min(residues)
    max_res = res_range[1] if res_range else max(residues)
    full_range = range(min_res, max_res + 1)
    missing_residues = sorted(set(full_range) - set(residues))

    colors = plt.cm.get_cmap('tab20').colors
    scale_factor = 0.4 if nucleus == "C" else 0.25
    fig_width = min(max(8, len(full_range) * scale_factor), 20)
    fig, ax = plt.subplots(figsize=(fig_width, 6))

    for k, g in groupby(enumerate(missing_residues), lambda ix: ix[0] - ix[1]):
        gap = list(map(itemgetter(1), g))
        ax.axvspan(gap[0] - 0.5, gap[-1] + 0.5, color='lightgrey', alpha=0.4, zorder=0)

    bar_width = max(0.7 / n_atoms, 0.3)
    for i, atom in enumerate(atoms):
        df_atom = df[df["Atom"] == atom]
        offset = i * bar_width - (bar_width * (n_atoms - 1) / 2)
        positions = df_atom["Residue"] + offset if n_atoms > 1 else df_atom["Residue"]
        ax.bar(
            positions,
            df_atom["CSP"],
            width=bar_width,
            color="steelblue" if nucleus == "N" else colors[i % len(colors)],
            edgecolor='black',
            linewidth=0.8,
            label=atom if nucleus == "C" else None
        )

    if nucleus == "C":
        ax.legend(title="Atom", bbox_to_anchor=(1.05, 1), loc='upper left')

    ymax = math.ceil(df["CSP"].max() * 10) / 10
    ax.set_ylim(0, ymax)
    ax.set_yticks(np.arange(0, ymax + 0.05, 0.1))
    ax.set_xlabel("Residue Number")
    ax.set_ylabel("Chemical Shift Perturbation (CSP)")
    ax.set_title("CSP vs Residue Number")
    ax.set_xticks(np.arange(min_res, max_res + 1, 5))
    ax.tick_params(axis='x', rotation=45)
    ax.grid(axis='x', linestyle='--', alpha=0.5)

    if annotate and thresholds:
        threshold1, threshold2 = thresholds
        ax.axhline(y=threshold1, color="red", linestyle="--", linewidth=1.5)
        ax.text(min_res, threshold1 + 0.01, f"1σ: {threshold1:.2f}", color="red", fontsize=10, va="bottom")
        ax.axhline(y=threshold2, color="firebrick", linestyle=":", linewidth=1.5)
        ax.text(min_res, threshold2 + 0.01, f"2σ: {threshold2:.2f}", color="firebrick", fontsize=10, va="bottom")

    plt.tight_layout()
    buf = BytesIO()
    fig.savefig(buf, format="png", dpi=300)
    plt.close()
    buf.seek(0)
    return buf
