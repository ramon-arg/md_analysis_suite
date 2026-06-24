"""
lie_comparison.py
-----------------
Compares mean EELEC and EVDW across clusters.
Produces two plots:
  - lie_comparison_stacked.png  : stacked bar chart
  - lie_comparison_grouped.png  : grouped bar chart (side by side)

Edit INPUT_FILES below to match your file paths and labels.

Usage:
    python3 lie_comparison.py
    python3 lie_comparison.py --out myprefix
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import argparse
import os
import sys

ELEC_COLOR    = "#D4AF37"   # gold   — electrostatic
VDW_COLOR     = "#4B0082"   # purple — van der Waals

# ─── INPUT FILES — edit paths and labels here ────────────────────────────────
INPUT_FILES = [
    ("05_lie_cluster1_prod1_interface.dat", "Cluster 1"),
    ("05_lie_cluster2_prod0_interface.dat", "Cluster 2"),
    ("05_lie_cluster3_prod1_interface.dat", "Cluster 3"),
    ("05_lie_cluster4_prod1_interface.dat", "Cluster 4"),
    ("05_lie_cluster5_prod1_interface.dat", "Cluster 5"),
]
# ─────────────────────────────────────────────────────────────────────────────


def parse_args():
    parser = argparse.ArgumentParser(description="Compare LIE means across clusters.")
    parser.add_argument("--out", type=str, default="lie_comparison",
                        help="Output filename prefix (default: lie_comparison)")
    return parser.parse_args()


def load_file(filepath):
    df = pd.read_csv(filepath, sep=r"\s+", comment="#",
                     names=["Frame", "EELEC", "EVDW"])
    return df


def apply_common_style(ax, title):
    ax.set_axisbelow(True)
    ax.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.6)
    ax.yaxis.set_major_locator(ticker.MultipleLocator(100))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(50))
    ax.axhline(0, color="#aaaaaa", linewidth=0.8)
    ax.set_ylabel("Energia Média (kcal/mol)", fontsize=9)
    ax.set_xlabel("Cluster", fontsize=9)
    ax.set_title(title, fontsize=10, fontweight="bold")
    ax.tick_params(labelsize=8)
    ax.legend(fontsize=8, bbox_to_anchor=(1.01, 1),
              loc="upper left", borderaxespad=0)


def plot_stacked(summary, outfile):
    fig, ax = plt.subplots(figsize=(max(6, len(summary) * 1.4), 5))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    x = np.arange(len(summary))

    bars_elec = ax.bar(x, summary["EELEC"], color=ELEC_COLOR,
                       edgecolor="white", linewidth=0.5,
                       label="EELEC", width=0.5)
    bars_vdw  = ax.bar(x, summary["EVDW"], bottom=summary["EELEC"],
                       color=VDW_COLOR, edgecolor="white", linewidth=0.5,
                       label="EVDW", width=0.5)

    for bar, val in zip(bars_elec, summary["EELEC"]):
        ax.text(bar.get_x() + bar.get_width() / 2,
                bar.get_y() + bar.get_height() / 2,
                f"{val:.1f}", ha="center", va="center",
                fontsize=7.5, color="black", fontweight="bold")

    for bar, val in zip(bars_vdw, summary["EVDW"]):
        ax.text(bar.get_x() + bar.get_width() / 2,
                bar.get_y() + bar.get_height() / 2,
                f"{val:.1f}", ha="center", va="center",
                fontsize=7.5, color="white", fontweight="bold")

    for i, row in summary.iterrows():
        ax.text(i, row["Total"] - 5, f"Total: {row['Total']:.1f}",
                ha="center", va="top", fontsize=7, color="#333333")

    ax.set_xticks(x)
    ax.set_xticklabels(summary["label"].tolist(), fontsize=9)
    apply_common_style(ax, "Comparação LIE (Empilhado) — Interface SOST-Romosozumab")

    plt.tight_layout()
    plt.savefig(outfile, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"  Stacked plot saved to '{outfile}'")


def plot_grouped(summary, outfile):
    fig, ax = plt.subplots(figsize=(max(6, len(summary) * 1.8), 5))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    x     = np.arange(len(summary))
    width = 0.35

    bars_elec = ax.bar(x - width / 2, summary["EELEC"], width,
                       color=ELEC_COLOR, edgecolor="white",
                       linewidth=0.5, label="EELEC")
    bars_vdw  = ax.bar(x + width / 2, summary["EVDW"],  width,
                       color=VDW_COLOR,  edgecolor="white",
                       linewidth=0.5, label="EVDW")

    for bar, val in zip(bars_elec, summary["EELEC"]):
        ax.text(bar.get_x() + bar.get_width() / 2,
                bar.get_height() - 8,
                f"{val:.1f}", ha="center", va="top",
                fontsize=7.5, color="black", fontweight="bold")

    for bar, val in zip(bars_vdw, summary["EVDW"]):
        ax.text(bar.get_x() + bar.get_width() / 2,
                bar.get_height() - 3,
                f"{val:.1f}", ha="center", va="top",
                fontsize=7.5, color="white", fontweight="bold")

    ax.set_xticks(x)
    ax.set_xticklabels(summary["label"].tolist(), fontsize=9)
    apply_common_style(ax, "Comparação LIE (Agrupado) — Interface SOST-Romosozumab")

    plt.tight_layout()
    plt.savefig(outfile, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"  Grouped plot saved to '{outfile}'")


def main():
    args = parse_args()

    print(f"\n{'='*55}")
    print("  LIE Cluster Comparison")
    print(f"{'='*55}")

    records = []
    for filepath, label in INPUT_FILES:
        if not os.path.exists(filepath):
            print(f"  WARNING: '{filepath}' not found, skipping.")
            continue

        df = load_file(filepath)
        mean_elec = df["EELEC"].mean()
        mean_vdw  = df["EVDW"].mean()
        records.append({
            "label": label,
            "EELEC": mean_elec,
            "EVDW":  mean_vdw,
            "Total": mean_elec + mean_vdw,
        })
        print(f"  {label}: EELEC={mean_elec:.2f}  EVDW={mean_vdw:.2f}  "
              f"Total={mean_elec + mean_vdw:.2f}  kcal/mol")

    if not records:
        print("ERROR: No valid files loaded.")
        sys.exit(1)

    summary = pd.DataFrame(records)
    print(f"\n  {'─'*45}")
    print(f"  {'Cluster':<20} {'EELEC':>10} {'EVDW':>10} {'Total':>10}")
    print(f"  {'─'*45}")
    for _, row in summary.iterrows():
        print(f"  {row['label']:<20} {row['EELEC']:>10.2f} "
              f"{row['EVDW']:>10.2f} {row['Total']:>10.2f}")
    print(f"  {'─'*45}\n")

    plot_stacked(summary, f"{args.out}_stacked.png")
    plot_grouped(summary, f"{args.out}_grouped.png")

    print(f"{'='*55}\n")


if __name__ == "__main__":
    main()
