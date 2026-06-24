"""
heatmap_contacts.py
-------------------
Generates a native contacts heatmap from CPPTRAJ nativecontacts output.

Tweaks vs original:
  - TotalFrac < 1.0 filtered out (not plotted)
  - Square cells (not rectangles)
  - Larger font for in-cell TotalFrac values
  - Column (SOST) and row (ROMO) gridlines colored by functional region
  - Portuguese axis labels and title
  - X-axis tick labels upright (no rotation)
  - Larger tick label fonts

Residue numbering:
  Complex residues 1-230   = ROMOSOZUMAB residues 1-230
  Complex residues 231-319 = SOST residues 1-89
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
import numpy as np
import argparse
import os
import sys

# ─── COMPLEX NUMBERING ────────────────────────────────────────────────────────
SOST_OFFSET = 230

# Functional region definitions
SOST_REGIONS = [
    (2,  25, '#333384', 'Loop 1'),
    (31, 54, '#aa2424', 'Loop 2'),
    (56, 85, '#267326', 'Loop 3'),
]

ROMO_REGIONS = [
    (26,  33,  '#CC6600', 'CDR-H1'),
    (51,  58,  '#CC6600', 'CDR-H2'),
    (97,  112, '#CC6600', 'CDR-H3'),
    (150, 156, '#4B0082', 'CDR-L1'),
    (173, 175, '#4B0082', 'CDR-L2'),
    (212, 220, '#4B0082', 'CDR-L3'),
]

ROMO_REGION_LABELS = {
    '#CC6600': 'CDR-H',
    '#4B0082': 'CDR-L',
}


def residue_region_color(res, regions, default="#dddddd"):
    """Return the color of the functional region a residue belongs to, or default."""
    for start, end, color, _ in regions:
        if start <= res <= end:
            return color
    return default


def parse_args():
    parser = argparse.ArgumentParser(description="Native contacts heatmap.")
    parser.add_argument("input", nargs="?", default="residue_contacts.dat",
                        help="Path to residue_contacts.dat")
    parser.add_argument("--label", type=str, default="",
                        help="Optional label appended to output filename")
    return parser.parse_args()


def load_data(filepath, threshold=5.0):
    df = pd.read_csv(filepath, sep=r"\s+", comment="#",
                     names=["Res1", "Res2", "TotalFrac", "Contacts"])
    before = len(df)
    df = df[df["TotalFrac"] >= threshold]
    print(f"  Loaded {before} pairs; {len(df)} retained after TotalFrac >= {threshold} filter")
    df["ROMO_local"] = df["Res1"]
    df["SOST_local"] = df["Res2"] - SOST_OFFSET
    return df


def plot_heatmap(df, label="", outfile="contact_heatmap.png"):
    romo_res = sorted(df["ROMO_local"].unique())
    sost_res = sorted(df["SOST_local"].unique())

    n_romo = len(romo_res)
    n_sost = len(sost_res)

    matrix = pd.DataFrame(0.0, index=romo_res, columns=sost_res)
    for _, row in df.iterrows():
        matrix.loc[row["ROMO_local"], row["SOST_local"]] = row["TotalFrac"]

    # ── Square cells: fix figure size so each cell is square ──
    cell_size = 0.75          # inches per cell
    fig_w = n_sost * cell_size + 3.5   # extra for colorbar + y-axis labels
    fig_h = n_romo * cell_size + 2.5   # extra for x-axis labels + title

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    fig.patch.set_facecolor("white")

    cmap = "YlOrRd"
    vmax = matrix.values.max()
    im = ax.imshow(matrix.values, aspect="equal", cmap=cmap,
                   interpolation="nearest", vmin=0, vmax=vmax)

    # ── Axis ticks ──
    ax.set_xticks(range(n_sost))
    ax.set_xticklabels([str(r) for r in sost_res],
                       rotation=0, ha="center", fontsize=11)
    ax.set_yticks(range(n_romo))
    ax.set_yticklabels([str(r) for r in romo_res], fontsize=11)

    ax.set_xlabel("Resíduos da esclerostina", fontsize=13, labelpad=10)
    ax.set_ylabel("Resíduos de Romosozumab",  fontsize=13, labelpad=10)
    ax.set_title(
        "Mapa de Contatos Nativos - TotalFrac por Par de Resíduo para Complexo II",
        fontsize=13, fontweight="bold", pad=14
    )

    # ── Colored bands centered on each cell by functional region ──
    band_lw = 1.2  # line thickness in points

    # Vertical bands — SOST regions (columns, x-axis)
    for j, res in enumerate(sost_res):
        color = residue_region_color(res, SOST_REGIONS)
        ax.axvline(x=j, color=color, linewidth=band_lw, alpha=1.0, zorder=0)

    # Horizontal bands — ROMO regions (rows, y-axis)
    for i, res in enumerate(romo_res):
        color = residue_region_color(res, ROMO_REGIONS)
        ax.axhline(y=i, color=color, linewidth=band_lw, alpha=1.0, zorder=0)

    # ── In-cell text (zorder=5 keeps numbers above lines) ──
    for i, r1 in enumerate(romo_res):
        for j, r2 in enumerate(sost_res):
            val = matrix.loc[r1, r2]
            if val > 0:
                text_color = "white" if val > vmax * 0.65 else "#333333"
                ax.text(j, i, f"{val:.1f}",
                        ha="center", va="center",
                        fontsize=11, fontweight="bold",
                        color=text_color, zorder=5)

    # ── Colorbar ──
    cbar = plt.colorbar(im, ax=ax, fraction=0.03, pad=0.04)
    cbar.set_label("TotalFrac (% frames)", fontsize=11)
    cbar.ax.tick_params(labelsize=10)

    # ── Single merged legend — bottom right inside the axes ──
    sost_patches = [
        mpatches.Patch(color=c, alpha=0.9, label=f"SOST {lbl}")
        for _, _, c, lbl in SOST_REGIONS
    ]
    romo_seen = {}
    romo_patches = []
    for _, _, c, lbl in ROMO_REGIONS:
        display = ROMO_REGION_LABELS.get(c, lbl)
        if display not in romo_seen:
            romo_patches.append(mpatches.Patch(color=c, alpha=0.9, label=f"ROMO {display}"))
            romo_seen[display] = True

    from matplotlib.lines import Line2D

    romo_line_seen = {}
    romo_lines = []
    for _, _, c, lbl in ROMO_REGIONS:
        display = ROMO_REGION_LABELS.get(c, lbl)
        if display not in romo_line_seen:
            romo_lines.append(Line2D([0], [0], color=c, linewidth=2, label=f"ROMO {display}"))
            romo_line_seen[display] = True

    sost_lines = [
        Line2D([0], [0], color=c, linewidth=2, label=f"SOST {lbl}")
        for _, _, c, lbl in SOST_REGIONS
    ]

    all_handles = romo_lines + sost_lines
    ax.legend(handles=all_handles,
              loc="lower left",
              bbox_to_anchor=(1.04, 0.0),
              bbox_transform=ax.transAxes,
              fontsize=9, title_fontsize=10,
              title="Regiões funcionais",
              framealpha=0.92,
              edgecolor="#aaaaaa")

    plt.tight_layout()
    plt.savefig(outfile, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"  Heatmap saved to '{outfile}'")


def main():
    args = parse_args()

    def outname(base):
        name, ext = os.path.splitext(base)
        return f"{name}_{args.label}{ext}" if args.label else base

    print(f"\n{'='*55}")
    print("  Native Contacts Heatmap")
    print(f"{'='*55}")

    if not os.path.exists(args.input):
        print(f"ERROR: File '{args.input}' not found.")
        sys.exit(1)

    df = load_data(args.input, threshold=5.0)

    if df.empty:
        print("No contacts with TotalFrac >= 1.0 found.")
        sys.exit(0)

    plot_heatmap(df, label=args.label, outfile=outname("new_contact_heatmap.png"))

    print(f"\n{'='*55}")
    print(f"  Done! → {outname('new_contact_heatmap.png')}")
    print(f"{'='*55}\n")


if __name__ == "__main__":
    main()
