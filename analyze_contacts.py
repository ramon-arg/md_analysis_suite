"""
analyze_contacts.py
-------------------
Parses CPPTRAJ nativecontacts residue_contacts.dat output.

Residue numbering:
  Complex residues 1-230   = ROMOSOZUMAB residues 1-230
  Complex residues 231-319 = SOST residues 1-89

Produces:
  - interface_residues.txt        : unique SOST and ROMO residues in contact
  - contact_barplot.png           : bar plot sorted by TotalFrac (interface residues only)
  - contact_barplot_regions.png   : bar plot over full sequence with highlighted regions
  - contact_heatmap.png           : heatmap of the contact map
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import argparse
import os
import sys

# ─── COMPLEX NUMBERING ────────────────────────────────────────────────────────
# In the complex:  Res1 = ROMOSOZUMAB (1-230), Res2 = SOST (231-319)
ROMO_COMPLEX_RANGE = (1, 230)
SOST_COMPLEX_RANGE = (231, 319)
SOST_OFFSET = 230   # complex_res - 230 = SOST local residue number

# Full sequence lengths (local numbering)
ROMO_TOTAL = 230
SOST_TOTAL = 89
# ─────────────────────────────────────────────────────────────────────────────

STYLE = "seaborn-v0_8-whitegrid"

# Colour palette matching reference image
ROMO_COLOR  = "#E8735A"   # salmon/coral  — Romosozumab bars
SOST_COLOR  = "#5B8DB8"   # steel blue    — SOST bars
SPINE_COLOR = "#cccccc"
GRID_COLOR  = "#e5e5e5"


def parse_args():
    parser = argparse.ArgumentParser(description="Analyze CPPTRAJ nativecontacts output.")
    parser.add_argument("input", nargs="?", default="residue_contacts.dat",
                        help="Path to residue_contacts.dat (default: residue_contacts.dat)")
    parser.add_argument("--threshold", type=float, default=0.0,
                        help="Minimum TotalFrac to include a contact (default: 0.0 = all)")
    parser.add_argument("--label", type=str, default="",
                        help="Optional label appended to output filenames (e.g. cluster1)")
    return parser.parse_args()


def load_data(filepath, threshold):
    df = pd.read_csv(filepath, sep=r"\s+", comment="#",
                     names=["Res1", "Res2", "TotalFrac", "Contacts"])
    print(f"  Loaded {len(df)} contact pairs from '{filepath}'")
    if threshold > 0.0:
        before = len(df)
        df = df[df["TotalFrac"] >= threshold]
        print(f"  After threshold (>={threshold}%): {len(df)} pairs retained (removed {before - len(df)})")

    # Add local residue numbers
    df["ROMO_local"] = df["Res1"]                      # already 1-230
    df["SOST_local"] = df["Res2"] - SOST_OFFSET        # 231->1 ... 319->89
    return df


def write_residue_txt(df, outfile="interface_residues.txt"):
    romo_res = sorted(df["ROMO_local"].unique().tolist())
    sost_res = sorted(df["SOST_local"].unique().tolist())

    with open(outfile, "w") as f:
        f.write("# Interface residues identified by CPPTRAJ nativecontacts\n")
        f.write("# ROMOSOZUMAB = chain A (complex residues 1-230, local 1-230)\n")
        f.write("# SOST        = chain B (complex residues 231-319, local 1-89)\n\n")
        sost_complex = [r + SOST_OFFSET for r in sost_res]
        f.write("ROMOSOZUMAB: "    + ",".join(str(r) for r in romo_res) + "\n")
        f.write("SOST (local): "   + ",".join(str(r) for r in sost_res) + "\n")
        f.write("SOST (complex): " + ",".join(str(r) for r in sost_complex) + "\n")

    print(f"\n  ROMOSOZUMAB interface residues ({len(romo_res)}): {romo_res}")
    print(f"  SOST interface residues ({len(sost_res)}): {sost_res}")
    print(f"  Written to '{outfile}'")
    return romo_res, sost_res


# ─── STYLE HELPERS ───────────────────────────────────────────────────────────
def apply_style(ax, title, xlabel, ylabel):
    ax.set_title(title, fontsize=13, fontweight="bold", pad=10, color="#222222")
    ax.set_xlabel(xlabel, fontsize=11, color="#333333")
    ax.set_ylabel(ylabel, fontsize=11, color="#333333")
    ax.tick_params(axis="both", labelsize=9, colors="#444444")
    ax.tick_params(axis="x", rotation=45)
    ax.yaxis.grid(True, color=GRID_COLOR, linewidth=0.8, linestyle="--")
    ax.set_axisbelow(True)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    for spine in ["left", "bottom"]:
        ax.spines[spine].set_color(SPINE_COLOR)
    ax.axhline(0, color="#aaaaaa", linewidth=0.8, linestyle="--")


# ─── PLOT 1: sorted bar plot (interface residues only) ───────────────────────
def plot_barplot(df, label="", outfile="contact_barplot.png"):
    plt.style.use(STYLE)
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    fig.patch.set_facecolor("white")

    # ── ROMOSOZUMAB ──
    romo_df = (df.groupby("ROMO_local")["TotalFrac"]
                 .sum()
                 .sort_values(ascending=False)
                 .reset_index())
    romo_df.columns = ["Residue", "TotalFrac"]

    bars = axes[0].bar(romo_df["Residue"].astype(str), romo_df["TotalFrac"],
                       color=ROMO_COLOR, edgecolor="white", linewidth=0.6, width=0.7)
    for bar, val in zip(bars, romo_df["TotalFrac"]):
        axes[0].text(bar.get_x() + bar.get_width() / 2,
                     bar.get_height() + 0.4,
                     f"{val:.1f}", ha="center", va="bottom",
                     fontsize=7, color="#333333")
    apply_style(axes[0],
                "Romosozumab — Interface Contact Frequency",
                "Residue Number", "Cumulative TotalFrac (% frames)")

    # ── SOST ──
    sost_df = (df.groupby("SOST_local")["TotalFrac"]
                 .sum()
                 .sort_values(ascending=False)
                 .reset_index())
    sost_df.columns = ["Residue", "TotalFrac"]

    bars2 = axes[1].bar(sost_df["Residue"].astype(str), sost_df["TotalFrac"],
                        color=SOST_COLOR, edgecolor="white", linewidth=0.6, width=0.7)
    for bar, val in zip(bars2, sost_df["TotalFrac"]):
        axes[1].text(bar.get_x() + bar.get_width() / 2,
                     bar.get_height() + 0.4,
                     f"{val:.1f}", ha="center", va="bottom",
                     fontsize=7, color="#333333")
    apply_style(axes[1],
                "SOST — Interface Contact Frequency",
                "Residue Number", "Cumulative TotalFrac (% frames)")

    title = f"Interface Residue Contact Frequency"
    if label:
        title += f" — {label}"
    fig.suptitle(title, fontsize=15, fontweight="bold", color="#111111", y=1.02)
    plt.tight_layout()
    plt.savefig(outfile, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"  Bar plot saved to '{outfile}'")


# ─── PLOT 2: full-sequence bar plot with highlighted regions ─────────────────
def plot_barplot_regions(df, label="", outfile="contact_barplot_regions.png"):
    plt.style.use(STYLE)
    fig, axes = plt.subplots(2, 1, figsize=(18, 10))
    fig.patch.set_facecolor("white")

    # Aggregate by local residue
    romo_agg = df.groupby("ROMO_local")["TotalFrac"].sum()
    sost_agg = df.groupby("SOST_local")["TotalFrac"].sum()

    # Full residue ranges
    romo_all = pd.Series(0.0, index=range(1, ROMO_TOTAL + 1))
    sost_all = pd.Series(0.0, index=range(1, SOST_TOTAL + 1))
    romo_all.update(romo_agg)
    sost_all.update(sost_agg)

    # ── ROMOSOZUMAB (ligand) ──
    ax = axes[0]
    ax.bar(romo_all.index, romo_all.values,
           color=ROMO_COLOR, edgecolor="white", linewidth=0.3, width=1.0)

    ax.axvspan(26,  33,  alpha=0.15, color='#CC6600', label='CDR-H')
    ax.axvspan(51,  58,  alpha=0.15, color='#CC6600')
    ax.axvspan(97,  112, alpha=0.15, color='#CC6600')
    ax.axvspan(150, 156, alpha=0.15, color='#4B0082', label='CDR-L')
    ax.axvspan(173, 175, alpha=0.15, color='#4B0082')
    ax.axvspan(212, 220, alpha=0.15, color='#4B0082')

    patch_H = mpatches.Patch(color='#CC6600', alpha=0.5, label='CDR-H')
    patch_L = mpatches.Patch(color='#4B0082', alpha=0.5, label='CDR-L')
    ax.legend(handles=[patch_H, patch_L], fontsize=9,
              framealpha=0.8, loc="upper right")
    apply_style(ax,
                "Romosozumab — Contact Frequency across Full Sequence",
                "Residue Number", "Cumulative TotalFrac (% frames)")
    ax.set_xlim(0.5, ROMO_TOTAL + 0.5)
    ax.tick_params(axis="x", rotation=0)

    # ── SOST (receptor) ──
    ax2 = axes[1]
    ax2.bar(sost_all.index, sost_all.values,
            color=SOST_COLOR, edgecolor="white", linewidth=0.3, width=1.0)

    ax2.axvspan(2,  25, alpha=0.15, color='#333384', label='Loop 1')
    ax2.axvspan(31, 54, alpha=0.15, color='#aa2424', label='Loop 2')
    ax2.axvspan(56, 85, alpha=0.15, color='#267326', label='Loop 3')

    patch_L1 = mpatches.Patch(color='#333384', alpha=0.5, label='Loop 1')
    patch_L2 = mpatches.Patch(color='#aa2424', alpha=0.5, label='Loop 2')
    patch_L3 = mpatches.Patch(color='#267326', alpha=0.5, label='Loop 3')
    ax2.legend(handles=[patch_L1, patch_L2, patch_L3], fontsize=9,
               framealpha=0.8, loc="upper right")
    apply_style(ax2,
                "SOST — Contact Frequency across Full Sequence",
                "Residue Number", "Cumulative TotalFrac (% frames)")
    ax2.set_xlim(0.5, SOST_TOTAL + 0.5)
    ax2.tick_params(axis="x", rotation=0)

    title = "Contact Frequency — Full Sequence with Functional Regions"
    if label:
        title += f" — {label}"
    fig.suptitle(title, fontsize=15, fontweight="bold", color="#111111", y=1.01)
    plt.tight_layout()
    plt.savefig(outfile, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"  Regions bar plot saved to '{outfile}'")


# ─── PLOT 3: heatmap ─────────────────────────────────────────────────────────
def plot_heatmap(df, label="", outfile="contact_heatmap.png"):
    plt.style.use(STYLE)

    romo_res = sorted(df["ROMO_local"].unique())
    sost_res = sorted(df["SOST_local"].unique())

    matrix = pd.DataFrame(0.0, index=romo_res, columns=sost_res)
    for _, row in df.iterrows():
        matrix.loc[row["ROMO_local"], row["SOST_local"]] = row["TotalFrac"]

    fig, ax = plt.subplots(figsize=(max(7, len(sost_res) * 0.75),
                                    max(5, len(romo_res) * 0.5)))
    fig.patch.set_facecolor("white")

    cmap = "YlOrRd"
    im = ax.imshow(matrix.values, aspect="auto", cmap=cmap,
                   interpolation="nearest", vmin=0)

    ax.set_xticks(range(len(sost_res)))
    ax.set_xticklabels([str(r) for r in sost_res], rotation=45,
                       ha="right", fontsize=9)
    ax.set_yticks(range(len(romo_res)))
    ax.set_yticklabels([str(r) for r in romo_res], fontsize=9)
    ax.set_xlabel("SOST Residues (local 1–89)", fontsize=11, labelpad=8)
    ax.set_ylabel("Romosozumab Residues (local 1–230)", fontsize=11, labelpad=8)

    title = "Native Contact Map — TotalFrac per Residue Pair"
    if label:
        title += f" — {label}"
    ax.set_title(title, fontsize=13, fontweight="bold", pad=12)

    vmax = matrix.values.max()
    for i, r1 in enumerate(romo_res):
        for j, r2 in enumerate(sost_res):
            val = matrix.loc[r1, r2]
            if val > 0:
                ax.text(j, i, f"{val:.1f}", ha="center", va="center",
                        fontsize=7,
                        color="white" if val > vmax * 0.65 else "#333333")

    cbar = plt.colorbar(im, ax=ax, fraction=0.03, pad=0.04)
    cbar.set_label("TotalFrac (% frames)", fontsize=10)
    plt.tight_layout()
    plt.savefig(outfile, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"  Heatmap saved to '{outfile}'")


# ─── MAIN ────────────────────────────────────────────────────────────────────
def main():
    args = parse_args()
    lbl = args.label

    def outname(base):
        name, ext = os.path.splitext(base)
        return f"{name}_{lbl}{ext}" if lbl else base

    print(f"\n{'='*55}")
    print("  CPPTRAJ Native Contacts Analyzer")
    print(f"{'='*55}")

    if not os.path.exists(args.input):
        print(f"ERROR: File '{args.input}' not found.")
        sys.exit(1)

    df = load_data(args.input, args.threshold)

    if df.empty:
        print("No contacts pass the threshold. Try lowering --threshold.")
        sys.exit(0)

    write_residue_txt(df, outname("interface_residues.txt"))
    plot_barplot(df, label=lbl, outfile=outname("contact_barplot.png"))
    plot_barplot_regions(df, label=lbl, outfile=outname("contact_barplot_regions.png"))
    plot_heatmap(df, label=lbl, outfile=outname("contact_heatmap.png"))

    print(f"\n{'='*55}")
    print("  Done! Output files:")
    for f in ["interface_residues.txt", "contact_barplot.png",
              "contact_barplot_regions.png", "contact_heatmap.png"]:
        print(f"    {outname(f)}")
    print(f"{'='*55}\n")


if __name__ == "__main__":
    main()
