# Supp Fig 8) Myeloid DA at patient level - per-donor abundance box plots of myeloid cell subsets (HBC down / mMac up in complications)
# R1 comment 3: Myeloid cell abundances sample level data
#=========================================================
# Box plots that show the abundance of every cell type in every donor split by conditions.
# I want to highlight that HBC go down during complications while mMacs go up during complications.
# Use file in "/projects/Pregnancy_placenta/Results/23_Supp_tables/DA/": AllMyeloid_comparisons_stats.csv, AllMyeloid_frequencies_per_pt.csv

import os
import pandas as pd
import scipy.stats
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

matplotlib.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']

os.chdir("/projects/Pregnancy_placenta/Results/25_Review/R1_C3_myeloid_abundnace/")

DA_dir = "/projects/Pregnancy_placenta/Results/23_Supp_tables/DA/"
# Patient_ID <-> Condition_group, many-to-many (a patient with both PrePE and FGR gets 2 rows)
code_dir="/projects/Pregnancy_placenta/Code_Publication/"
PATIENT_CONDITION_FILE = code_dir + "AuxFiles/patient_condition_list.csv"

# Controls first, then complications, so the HBC-down / mMac-up trend reads left-to-right
CONDITION_ORDER = ["Cont_labored", "Cont_unlabored", "Cont_PTB", "PTB", "DM1", "FGR", "PrePE", "TermPE"]
condition_colors = dict(zip(CONDITION_ORDER, ["#868686FF", "#868686FF", "#868686FF", "#CD534CFF", "#CD534CFF", "#CD534CFF", "#CD534CFF", "#CD534CFF"]))

# Same 13 standard case/control pairs used throughout the DA pipeline (24_Differential_abundance.py,
# Publication_code/Figure_3.py) - (case, control), case listed first to match the "Comparison"
# column in AllMyeloid_comparisons_stats.csv ("-".join([case, control]))
COMPARISONS = [("Cont_labored", "Cont_unlabored"), ("PTB", "Cont_PTB"), ("DM1", "Cont_labored"),
               ("FGR", "Cont_labored"), ("PrePE", "Cont_labored"), ("TermPE", "Cont_labored"),
               ("DM1", "Cont_unlabored"), ("FGR", "Cont_unlabored"), ("PrePE", "Cont_unlabored"),
               ("TermPE", "Cont_unlabored"), ("PTB", "Cont_labored"), ("PrePE", "Cont_PTB"),
               ("PTB", "Cont_unlabored")]

freq = pd.read_csv(DA_dir+"AllMyeloid_frequencies_per_pt.csv")
stats = pd.read_csv(DA_dir+"AllMyeloid_comparisons_stats.csv")
pt_cond = pd.read_csv(PATIENT_CONDITION_FILE)

cell_cols = [c for c in freq.columns if c != "Patient_ID"]
long_freq = pt_cond.merge(freq, on="Patient_ID", how="left")


# Standard significance stars: *** p<0.001, ** p<0.01, * p<0.05 (callers already filter to p<=0.05).
def pval_to_asterisks(pval):
    if pval < 0.001:
        return "***"
    if pval < 0.01:
        return "**"
    return "*"


# Draws a stacked significance bracket above the plot for every (case, control) pair in sig_pairs.
# Pairs are already filtered to significant (p<=0.05) by significant_pairs_from_stats /
# significant_pairs_fresh, so the label is just the star tier, not the raw p-value.
def add_significance_brackets(ax, sub_df, y_col, sig_pairs, lw=1):
    if not sig_pairs:
        return
    y_max = sub_df[y_col].max()
    y_range = y_max - sub_df[y_col].min()
    step = max(y_range, 1e-6) * 0.1
    for i, (case, control, pval) in enumerate(sig_pairs):
        x1, x2 = CONDITION_ORDER.index(case), CONDITION_ORDER.index(control)
        y = y_max + step * (i + 1)
        ax.plot([x1, x1, x2, x2], [y, y + step * 0.2, y + step * 0.2, y], lw=lw, c="black")
        ax.text((x1 + x2) / 2, y + step * 0.2, pval_to_asterisks(pval), ha="center", va="center", fontsize=9)
    ax.set_ylim(top=y_max + step * (len(sig_pairs) + 1.5))


def significant_pairs_from_stats(cell, threshold=0.05):
    sig = []
    for case, control in COMPARISONS:
        row = stats.loc[(stats["Cluster"] == cell) & (stats["Comparison"] == f"{case}-{control}")]
        if row.empty or pd.isna(row["P-val"].iloc[0]) or row["P-val"].iloc[0] > threshold:
            continue
        sig.append((case, control, row["P-val"].iloc[0]))
    return sig


def significant_pairs_fresh(cell, threshold=0.05):
    sig = []
    for case, control in COMPARISONS:
        case_vals = long_freq.loc[long_freq["Condition_group"] == case, cell].dropna()
        cont_vals = long_freq.loc[long_freq["Condition_group"] == control, cell].dropna()
        if len(case_vals) < 2 or len(cont_vals) < 2:
            continue
        pval = scipy.stats.mannwhitneyu(case_vals, cont_vals).pvalue
        if pval <= threshold:
            sig.append((case, control, pval))
    return sig


# Draws one cell type's box+stripplot onto an existing ax - shared by the HBC/mMac summary
# figure and the all-in-one grid figure (boxplot_grid) below, so both stay visually consistent.
def draw_boxplot(ax, cell, sig_pairs, show_legend=True, fontsize=10, dot_size=4, show_xlabels=True,
                  title_fontsize=None, title_bold=False, box_linewidth=None, ytick_fontsize=None,
                  ylabel_text=None, bracket_lw=1):
    sub = long_freq[["Patient_ID", "Condition_group", cell]].dropna(subset=[cell]).copy()

    box_kwargs = {} if box_linewidth is None else {"linewidth": box_linewidth}
    sns.boxplot(data=sub, x="Condition_group", y=cell, hue="Condition_group", width=0.75, fliersize=0,
                ax=ax, palette=condition_colors, order=CONDITION_ORDER, hue_order=CONDITION_ORDER, **box_kwargs)
    sns.stripplot(data=sub, x="Condition_group", y=cell, size=dot_size, linewidth=0.25, ax=ax, color="black",
                  order=CONDITION_ORDER)
    add_significance_brackets(ax, sub, cell, sig_pairs, lw=bracket_lw)

    ax.set_ylabel(ylabel_text or "Frequency of " + cell + " (%)", fontsize=fontsize)
    ax.set_xlabel("")
    ax.tick_params(axis="x", labelrotation=90, labelsize=fontsize, labelbottom=show_xlabels)
    if ytick_fontsize is not None:
        ax.tick_params(axis="y", labelsize=ytick_fontsize)
    ax.set_title(cell, fontsize=title_fontsize or fontsize, fontweight="bold" if title_bold else "normal")
    if not show_legend and ax.get_legend() is not None:
        ax.get_legend().remove()


A4_PORTRAIT = (8.27, 11.69)  # inches

# All cell types as subplots in one big grid figure, so they can all be compared side by side
# at a glance. Page is a fixed A4 portrait sheet; margin reserves a 10% border on every side so
# the grid itself only covers 80% of the page (width and height).
def boxplot_grid(cells, filename, sig_pairs_fn=significant_pairs_from_stats, n_rows=7, n_cols=5,
                  figsize=A4_PORTRAIT, margin=0.1):
    assert len(cells) <= n_rows * n_cols, f"{len(cells)} cell types don't fit in a {n_rows}x{n_cols} grid"
    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)
    axes = axes.flatten()

    # Last filled row per column, so only the bottom-most panel in each column keeps x labels
    # (the grid's last row may be partially filled, so that isn't always n_rows - 1).
    bottom_row = {}
    for i in range(len(cells)):
        row, col = divmod(i, n_cols)
        bottom_row[col] = row

    for i, (ax, cell) in enumerate(zip(axes, cells)):
        row, col = divmod(i, n_cols)
        draw_boxplot(ax, cell, sig_pairs_fn(cell), show_legend=False, fontsize=7, dot_size=2,
                     show_xlabels=(row == bottom_row[col]), title_fontsize=9, title_bold=True,
                     box_linewidth=0.75, ytick_fontsize=5, ylabel_text="Frequency", bracket_lw=0.5)
    for ax in axes[len(cells):]:  # hide unused panels when len(cells) doesn't fill the grid
        ax.axis("off")
    fig.subplots_adjust(left=margin, right=1 - margin, top=1 - margin, bottom=margin, wspace=0.6, hspace=0.7)
    fig.savefig(filename)
    plt.close(fig)


# 1) Every individual myeloid cell type,all 8 condition groups, every donor as a dot
#    Combined into one big grid figure so all cell types sit side by side
boxplot_grid(cell_cols, "AllMyeloid_boxplots_grid.pdf")

# 2) Lineage-level summary: total HBC vs total mMac fraction of myeloid cells, to highlight the
# opposite trend during complications (HBC down, mMac up) that isn't obvious subset-by-subset
hbc_cols = [c for c in cell_cols if c.startswith("HBC")]
mMac_cols = [c for c in cell_cols if c.startswith("mMac")]
long_freq["HBC"] = long_freq[hbc_cols].sum(axis=1)
long_freq["mMac"] = long_freq[mMac_cols].sum(axis=1)

fig, axes = plt.subplots(1, 2, figsize=(8, 4))
for ax, cell in zip(axes, ["HBC", "mMac"]):
    draw_boxplot(ax, cell, significant_pairs_fresh(cell))
y_top = max(ax.get_ylim()[1] for ax in axes)
y_bottom = min(ax.get_ylim()[0] for ax in axes)
for ax in axes:
    ax.set_ylim(y_bottom, y_top)
plt.tight_layout()
fig.savefig("AllMyeloid_boxplots_HBC_vs_mMac_summary.pdf")
plt.close(fig)
