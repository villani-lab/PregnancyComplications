import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def _padj_to_ast(padj):
    """Convert adjusted p-values to significance stars."""
    if pd.isna(padj):
        return ""
    if padj < 0.001:
        return "***"
    if padj < 0.01:
        return "**"
    if padj < 0.05:
        return "*"
    return "n.s."


def _infer_comparison_label(filepath: str) -> str:
    """Infer a readable comparison label from the DEG workbook filename."""
    base = os.path.splitext(os.path.basename(filepath))[0]
    prefixes = (
        "DEGs_HBC_V5_Level2_",
        "DEGs_PAMM_V5_Level1.5_",
        "DEGs_Mon_V5_Level1.5_",
        "DEGs_DC_V5_Level1.5_",
        "DEGs_Troph_V6_Level3_",
        "DEGs_Endo_",
    )
    for prefix in prefixes:
        if base.startswith(prefix):
            return base[len(prefix):]
    return base


def build_gene_stats_table(file_map_by_condition,
                           genes_orig,
                           condition_order=None,
                           condition_labels=None,
                           lfc_col="log2FoldChange",
                           padj_col="padj",
                           include_missing=False,
                           save_path=None):
    """Build a long-format gene statistics table from DEG workbooks.

    Parameters mirror the paired dotplot inputs so the same setup can be reused.
    Returns a DataFrame with one row per condition x celltype x gene.
    """
    if condition_order is None:
        condition_order = list(file_map_by_condition.keys())
    if condition_labels is None:
        condition_labels = {}
    if len(condition_order) == 0:
        raise ValueError("condition_order must contain at least one condition.")

    grouped_input = hasattr(genes_orig, "values") and hasattr(genes_orig, "keys")
    if grouped_input:
        gene_group_lookup = {}
        genes_list = []
        for group_name, group_genes in genes_orig.items():
            for gene in group_genes:
                genes_list.append(gene)
                gene_group_lookup[gene] = group_name
    else:
        genes_list = list(genes_orig)
        gene_group_lookup = {gene: None for gene in genes_list}

    genes_prefixed = ["G:" + g for g in genes_list]
    convert_clusters = {'HBC_4_IL1B_MMP19': 'HBC_4', 'mMac_upper': 'mMac_4_6_7'}
    rows = []

    for condition in condition_order:
        file_map = file_map_by_condition[condition]
        condition_label = condition_labels.get(condition, condition)

        for original_celltype, filepath in file_map.items():
            celltype = convert_clusters.get(original_celltype, original_celltype)

            if filepath is None or not os.path.exists(filepath):
                if include_missing:
                    for gene in genes_list:
                        rows.append({
                            "condition": condition_label,
                            "celltype": celltype,
                            "gene_group": gene_group_lookup.get(gene),
                            "gene_name": gene,
                            "log2fc": np.nan,
                            "padj": np.nan,
                            "ast": "",
                            "comparison": condition_label,
                            "source_file": filepath,
                        })
                continue

            try:
                df_sheet = read_celltype_sheet(filepath, sheet=original_celltype)
                sub = extract_gene_stats(df_sheet, genes_prefixed, lfc_col, padj_col)
            except Exception as e:
                print(f"[WARN] Could not read stats for '{original_celltype}' from {filepath}: {e}")
                if include_missing:
                    for gene in genes_list:
                        rows.append({
                            "condition": condition_label,
                            "celltype": celltype,
                            "gene_group": gene_group_lookup.get(gene),
                            "gene_name": gene,
                            "log2fc": np.nan,
                            "padj": np.nan,
                            "ast": "",
                            "comparison": condition_label,
                            "source_file": filepath,
                        })
                continue

            comparison_label = _infer_comparison_label(filepath)
            for gene, gene_prefixed in zip(genes_list, genes_prefixed):
                if gene_prefixed in sub.index:
                    log2fc = float(sub.loc[gene_prefixed, lfc_col])
                    padj = float(sub.loc[gene_prefixed, padj_col])
                else:
                    log2fc = np.nan
                    padj = np.nan
                    if not include_missing:
                        continue

                rows.append({
                    "condition": condition_label,
                    "celltype": celltype,
                    "gene_group": gene_group_lookup.get(gene),
                    "gene_name": gene,
                    "log2fc": log2fc,
                    "padj": padj,
                    "ast": _padj_to_ast(padj),
                    "comparison": comparison_label,
                    "source_file": filepath,
                })

    out = pd.DataFrame(rows)
    if not out.empty:
        out["condition"] = pd.Categorical(
            out["condition"],
            categories=[condition_labels.get(c, c) for c in condition_order],
            ordered=True,
        )
        out["gene_name"] = pd.Categorical(out["gene_name"], categories=genes_list, ordered=True)
        if grouped_input:
            group_order = list(genes_orig.keys())
            out["gene_group"] = pd.Categorical(out["gene_group"], categories=group_order, ordered=True)
        out = out.sort_values(["condition", "celltype", "gene_group", "gene_name"], na_position="last").reset_index(drop=True)

    if save_path is not None:
        if save_path.lower().endswith(".xlsx"):
            out.to_excel(save_path, index=False)
        else:
            out.to_csv(save_path, index=False)

    return out
# Lineage separators config
# You can adjust how lineage is inferred here if needed.
def lineage_of_celltype(ct: str) -> str:
    if ct.startswith("HBC"):
        return "HBC"
    if ct.startswith("mMac"):
        return "mMac"
    if ct.startswith("VCT"):
        return "VCT"
    if ct.startswith("fPV"):
        return "PV"
    if (ct=="mCD4_eff") | (ct=="mCD8_eff_GZMK"):
        return "T"
    return "Other"

# ----------------------------
# Helpers
# ----------------------------
def read_celltype_sheet(filepath: str, sheet: str) -> pd.DataFrame:
    """Read a single sheet (celltype) from an Excel file."""
    return pd.read_excel(filepath, sheet_name=sheet, engine="openpyxl", index_col=0)

def extract_gene_stats(df: pd.DataFrame, genes: list,lfc_col, padj_col) -> pd.DataFrame:
    """Return DF indexed by gene with columns log2FoldChange and padj."""
    df2 = df.copy()
    df2.index = df2.index.astype(str)
    sub = df2.loc[df2.index.intersection(genes), [lfc_col, padj_col]].copy()
    sub.index.name = "gene"

    # Ensure expected columns
    if lfc_col not in sub.columns or padj_col not in sub.columns:
        missing_cols = [c for c in [lfc_col, padj_col] if c not in sub.columns]
        raise ValueError(f"Missing expected columns {missing_cols} in sheet.")
    return sub[[lfc_col, padj_col]]

def padj_to_size(padj_vals: np.ndarray,padj_floor, padj_cap, max_neglog10, size_max, size_min) -> np.ndarray:
    """Map padj -> dot size using -log10(padj), capped & scaled."""
    p = np.array(padj_vals, dtype=float)
    p = np.clip(p, padj_floor, padj_cap)
    neglog = -np.log10(p)
    neglog = np.clip(neglog, 0, max_neglog10)
    # scale 0..max_neglog10 -> size_min..size_max
    return size_min + (neglog / max_neglog10) * (size_max - size_min)


def _collect_dotplot_matrices(file_map, genes, lfc_col="log2FoldChange", padj_col="padj"):
    """Collect log2FC and padj matrices for a set of genes across cell types."""
    celltypes = list(file_map.keys())
    log2fc = pd.DataFrame(index=genes, columns=celltypes, dtype=float)
    padj = pd.DataFrame(index=genes, columns=celltypes, dtype=float)

    for ct in celltypes:
        fp = file_map.get(ct)
        if fp is None:
            print(f"[WARN] No file mapping for celltype: {ct}")
            continue
        if not os.path.exists(fp):
            print(f"[WARN] File not found for {ct}: {fp}")
            continue
        try:
            df_sheet = read_celltype_sheet(fp, sheet=ct)
        except Exception as e:
            print(f"[WARN] Could not read sheet '{ct}' from {fp}: {e}")
            continue
        try:
            sub = extract_gene_stats(df_sheet, genes, lfc_col, padj_col)
        except Exception as e:
            print(f"[WARN] Could not extract stats for '{ct}' from {fp}: {e}")
            continue
        for g in genes:
            if g in sub.index:
                log2fc.loc[g, ct] = float(sub.loc[g, lfc_col])
                padj.loc[g, ct] = float(sub.loc[g, padj_col])

    log2fc.index = [g[2:] for g in log2fc.index]
    padj.index = [g[2:] for g in padj.index]

    convert_clusters = {'HBC_4_IL1B_MMP19': 'HBC_4', 'mMac_upper': 'mMac_4_6_7'}
    log2fc = log2fc.rename(columns=convert_clusters)
    padj = padj.rename(columns=convert_clusters)
    celltypes = [convert_clusters.get(x, x) for x in celltypes]

    return log2fc, padj, celltypes


def plot_FC_dotplot(file_map,
                    gene_groups,
                    title="",
                    pdf_name="dot_plot.pdf",
                    save=True,
                    clip_fc=True,
                    v_clip=2.5,
                    sig_padj_thresh=0.15,
                    sig_abs_lfc_thresh=0.58,
                    show_gene_groups=True,
                    group_label_x=1.08,
                    group_label_ha="left",
                    group_label_fontsize=12,
                    draw_group_separators=True):
    celltypes = list(file_map.keys())
    genes_orig = [g for grp in gene_groups.values() for g in grp]
    genes = ["G:" + g for grp in gene_groups.values() for g in grp]

    lfc_col = "log2FoldChange"
    padj_col = "padj"

    padj_floor = 1e-300
    padj_cap = 1.0
    size_min = 40
    size_max = 600
    max_neglog10 = 20
    sig_alpha = sig_padj_thresh
    sig_edgecolor = "black"
    sig_linewidth = 1.6
    ns_edgecolor = "none"
    ns_linewidth = 0.0

    print("1) Collect DESeq2 stats into matrices")
    log2fc, padj, celltypes = _collect_dotplot_matrices(file_map, genes, lfc_col=lfc_col, padj_col=padj_col)

    print("2) Plot dotplot")
    # choose color scale
    if clip_fc:
        v = float(v_clip)
        vmin, vmax = -v, v
    else:
        v = float(np.nanmax(np.abs(log2fc.values)))
        v = 1.0 if (v == 0 or np.isnan(v)) else v
        vmin, vmax = -v, v

    n_ct = len(celltypes)
    n_g = len(genes_orig)
    fig_w = max(3.0, 0.5 * n_ct + 1)
    print("figure width ",str(fig_w))
    fig_h = max(4.0, 0.35 * n_g + 2.0)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    xs, ys, colors, sizes, padj_vals, fc_vals = [], [], [], [], [], []
    for i, g in enumerate(genes_orig):
        for j, ct in enumerate(celltypes):
            fc = log2fc.loc[g, ct]
            p = padj.loc[g, ct]
            if np.any(np.isnan(fc)) or np.any(np.isnan(p)):
                continue
            xs.append(j)
            ys.append(i)
            colors.append(np.clip(fc, -v_clip, v_clip) if clip_fc else fc)
            sizes.append(p)
            padj_vals.append(p)
            fc_vals.append(fc)

    sizes = padj_to_size(np.array(sizes, dtype=float), padj_floor, padj_cap, max_neglog10, size_max, size_min)
    padj_vals = np.array(padj_vals, dtype=float)
    fc_vals = np.array(fc_vals, dtype=float)

    print("3) Split significant to non-significant")
    sig_mask = (padj_vals < sig_alpha) & (np.abs(fc_vals) > sig_abs_lfc_thresh)
    xs = np.array(xs)
    ys = np.array(ys)
    colors = np.array(colors)

    sca = ax.scatter(xs[~sig_mask], ys[~sig_mask],
                     s=sizes[~sig_mask], c=colors[~sig_mask],
                     cmap="bwr", vmin=vmin, vmax=vmax,
                     edgecolors=ns_edgecolor, linewidths=ns_linewidth)

    ax.scatter(xs[sig_mask], ys[sig_mask],
               s=sizes[sig_mask], c=colors[sig_mask],
               cmap="bwr", vmin=vmin, vmax=vmax,
               edgecolors=sig_edgecolor, linewidths=sig_linewidth)

    ax.set_xticks(range(len(celltypes)))
    ax.set_xticklabels(celltypes, rotation=90, ha="right")
    ax.set_yticks(range(len(genes_orig)))
    ax.set_yticklabels(genes_orig)
    ax.invert_yaxis()
    ax.grid(False)

    ax.set_title(title)

    r_pts = np.sqrt(size_max) / 2.0
    r_in = r_pts / 72.0
    x_pad = (r_in / fig_w) * n_ct * 1.2
    y_pad = (r_in / fig_h) * n_g * 1.2
    ax.set_xlim(-0.5 - x_pad, (n_ct - 1) + 0.5 + x_pad)
    ax.set_ylim((n_g - 1) + 0.5 + y_pad, -0.5 - y_pad)

    if show_gene_groups:
        group_names = list(gene_groups.keys())
        group_sizes = [len(v) for v in gene_groups.values()]
        boundaries = np.cumsum(group_sizes)[:-1]

        if draw_group_separators:
            for b in boundaries:
                ax.axhline(b - 0.5, color="0.75", linewidth=0.8, zorder=0)

        starts = np.r_[0, np.cumsum(group_sizes)[:-1]]
        ends = np.cumsum(group_sizes) - 1
        mids = (starts + ends) / 2.0
        for name, ymid in zip(group_names, mids):
            ax.text(group_label_x, ymid, name, transform=ax.get_yaxis_transform(),
                    va="center", ha=group_label_ha, fontsize=group_label_fontsize,
                    fontweight="bold", clip_on=False)

    legend_padj = [0.15, 0.05, 0.01, 1e-3, 1e-5]
    for p in legend_padj:
        ax.scatter([], [], s=padj_to_size([p], padj_floor, padj_cap, max_neglog10, size_max, size_min)[0],
                   c="gray", label=f"padj={p:g}")
    ax.scatter([], [], s=200, facecolors="white", edgecolors="black", linewidths=sig_linewidth,
               label=f"sig: padj < {sig_padj_thresh:g}, |log2FC| > {sig_abs_lfc_thresh:g}")

    plt.subplots_adjust(left=0.22, right=0.70, top=0.88, bottom=0.25)
    ax_pos = ax.get_position()
    cbar_width = 0.012
    cbar_height = ax_pos.height * 0.26
    cbar_x = ax_pos.x1 + 0.01
    cbar_y = ax_pos.y0 + ax_pos.height * 0.08
    cax = fig.add_axes([cbar_x, cbar_y, cbar_width, cbar_height])
    cbar = fig.colorbar(sca, cax=cax)
    if clip_fc:
        cbar.set_label(f"log2FoldChange (clipped to ±{v_clip:g})")
    else:
        cbar.set_label("log2FoldChange")
    cbar.ax.tick_params(labelsize=9)

    legend_x = (cbar_x + cbar_width + 0.05 - ax_pos.x0) / ax_pos.width
    legend_y = ((ax_pos.y0 + ax_pos.height * 0.02) - ax_pos.y0) / ax_pos.height
    ax.legend(title="Legend", loc="upper left", bbox_to_anchor=(legend_x, legend_y),
              frameon=False, borderaxespad=0.0)

    if save:
        fig.savefig(pdf_name)

    return fig, ax




def plot_FC_dotplot_V2(file_map,
                    genes_orig,
                    title="",
                    pdf_name="dot_plot.pdf",
                    save=True,
                    clip_fc=True,
                    v_clip=2.5,
                    sig_padj_thresh=0.15,
                    sig_abs_lfc_thresh=0.58):
    celltypes = list(file_map.keys())
    genes = ["G:" + g for g in genes_orig]

    lfc_col = "log2FoldChange"
    padj_col = "padj"

    padj_floor = 1e-300
    padj_cap = 1.0
    size_min = 40
    size_max = 600
    max_neglog10 = 20
    sig_alpha = sig_padj_thresh
    sig_edgecolor = "black"
    sig_linewidth = 1.6
    ns_edgecolor = "none"
    ns_linewidth = 0.0

    print("1) Collect DESeq2 stats into matrices")
    log2fc, padj, celltypes = _collect_dotplot_matrices(file_map, genes, lfc_col=lfc_col, padj_col=padj_col)

    print("2) Plot dotplot")
    # choose color scale
    if clip_fc:
        v = float(v_clip)
        vmin, vmax = -v, v
    else:
        v = float(np.nanmax(np.abs(log2fc.values)))
        v = 1.0 if (v == 0 or np.isnan(v)) else v
        vmin, vmax = -v, v

    n_ct = len(celltypes)
    n_g = len(genes_orig)
    fig_w = max(3.0, 0.5 * n_ct + 1)
    print("figure width ",str(fig_w))
    fig_h = max(4.0, 0.35 * n_g + 2.0)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    xs, ys, colors, sizes, padj_vals, fc_vals = [], [], [], [], [], []
    for i, g in enumerate(genes_orig):
        for j, ct in enumerate(celltypes):
            fc = log2fc.loc[g, ct]
            p = padj.loc[g, ct]
            if np.any(np.isnan(fc)) or np.any(np.isnan(p)):
                continue
            xs.append(j)
            ys.append(i)
            colors.append(np.clip(fc, -v_clip, v_clip) if clip_fc else fc)
            sizes.append(p)
            padj_vals.append(p)
            fc_vals.append(fc)

    sizes = padj_to_size(np.array(sizes, dtype=float), padj_floor, padj_cap, max_neglog10, size_max, size_min)
    padj_vals = np.array(padj_vals, dtype=float)
    fc_vals = np.array(fc_vals, dtype=float)

    print("3) Split significant to non-significant")
    sig_mask = (padj_vals < sig_alpha) & (np.abs(fc_vals) > sig_abs_lfc_thresh)
    xs = np.array(xs)
    ys = np.array(ys)
    colors = np.array(colors)

    sca = ax.scatter(xs[~sig_mask], ys[~sig_mask],
                     s=sizes[~sig_mask], c=colors[~sig_mask],
                     cmap="bwr", vmin=vmin, vmax=vmax,
                     edgecolors=ns_edgecolor, linewidths=ns_linewidth)

    ax.scatter(xs[sig_mask], ys[sig_mask],
               s=sizes[sig_mask], c=colors[sig_mask],
               cmap="bwr", vmin=vmin, vmax=vmax,
               edgecolors=sig_edgecolor, linewidths=sig_linewidth)

    ax.set_xticks(range(len(celltypes)))
    ax.set_xticklabels(celltypes, rotation=90, ha="right")
    ax.set_yticks(range(len(genes_orig)))
    ax.set_yticklabels(genes_orig)
    ax.invert_yaxis()
    ax.grid(False)

    ax.set_title(title + "\ncolor = log2FoldChange, size = padj (via -log10(padj))")

    r_pts = np.sqrt(size_max) / 2.0
    r_in = r_pts / 72.0
    x_pad = (r_in / fig_w) * n_ct * 1.2
    y_pad = (r_in / fig_h) * n_g * 1.2
    ax.set_xlim(-0.5 - x_pad, (n_ct - 1) + 0.5 + x_pad)
    ax.set_ylim((n_g - 1) + 0.5 + y_pad, -0.5 - y_pad)

    cbar = plt.colorbar(sca, ax=ax, shrink=0.3, pad=0.02, aspect=25)
    if clip_fc:
        cbar.set_label(f"log2FoldChange (clipped to ±{v_clip:g})")
    else:
        cbar.set_label("log2FoldChange")
    cbar.ax.tick_params(labelsize=9)

    legend_padj = [0.15, 0.05, 0.01, 1e-3, 1e-5]
    for p in legend_padj:
        ax.scatter([], [], s=padj_to_size([p], padj_floor, padj_cap, max_neglog10, size_max, size_min)[0],
                   c="gray", label=f"padj={p:g}")
    ax.scatter([], [], s=200, facecolors="white", edgecolors="black", linewidths=sig_linewidth,
               label=f"sig: padj < {sig_padj_thresh:g}, |log2FC| > {sig_abs_lfc_thresh:g}")
    ax.legend(title="Legend", loc="upper left", bbox_to_anchor=(1.02, 0.00), frameon=False, borderaxespad=0.0)

    plt.subplots_adjust(left=0.28, right=0.78, top=0.88, bottom=0.25)
    if save:
        fig.savefig(pdf_name)

    return fig, ax


def plot_FC_dotplot_paired_conditions(file_map_by_condition,
                                      genes_orig,
                                      title="",
                                      pdf_name="dot_plot.pdf",
                                      save=True,
                                      clip_fc=True,
                                      v_clip=2.5,
                                      sig_padj_thresh=0.15,
                                      sig_abs_lfc_thresh=0.58,
                                      condition_order=None,
                                      condition_labels=None,
                                      x_layout="interleaved",
                                      condition_gap=0.4,
                                      facet_label_y=1.01,
                                      show_title=False,
                                      show_gene_groups=None,
                                      show_condition_headers=True,
                                      group_label_x=1.08,
                                      group_label_ha="left",
                                      group_label_fontsize=12,
                                      draw_group_separators=True,
                                      separator_color="#d3d3d3",
                                      separator_linewidth=1.5):
    """Plot multiple conditions in one dotplot.

    genes_orig can be either:
      - a flat list of gene symbols
      - an ordered dict-like mapping of group name -> list of gene symbols

    x_layout:
      - "interleaved": for each cell type, plot all conditions side by side
      - "condition_blocks": plot all cell types for the first condition, then the next
    """
    if condition_order is None:
        condition_order = list(file_map_by_condition.keys())
    if condition_labels is None:
        condition_labels = {}
    if len(condition_order) == 0:
        raise ValueError("condition_order must contain at least one condition.")

    grouped_input = hasattr(genes_orig, "values") and hasattr(genes_orig, "keys")
    if grouped_input:
        gene_groups = genes_orig
        genes_orig = [g for grp in gene_groups.values() for g in grp]
    else:
        gene_groups = None
        genes_orig = list(genes_orig)

    if show_gene_groups is None:
        show_gene_groups = grouped_input

    genes = ["G:" + g for g in genes_orig]
    lfc_col = "log2FoldChange"
    padj_col = "padj"

    padj_floor = 1e-300
    padj_cap = 1.0
    size_min = 40
    size_max = 600
    max_neglog10 = 20
    sig_alpha = sig_padj_thresh
    sig_edgecolor = "black"
    sig_linewidth = 1.6
    ns_edgecolor = "none"
    ns_linewidth = 0.0

    print("1) Collect DESeq2 stats into matrices for paired dotplot")
    stats_by_condition = {}
    celltype_order = []
    for condition in condition_order:
        log2fc, padj, celltypes = _collect_dotplot_matrices(
            file_map_by_condition[condition], genes, lfc_col=lfc_col, padj_col=padj_col
        )
        stats_by_condition[condition] = {"log2fc": log2fc, "padj": padj, "celltypes": celltypes}
        for ct in celltypes:
            if ct not in celltype_order:
                celltype_order.append(ct)

    if x_layout not in {"interleaved", "condition_blocks"}:
        raise ValueError("x_layout must be either 'interleaved' or 'condition_blocks'.")

    print("2) Build x-axis order")
    x_positions = {}
    xticks = []
    xticklabels = []
    group_centers = []
    n_conditions = len(condition_order)
    separator_positions = []

    if x_layout == "interleaved":
        step = n_conditions + condition_gap

        for cell_idx, ct in enumerate(celltype_order):
            start = cell_idx * step
            group_positions = []
            for cond_idx, condition in enumerate(condition_order):
                pos = start + cond_idx
                x_positions[(ct, condition)] = pos
                xticks.append(pos)
                label = condition_labels.get(condition, condition)
                xticklabels.append(label)
                group_positions.append(pos)
            group_centers.append(np.mean(group_positions))

        for cell_idx in range(len(celltype_order) - 1):
            separator_positions.append(cell_idx * step + (n_conditions - 0.5))
    else:
        current_pos = 0.0
        for cond_idx, condition in enumerate(condition_order):
            stats = stats_by_condition[condition]
            block_positions = []
            for ct in stats["celltypes"]:
                pos = current_pos
                x_positions[(ct, condition)] = pos
                xticks.append(pos)
                xticklabels.append(ct)
                block_positions.append(pos)
                current_pos += 1.0

            if block_positions:
                group_centers.append(np.mean(block_positions))

            if cond_idx < (n_conditions - 1):
                separator_positions.append(current_pos - 0.5 + (condition_gap / 2.0))
                current_pos += condition_gap

    print("3) Plot paired dotplot")
    if clip_fc:
        v = float(v_clip)
        vmin, vmax = -v, v
    else:
        max_abs = []
        for condition in condition_order:
            max_abs.append(np.nanmax(np.abs(stats_by_condition[condition]["log2fc"].values)))
        v = float(np.nanmax(max_abs))
        v = 1.0 if (v == 0 or np.isnan(v)) else v
        vmin, vmax = -v, v

    n_slots = len(xticks)
    n_g = len(genes_orig)
    fig_w = max(4.5, 0.55 * n_slots + 1.5)
    fig_h = max(4.0, 0.35 * n_g + 2.0)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    xs, ys, colors, sizes, padj_vals, fc_vals = [], [], [], [], [], []
    for i, gene in enumerate(genes_orig):
        for condition in condition_order:
            stats = stats_by_condition[condition]
            for ct in stats["celltypes"]:
                fc = stats["log2fc"].loc[gene, ct]
                p = stats["padj"].loc[gene, ct]
                if np.any(np.isnan(fc)) or np.any(np.isnan(p)):
                    continue
                xs.append(x_positions[(ct, condition)])
                ys.append(i)
                colors.append(np.clip(fc, -v_clip, v_clip) if clip_fc else fc)
                sizes.append(p)
                padj_vals.append(p)
                fc_vals.append(fc)

    sizes = padj_to_size(np.array(sizes, dtype=float), padj_floor, padj_cap, max_neglog10, size_max, size_min)
    padj_vals = np.array(padj_vals, dtype=float)
    fc_vals = np.array(fc_vals, dtype=float)
    sig_mask = (padj_vals < sig_alpha) & (np.abs(fc_vals) > sig_abs_lfc_thresh)
    xs = np.array(xs)
    ys = np.array(ys)
    colors = np.array(colors)

    sca = ax.scatter(xs[~sig_mask], ys[~sig_mask],
                     s=sizes[~sig_mask], c=colors[~sig_mask],
                     cmap="bwr", vmin=vmin, vmax=vmax,
                     edgecolors=ns_edgecolor, linewidths=ns_linewidth)

    ax.scatter(xs[sig_mask], ys[sig_mask],
               s=sizes[sig_mask], c=colors[sig_mask],
               cmap="bwr", vmin=vmin, vmax=vmax,
               edgecolors=sig_edgecolor, linewidths=sig_linewidth)

    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels, rotation=90, ha="right")
    ax.set_yticks(range(len(genes_orig)))
    ax.set_yticklabels(genes_orig)
    ax.invert_yaxis()
    ax.grid(False)

    r_pts = np.sqrt(size_max) / 2.0
    r_in = r_pts / 72.0
    x_pad = (r_in / fig_w) * max(n_slots, 1) * 1.2
    y_pad = (r_in / fig_h) * n_g * 1.2
    ax.set_xlim(min(xticks) - 0.5 - x_pad, max(xticks) + 0.5 + x_pad)
    ax.set_ylim((n_g - 1) + 0.5 + y_pad, -0.5 - y_pad)

    if x_layout == "interleaved":
        for center, ct in zip(group_centers, celltype_order):
            ax.text(center, facet_label_y, ct, transform=ax.get_xaxis_transform(),
                    va="bottom", ha="center", fontsize=11, fontweight="bold",
                    clip_on=False)
    elif show_condition_headers:
        for center, condition in zip(group_centers, condition_order):
            label = condition_labels.get(condition, condition)
            ax.text(center, facet_label_y, label, transform=ax.get_xaxis_transform(),
                    va="bottom", ha="center", fontsize=11, fontweight="bold",
                    clip_on=False)

    if show_gene_groups and gene_groups is not None:
        group_names = list(gene_groups.keys())
        group_sizes = [len(v) for v in gene_groups.values()]
        boundaries = np.cumsum(group_sizes)[:-1]

        for b in boundaries:
            ax.axhline(b - 0.5, color="0.75", linewidth=0.8, zorder=0)

        starts = np.r_[0, np.cumsum(group_sizes)[:-1]]
        ends = np.cumsum(group_sizes) - 1
        mids = (starts + ends) / 2.0
        for name, ymid in zip(group_names, mids):
            ax.text(group_label_x, ymid, name, transform=ax.get_yaxis_transform(),
                    va="center", ha=group_label_ha, fontsize=group_label_fontsize,
                    fontweight="bold", clip_on=False)

    if draw_group_separators:
        for separator_x in separator_positions:
            ax.axvline(separator_x, color=separator_color, linewidth=separator_linewidth)

    cbar = plt.colorbar(sca, ax=ax, shrink=0.3, pad=0.02, aspect=25)
    if clip_fc:
        cbar.set_label(f"log2FoldChange (clipped to ±{v_clip:g})")
    else:
        cbar.set_label("log2FoldChange")
    cbar.ax.tick_params(labelsize=9)

    legend_padj = [0.15, 0.05, 0.01, 1e-3, 1e-5]
    for p in legend_padj:
        ax.scatter([], [], s=padj_to_size([p], padj_floor, padj_cap, max_neglog10, size_max, size_min)[0],
                   c="gray", label=f"padj={p:g}")
    ax.scatter([], [], s=200, facecolors="white", edgecolors="black", linewidths=sig_linewidth,
               label=f"sig: padj < {sig_padj_thresh:g}, |log2FC| > {sig_abs_lfc_thresh:g}")
    ax.legend(title="Legend", loc="upper left", bbox_to_anchor=(1.02, 0.00), frameon=False, borderaxespad=0.0)

    if show_title and title:
        fig.suptitle(title, y=0.98)
        top_margin = 0.88
    else:
        top_margin = 0.94
    plt.subplots_adjust(left=0.28, right=0.78, top=top_margin, bottom=0.32)
    if save:
        fig.savefig(pdf_name)

    return fig, ax
