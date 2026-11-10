import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scanpy as sc
from scipy import sparse
# Read amrkers
def read_markers_dotplot(file_path, lineage):
    df = pd.read_excel(file_path)
    df_lin = df.loc[df["Lineage"] == lineage,].copy()
    strip_cat = pd.DataFrame({"markers": df_lin["Gene/ Protein name"],
                                    "types": df_lin["Type"],
                                    "blocks": df_lin["Block"]})
    return(strip_cat)

def dot_plot_with_umap_stack(
    rna_anndata,
    prot_anndata=None,
    groupby="Final_cluster",
    umap_title=None,
    append_total_cells_to_umap_title=True,
    umap_title_fontsize=14,
    umap_title_fontweight="bold",
    strip_categories=None,
    horizontal_strip=True,
    size_exponent=1.3,
    save_fig=None,
    figsize=(15, 5),
    legend_loc=[.97, 0.16, 0.2, 0.1],
    cbar_loc=[.995, 0.3, 0.015, 0.61],
    cbar2_loc=[.997, 0.3, 0.015, 0.61],
    umap_loc=[-.375, .15, .8, .8],
    standard_scale=False,   # allowed: "var", "group", False/None
    title="average of scaled\nlog10CPM",
    panel_width=1,
    ytick_fontsize=10,
    cluster_number="Cluster number",
    shown_categories=None,
    legend_fontsize=10,
    legend_fontoutline=2,
    umap_dotsize=5,
    xtick_fontsize=8,
    category_size=5,
    legend_dotsize=40,
    legend_length=60,
    dpi = 200,
    umap_aspect=.8,
    umap_legend_loc="on data",
    cbar_orientation="vertical",   # "vertical" or "horizontal"
):
    """
    Dotplot with UMAP on the left and stacked feature blocks on the right.

    strip_categories must be a DataFrame with columns:
        - markers
        - types   (values: "RNA" or "Protein")
        - blocks
    """

    def _to_dense_array(X):
        if sparse.issparse(X):
            return X.toarray()
        return np.asarray(X)

    def _minmax_scale_df(df, axis):
        # axis=0 => scale each column
        # axis=1 => scale each row
        if axis == 0:
            out = df.copy()
            for col in out.columns:
                col_min = out[col].min()
                col_max = out[col].max()
                if col_max == col_min:
                    out[col] = 0.0
                else:
                    out[col] = (out[col] - col_min) / (col_max - col_min)
            return out
        elif axis == 1:
            out = df.copy()
            row_min = out.min(axis=1)
            row_max = out.max(axis=1)
            denom = (row_max - row_min).replace(0, np.nan)
            out = out.sub(row_min, axis=0).div(denom, axis=0).fillna(0)
            return out
        else:
            raise ValueError("axis must be 0 or 1")

    def get_mean_pct(anndata, cluster, standard_scale=None, expression_cutoff=0):
        if cluster not in anndata.obs.columns:
            raise ValueError(f"'{cluster}' not found in anndata.obs")

        if not isinstance(anndata.obs[cluster].dtype, pd.CategoricalDtype):
            anndata = anndata.copy()
            anndata.obs[cluster] = pd.Categorical(anndata.obs[cluster])

        cats = list(anndata.obs[cluster].cat.categories)

        X = _to_dense_array(anndata.X)
        expr_df = pd.DataFrame(X, index=anndata.obs_names, columns=anndata.var_names)

        # mean expression by cluster
        mean_df = expr_df.groupby(anndata.obs[cluster], observed=False).mean()
        mean_df = mean_df.reindex(cats)

        if standard_scale == "var":
            mean_df = _minmax_scale_df(mean_df, axis=0)
        elif standard_scale == "group":
            mean_df = _minmax_scale_df(mean_df, axis=1)
        elif standard_scale in [False, None]:
            pass
        else:
            raise ValueError("standard_scale must be one of: 'var', 'group', False, None")

        pbdt = mean_df.unstack().reset_index()
        pbdt.columns = ["Gene", cluster, "mean"]

        # percent expressing by cluster
        pct_df = (expr_df > expression_cutoff).groupby(anndata.obs[cluster], observed=False).mean() * 100
        pct_df = pct_df.reindex(cats)

        pct = pct_df.unstack().reset_index()
        pct.columns = ["Gene", cluster, "pct_raw"]

        # merge
        pbdt = pbdt.merge(pct, on=["Gene", cluster], how="left")

        # bin pct values for dot size
        bins = [0, 1, 10, 20, 40, 60, 80, 100]
        pbdt["pct"] = pd.cut(
            pbdt["pct_raw"].clip(lower=0, upper=100),
            bins=bins,
            include_lowest=True
        ).apply(lambda x: x.left if pd.notnull(x) else 0)

        # ordering
        pbdt[cluster] = pd.Categorical(
            pbdt[cluster],
            categories=cats[::-1],
            ordered=True
        )
        pbdt["Gene"] = pd.Categorical(
            pbdt["Gene"],
            categories=list(anndata.var_names),
            ordered=True
        )

        pbdt = pbdt.sort_values([cluster, "Gene"]).copy()
        return pbdt

    # ----------------------------
    # Input validation
    # ----------------------------
    if rna_anndata is None and prot_anndata is None:
        raise ValueError("At least one of rna_anndata or prot_anndata must be provided.")

    if strip_categories is None:
        raise ValueError(
            "strip_categories must be provided as a DataFrame with columns "
            "['markers', 'types', 'blocks']"
        )

    required_cols = {"markers", "types", "blocks"}
    if not required_cols.issubset(strip_categories.columns):
        raise ValueError(
            f"strip_categories must contain columns {required_cols}, "
            f"but has {list(strip_categories.columns)}"
        )

    strip_categories = strip_categories.copy()

    if strip_categories[["markers", "types"]].duplicated().any():
        raise ValueError("Duplicate marker/type pairs found in strip_categories")

    valid_types = {"Protein", "RNA"}
    if not set(strip_categories["types"].unique()).issubset(valid_types):
        raise ValueError("strip_categories['types'] must only contain 'Protein' and/or 'RNA'")

    if cbar_orientation not in {"vertical", "horizontal"}:
        raise ValueError("cbar_orientation must be either 'vertical' or 'horizontal'")

    # ----------------------------
    # Pick object used for categories/UMAP defaults
    # ----------------------------
    base_anndata = rna_anndata if rna_anndata is not None else prot_anndata

    if groupby not in base_anndata.obs.columns:
        raise ValueError(f"'{groupby}' not found in anndata.obs")

    # Ensure categorical
    if rna_anndata is not None and not isinstance(rna_anndata.obs[groupby].dtype, pd.CategoricalDtype):
        rna_anndata = rna_anndata.copy()
        rna_anndata.obs[groupby] = pd.Categorical(rna_anndata.obs[groupby])

    if prot_anndata is not None and not isinstance(prot_anndata.obs[groupby].dtype, pd.CategoricalDtype):
        prot_anndata = prot_anndata.copy()
        prot_anndata.obs[groupby] = pd.Categorical(prot_anndata.obs[groupby])

    if not isinstance(base_anndata.obs[groupby].dtype, pd.CategoricalDtype):
        base_anndata = base_anndata.copy()
        base_anndata.obs[groupby] = pd.Categorical(base_anndata.obs[groupby])

    full_categories = list(base_anndata.obs[groupby].cat.categories)

    # ----------------------------
    # shown_categories handling
    # ----------------------------
    if shown_categories is None:
        categories_to_plot = full_categories
        umap_anndata = base_anndata.copy()
    else:
        shown_categories = list(shown_categories)

        bad = [x for x in shown_categories if x not in full_categories]
        if len(bad) > 0:
            raise ValueError(f"These shown_categories are not in {groupby}: {bad}")

        categories_to_plot = shown_categories
        print(f"dropping {len([x for x in full_categories if x not in shown_categories])} categories")

        umap_anndata = base_anndata.copy()

        if rna_anndata is not None:
            rna_anndata = rna_anndata[rna_anndata.obs[groupby].isin(categories_to_plot)].copy()
            rna_anndata.obs[groupby] = pd.Categorical(
                rna_anndata.obs[groupby],
                categories=categories_to_plot,
                ordered=True
            )

        if prot_anndata is not None:
            prot_anndata = prot_anndata[prot_anndata.obs[groupby].isin(categories_to_plot)].copy()
            prot_anndata.obs[groupby] = pd.Categorical(
                prot_anndata.obs[groupby],
                categories=categories_to_plot,
                ordered=True
            )

    # also reorder base object if needed
    anndata_to_use = rna_anndata if rna_anndata is not None else prot_anndata
    if anndata_to_use is None:
        raise ValueError("No data left after filtering shown_categories.")

    if not isinstance(anndata_to_use.obs[groupby].dtype, pd.CategoricalDtype):
        anndata_to_use = anndata_to_use.copy()
        anndata_to_use.obs[groupby] = pd.Categorical(anndata_to_use.obs[groupby])

    n_row = len(categories_to_plot)

    # ----------------------------
    # Marker extraction
    # ----------------------------
    pb_list = []

    if prot_anndata is not None:
        protein = strip_categories.loc[strip_categories["types"] == "Protein", "markers"]
        if len(protein) > 0:
            missing_prot = protein[~protein.isin(prot_anndata.var_names)]
            if len(missing_prot) > 0:
                raise ValueError(
                    f"These protein markers are missing from prot_anndata.var_names: {list(missing_prot)}"
                )

            prot_pb = get_mean_pct(
                prot_anndata[:, protein.tolist()].copy(),
                groupby,
                standard_scale=standard_scale
            )
            prot_pb["feature"] = "protein"
            prot_pb["block"] = prot_pb["Gene"].astype(str).replace(
                dict(zip(
                    strip_categories.loc[strip_categories["types"] == "Protein", "markers"],
                    strip_categories.loc[strip_categories["types"] == "Protein", "blocks"]
                ))
            )
            pb_list.append(prot_pb)

    if rna_anndata is not None:
        genes = strip_categories.loc[strip_categories["types"] == "RNA", "markers"]
        if len(genes) > 0:
            missing_rna = genes[~genes.isin(rna_anndata.var_names)]
            if len(missing_rna) > 0:
                raise ValueError(
                    f"These RNA markers are missing from rna_anndata.var_names: {list(missing_rna)}"
                )

            rna_pb = get_mean_pct(
                rna_anndata[:, genes.tolist()].copy(),
                groupby,
                standard_scale=standard_scale
            )
            rna_pb["feature"] = "RNA"
            rna_pb["block"] = rna_pb["Gene"].astype(str).replace(
                dict(zip(
                    strip_categories.loc[strip_categories["types"] == "RNA", "markers"],
                    strip_categories.loc[strip_categories["types"] == "RNA", "blocks"]
                ))
            )
            pb_list.append(rna_pb)

    if len(pb_list) == 0:
        raise ValueError("No markers matched the provided AnnData objects.")

    pb_all = pd.concat(pb_list, axis=0, ignore_index=True)

    # ----------------------------
    # Shared color ranges across all RNA / protein markers
    # ----------------------------
    rna_vmin = None
    rna_vmax = None
    prot_vmin = None
    prot_vmax = None

    if (pb_all["feature"] == "RNA").any():
        rna_vmin = pb_all.loc[pb_all["feature"] == "RNA", "mean"].min()
        rna_vmax = pb_all.loc[pb_all["feature"] == "RNA", "mean"].max()

    if (pb_all["feature"] == "protein").any():
        prot_vmin = pb_all.loc[pb_all["feature"] == "protein", "mean"].min()
        prot_vmax = pb_all.loc[pb_all["feature"] == "protein", "mean"].max()

    # ----------------------------
    # Colors
    # ----------------------------
    groupby_color_key = groupby + "_colors"
    cluster_color_key = cluster_number + "_colors"

    if groupby_color_key not in anndata_to_use.uns and cluster_color_key not in anndata_to_use.uns:
        raise ValueError(
            f"Could not find '{groupby_color_key}' or '{cluster_color_key}' in anndata.uns"
        )

    color_for_umap = cluster_number if cluster_number in umap_anndata.obs.columns else groupby
    color_key_for_legend = cluster_color_key if cluster_color_key in umap_anndata.uns else groupby_color_key
    legend_colors = umap_anndata.uns[color_key_for_legend]
    if len(legend_colors) < len(full_categories):
        legend_colors = (legend_colors * (len(full_categories) // len(legend_colors) + 1))[:len(full_categories)]

    umap_palette = None
    if shown_categories is not None and color_for_umap == groupby:
        shown_categories_set = set(categories_to_plot)
        umap_palette = [
            color if cat in shown_categories_set else "#D3D3D3"
            for cat, color in zip(full_categories, legend_colors)
        ]

    # ----------------------------
    # Block order
    # ----------------------------
    block_order = strip_categories["blocks"].drop_duplicates().tolist()
    block_length = strip_categories["blocks"].value_counts().reindex(block_order)

    # ----------------------------
    # Plot setup
    # ----------------------------
    sns.set_style("white")
    fig, ax = plt.subplots(
        1,
        len(block_length) + 1,
        figsize=figsize,
        width_ratios=np.hstack([legend_length, (block_length.values * 6)])
    )
    if len(block_length) + 1 == 2:
        ax = np.array(ax).reshape(-1)

    if cbar_orientation == "horizontal":
        plt.subplots_adjust(wspace=.05, bottom=0.30)
    else:
        plt.subplots_adjust(wspace=.05)

    margin = 0.5

    # ----------------------------
    # UMAP
    # ----------------------------
    umap_ax = fig.add_axes(umap_loc, aspect="equal")
    sc.pl.umap(
        umap_anndata,
        color=color_for_umap,
        palette=umap_palette,
        size=umap_dotsize,
        return_fig=False,
        show=False,
        ax=umap_ax,
        legend_loc=umap_legend_loc,
        legend_fontsize=legend_fontsize,
        legend_fontoutline=legend_fontoutline
    )
    if (
        shown_categories is not None
        and color_for_umap == groupby
        and umap_legend_loc == "on data"
    ):
        shown_categories_set = set(categories_to_plot)
        for text in list(umap_ax.texts):
            if text.get_text() not in shown_categories_set:
                text.remove()

    base_umap_title = umap_title if umap_title is not None else str(color_for_umap)
    if append_total_cells_to_umap_title:
        total_cells = rna_anndata.n_obs if rna_anndata is not None else base_anndata.n_obs
        umap_ax.set_title(
            f"{base_umap_title} - {total_cells} total cells",
            fontsize=umap_title_fontsize,
            fontweight=umap_title_fontweight
        )
    else:
        umap_ax.set_title(
            base_umap_title,
            fontsize=umap_title_fontsize,
            fontweight=umap_title_fontweight
        )
    umap_ax.set_aspect(aspect=umap_aspect, adjustable="box")

    cat_to_color = dict(zip(full_categories, legend_colors))
    plot_colors = [cat_to_color[c] for c in categories_to_plot[::-1]]

    sns.scatterplot(
        x=np.repeat(1, len(categories_to_plot)),
        y=pd.Series(np.arange(0, len(categories_to_plot)) + margin),
        c=plot_colors,
        s=legend_dotsize,
        linewidth=.3,
        edgecolor="black",
        ax=ax[0],
        clip_on=False
    )

    ax[0].set_xticks(np.arange(-2, legend_length, 1))
    ax[0].set_xlim(-2, legend_length)
    ax[0].set_ylim(-margin, n_row + margin)
    ax[0].set_yticks(np.arange(0, n_row, 1))

    for k, cat in enumerate(categories_to_plot[::-1]):
        ax[0].text(
            5, k + margin, str(cat),
            horizontalalignment="left",
            verticalalignment="center",
            size=category_size
        )

    ax[0].set_yticklabels([])
    ax[0].set_xticklabels([])
    ax[0].spines['top'].set_visible(False)
    ax[0].spines['right'].set_visible(False)
    ax[0].spines['left'].set_visible(False)
    ax[0].spines['bottom'].set_visible(False)

    if horizontal_strip:
        yticks = ax[0].get_yticks()
        if len(yticks) > 0:
            yticks2 = np.append(yticks, yticks[-1] + 1)
            for y0, y1 in zip(yticks2[::2], np.append(yticks2[1::2], n_row + margin)):
                ax[0].axhspan(y0, y1, color='black', xmax=1.04, alpha=0.1, zorder=0, clip_on=False)
        ax[0].set_ylim(-1 * margin, n_row + margin)
        ax[0].set_yticks(np.arange(0, n_row, 1))

    # ----------------------------
    # Colorbars
    # ----------------------------
    cbar_ax = None
    cbar_ax2 = None

    has_prot = (prot_anndata is not None) and (strip_categories["types"] == "Protein").any()
    has_rna = (rna_anndata is not None) and (strip_categories["types"] == "RNA").any()

    if has_prot and has_rna:
        cbar_ax = fig.add_axes(cbar_loc)
        cbar_ax2 = fig.add_axes(cbar2_loc)
    elif has_rna:
        cbar_ax2 = fig.add_axes(cbar_loc)
    elif has_prot:
        cbar_ax = fig.add_axes(cbar_loc)

    # ----------------------------
    # Dot size legend
    # ----------------------------
    legend_scale = 0.8
    legend_ax = fig.add_axes(legend_loc)
    legend = pd.DataFrame({
        "size": [100, 80, 60, 40, 20],
        "x": np.arange(1, 6, 1) ** legend_scale,
        "y": np.repeat(1, 5)
    })

    legend_ax.scatter(
        legend["x"],
        legend["y"],
        s=legend["size"] ** size_exponent,
        c="grey",
        edgecolors="black",
        linewidth=.5
    )
    legend_ax.tick_params(which='major', length=7)
    legend_ax.set_xticks(np.arange(1, 6, 1) ** legend_scale)
    legend_ax.set_xlim([-1 * margin, 4])
    legend_ax.set_yticklabels([])
    legend_ax.set_xticklabels(legend["size"], size=8)
    legend_ax.spines['top'].set_visible(False)
    legend_ax.spines['right'].set_visible(False)
    legend_ax.spines['left'].set_visible(False)
    legend_ax.spines['bottom'].set_visible(False)
    legend_ax.set_title("Fraction in group\nexpress genes (%)", loc="center", size=10)

    # ----------------------------
    # Dotplots by block
    # ----------------------------
    first_rna_im = None
    first_prot_im = None

    for k, block_name in enumerate(block_order):
        block_pb = pb_all.loc[pb_all["block"] == block_name].copy()

        gene_order = strip_categories.loc[
            strip_categories["blocks"] == block_name, "markers"
        ].tolist()

        # remove "G:" and "P:" prefixes for display only
        gene_labels = [g.replace("G:", "").replace("P:", "") for g in gene_order]

        block_pb["Gene"] = pd.Categorical(
            block_pb["Gene"].astype(str),
            categories=gene_order,
            ordered=True
        )
        block_pb = block_pb.sort_values([groupby, "Gene"]).copy()

        ax[k + 1].set_yticks(np.arange(margin, n_row + margin, 1))
        ax[k + 1].set_ylim(-1 * margin, n_row + margin)

        feature_type = block_pb["feature"].iloc[0]

        if feature_type == "RNA":
            im = ax[k + 1].scatter(
                block_pb["Gene"],
                block_pb[groupby],
                c=block_pb["mean"],
                s=block_pb["pct"].astype(float) ** size_exponent,
                cmap="Reds",
                vmin=rna_vmin,
                vmax=rna_vmax,
                edgecolors="black",
                linewidth=0.3
            )
            if first_rna_im is None:
                first_rna_im = im
        else:
            im = ax[k + 1].scatter(
                block_pb["Gene"],
                block_pb[groupby],
                c=block_pb["mean"],
                s=block_pb["pct"].astype(float) ** size_exponent,
                cmap="Blues",
                vmin=prot_vmin,
                vmax=prot_vmax,
                edgecolors="black",
                linewidth=0.3
            )
            if first_prot_im is None:
                first_prot_im = im

        if horizontal_strip:
            yticks = ax[k + 1].get_yticks()
            if len(yticks) > 0:
                yticks2 = np.append(yticks, yticks[-1] + 1)
                for y0, y1 in zip(yticks2[::2] - 1, np.append(yticks2[1::2] - 1, n_row)):
                    ax[k + 1].axhspan(y0, y1, color='black', alpha=0.1, zorder=0)
                ax[k + 1].set_yticks(yticks)

        ax[k + 1].set_xticks(np.arange(0, len(gene_order), 1))
        ax[k + 1].set_xticklabels(gene_labels, rotation=90, size=xtick_fontsize) # modified
        ax[k + 1].set_yticklabels([], rotation=90)
        ax[k + 1].set_xlim(ax[k + 1].get_xlim()[0] - .5, ax[k + 1].get_xlim()[1] + .5)
        ax[k + 1].set_ylim(-1, n_row)

    # ----------------------------
    # Attach colorbars
    # ----------------------------
    prot_cbar = None
    rna_cbar = None

    if cbar_ax2 is not None and first_rna_im is not None:
        rna_cbar = fig.colorbar(first_rna_im, cax=cbar_ax2, orientation=cbar_orientation)

    if cbar_ax is not None and first_prot_im is not None:
        prot_cbar = fig.colorbar(first_prot_im, cax=cbar_ax, orientation=cbar_orientation)

    # ----------------------------
    # Colorbar text formatting
    # ----------------------------
    if cbar_orientation == "vertical":
        if prot_cbar is not None and rna_cbar is not None:
            cbar_ax.set_title(title, pad=12, fontsize=10)
            cbar_ax.set_xlabel("CITE", fontsize=10, labelpad=4)
            cbar_ax2.set_xlabel("RNA", fontsize=10, labelpad=4)
        elif rna_cbar is not None:
            cbar_ax2.set_title(title, pad=12, fontsize=10)
            cbar_ax2.set_xlabel("RNA", fontsize=10, labelpad=4)
        elif prot_cbar is not None:
            cbar_ax.set_title(title, pad=12, fontsize=10)
            cbar_ax.set_xlabel("CITE", fontsize=10, labelpad=4)

    else:  # horizontal
        if prot_cbar is not None:
            prot_cbar.set_label("CITE", fontsize=10, labelpad=2)
            prot_cbar.ax.xaxis.set_label_position("bottom")
            prot_cbar.ax.tick_params(axis='x', labelsize=9, length=3)

        if rna_cbar is not None:
            rna_cbar.set_label("RNA", fontsize=10, labelpad=2)
            rna_cbar.ax.xaxis.set_label_position("bottom")
            rna_cbar.ax.tick_params(axis='x', labelsize=9, length=3)

        if cbar_ax is not None:
            cbar_ax.text(
                0.5, 1.65, title,
                transform=cbar_ax.transAxes,
                ha="center",
                va="bottom",
                fontsize=10
            )

        if cbar_ax2 is not None:
            cbar_ax2.text(
                0.5, 1.65, title,
                transform=cbar_ax2.transAxes,
                ha="center",
                va="bottom",
                fontsize=10
            )

    # ----------------------------
    # Panel borders
    # ----------------------------
    for axis_name in ['top', 'bottom', 'left', 'right']:
        for a in ax:
            a.spines[axis_name].set_linewidth(panel_width)

        if cbar_ax is not None:
            cbar_ax.spines[axis_name].set_linewidth(panel_width)
        if cbar_ax2 is not None:
            cbar_ax2.spines[axis_name].set_linewidth(panel_width)

    if save_fig is not None:
        fig.savefig(save_fig, bbox_inches='tight', dpi=dpi, format='pdf')

    return fig, ax
