import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import pandas as pd
import numpy as np
import muon as mu
import scanpy as sc
import scipy
import adjustText
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Patch
from pandas import pivot_table
from matplotlib.patches import Patch
import statsmodels
from statannotations.Annotator import Annotator
from scipy import sparse
from matplotlib.ticker import MaxNLocator
import math
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
from adpbulk.adpbulk import ADPBulk
import os
import re
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']

# default lineage palette you provided
_DEF_LINEAGES = ['Trophoblasts','EVT','fPV','fFB','dS','Endothelial','HBC','Myeloid','T','NK','B_Plasma','HSC','Other']
_DEF_COLORS = ["#60B364","#B2DF8A",  # trophoblasts
               "#08306b","#2b8cbe","#c6dbef","#8F6FE8",  # stromal/endo
               "#800026","#f03b20","#fd8d3c","#fed976","#fa9fb5",  # immune
               "grey","lightgrey"]
_DEF_PALETTE = dict(zip(_DEF_LINEAGES, _DEF_COLORS))
sc.set_figure_params(vector_friendly=True, dpi_save=400)



# This function runs DA for each condition group:
#  1. Computes frequencies per patient per cluster
#  2. Computes comparisons between condition groups
#  3. Plots box plots per comparison
#
# mdata - a muon data object that includes the following columns in obs: <cluster_col>, "Sample", "Condition"
# cluster_col - the name of the column in data.obs that contains the lineage/ clusters to plot
# comparisons - a list of lists with the two conditions to compare, the first in the list should be the case and the second would be the control
# samp_condition_file - a file that lists which samples are part of which condition groups
# lineage - PDF prefix for the file in which the plot will be saved
def DA_condition_groups(mdata, cluster_col, comparisons, samp_condition_file, lineage):
    samp_cond = pd.read_csv(samp_condition_file)
    if (not os.path.isfile(lineage + "_comparisons_stats.csv")):
        print("1) Computing frequencies...")
        # Computing frequencies
        df_freq = pd.crosstab(mdata["GEX"].obs.loc[:, "Patient_ID"], mdata["GEX"].obs.loc[:, cluster_col])
        df_freq = df_freq.div(df_freq.sum(axis=1), axis=0) * 100.0
        df_freq.to_csv(lineage + "_frequencies_per_pt.csv")

        print("2) Computing comparisons stats...")
        comp_res = pd.DataFrame()
        for dis_pair in comparisons:  # dis_pair =["PTB","Cont_PTB"]
            pair_res = pd.DataFrame()
            for clust in df_freq.columns:  # clust="1"
                case_samp = samp_cond.loc[samp_cond["Condition_group"] == dis_pair[0], "Patient_ID"].to_list()
                cont_samp = samp_cond.loc[samp_cond["Condition_group"] == dis_pair[1], "Patient_ID"].to_list()
                data_case = df_freq.loc[case_samp, clust]
                data_cont = df_freq.loc[cont_samp, clust]

                pval = scipy.stats.mannwhitneyu(x=data_case, y=data_cont).pvalue
                res = pd.DataFrame([[clust, "-".join(dis_pair), np.mean(data_case), np.mean(data_cont),
                                     np.log(np.mean(data_case) / np.mean(data_cont)), pval]],
                                   columns=["Cluster", "Comparison", "mean_case", "mean_cont", "log(FC)", "P-val"])
                pair_res = pd.concat([pair_res, res])
            pair_res["qval"] = statsmodels.stats.multitest.fdrcorrection(pair_res["P-val"])[1]
            comp_res = pd.concat([comp_res, pair_res])

        comp_res.sort_values(by="qval", inplace=True)
        comp_res.to_csv(lineage + "_comparisons_stats.csv", index=False)
    else:
        comp_res = pd.read_csv(lineage + "_comparisons_stats.csv")
        df_freq = pd.read_csv(lineage + "_frequencies_per_pt.csv", index_col=0)







# This function create a pseudobulk count matrix and associated metadata for downstream task
# such as DESeq2 or WGCNA etc
#==========================================================================================================#
# data: an anndata, with .X contains the count to be used for aggregation
# sample_label: the label for the samples in data.obs
# cluster_label: the label for the clusters in data.obs on which the cells will be groups and counted as pseudo-bulk
# condition_label: the label for the condition in data.obs
# other_meta: other metadata columns that appear in obs to include in the metadata csv output file
# title: title of the cound and metadata files
# perClust: whether we are aggregating count for each Cluster or for all cells in the data object, default to True
# doublet_clust_to_rm: name of doublet clusters to remove before computing pseudo-bulk
# doublet_column: the name of the column in obs in which <doublet_clust_to_rm> is stored
#
# In the PB files that "sample" names will be constructed as so: <sample_label>:<condition_label>:<cluster_label>
# These names will be later deconstructed in the R code before DESeq2. It is therefore important that condition, sample and cluster names
# do not include the character ':' or other special characters such as space of backslah.
# Names of these features are slightly changed in this funciton to allow smooth run of DESeq2.
#=========================================================================================================#
def create_PB_files(data, sample_label, cluster_label, condition_label, other_meta, title, perClust=True, doublet_clust_to_rm=[], doublet_column=""):
    assert sample_label  in  data.obs.columns
    assert cluster_label  in data.obs.columns
    assert condition_label  in data.obs.columns

    # Making sure the raw layer contains only integers
    count_df_small = pd.DataFrame.sparse.from_spmatrix(data.raw.X[1:1000,1:1000])
    if not np.all(count_df_small.astype(int) == count_df_small):
        warnings.warn("Error: The data in data.raw is not integres! ")
    if not ("raw" in data.layers):
        data.layers["raw"] = data.raw.X.copy()


    # Removing doublets before computing pseudobulk
    if doublet_clust_to_rm!=[]:
        print("0) Removing doublets before computing pseudobulk")
        data = data[~data.obs[doublet_column].isin(doublet_clust_to_rm),:].copy()

    # Changing the cluster names to not include special characters - importnat for later code in R (DE_all_pairs.r)
    data.obs[cluster_label] = data.obs[cluster_label].str.replace("\\","_")
    data.obs[cluster_label] = data.obs[cluster_label].str.replace(" ", "_")
    data.obs[cluster_label] = data.obs[cluster_label].str.replace("/", "_")
    data.obs[cluster_label] = data.obs[cluster_label].str.replace("__", "_")
    data.obs[cluster_label] = data.obs[cluster_label].str.replace(":", "_")
    data.obs[condition_label] = data.obs[condition_label].str.replace("Term PET/GHTN", "TermPET")

    # Generating the pseudobulk
    if not perClust:
        print("1) Generating pseudobulk counts for all cells in the obejct")
        adpb = ADPBulk(data, [sample_label,condition_label], method="sum", layer="raw")
        pseudobulk_matrix = sc.AnnData(adpb.fit_transform())
        meta=adpb.get_meta()
        meta.index=meta[sample_label]
    else:
        print("1) Generating pseudobulk anndata per cluster")
        adpb = ADPBulk(data, [sample_label,condition_label,cluster_label], method="sum", layer="raw")
        pseudobulk_matrix = sc.AnnData(adpb.fit_transform()) # group indices and sums counts
        meta=adpb.get_meta()
        meta.index=meta[sample_label].astype(str)+":"+meta[condition_label].astype(str)+":"+meta[cluster_label].astype(str)

    del meta["SampleName"]
    meta["id"]=meta[sample_label]+"_"+meta[condition_label]
    print("2) Generating pseudobulk matrix")
    count_mtx = pd.DataFrame(pseudobulk_matrix.X.T, index=data.var_names,columns=meta.index)
    count_mtx.to_csv(title+"_pb_counts.csv")

    print("3) Generate metadata dataframe")
    if other_meta is not None:
        meta_others=data.obs[other_meta+[sample_label]].drop_duplicates()
        meta=pd.merge(meta,meta_others,how="left",on=sample_label)
    meta=meta.rename({sample_label:"sample",cluster_label:"cluster",condition_label:"condition"},axis=1)

    meta.to_csv(title+"_pb_meta.csv")

    # Saving the amount of cells per sample and per cluster we combined to a pseudo-bulk to allow later to filter clusters/ samples with too few cells
    df = data.obs[[sample_label, cluster_label]].value_counts().reset_index()
    df = df.rename(columns={sample_label:"sample", cluster_label:"cluster"})
    df.to_csv(title+"_pb_cells_per_samp_clust.csv")



# Wrapper barplot grid: Only DE bars
def multiDE_BarGrid(
    comparisons,
    pval_adj_thresh,
    logFC_thresh,
    pdf_name="DA_DE_grid.pdf",
    path_to_stats=".",
    DEGs_transparency_thresh=50, # A threshshold below which the bar will be partialy transparent
    onlyMinDEGsPop=False,       # include only populations with a minimum DEGs amount
    mimDEGs_thresh=0,
    share_y=False,
    xticks_fontsize=5,
    width=10,
    highlight_deg_diffs=False,      # turn on/off. Will outline populations that differ between 2 conditions in at lease <min_DEGs_diff> DEGs per direction
    min_DEGs_diff=50,               # required diff to add border
    border_lw=1.2,                   # border thickness
    global_order = []  # order of cell subsets on x axis
):
    lineage_order = ['Trophoblasts', 'EVT', 'Myeloid', 'T', 'NK','B_Plasma', 'fPV', 'fFB', 'dS', 'Endothelial']
    n_comp = len(comparisons)

    # ---- Read all DE tables first ----
    dfs = []
    for comp in comparisons:
        de_fname = f"{path_to_stats}/DE_{comp}_adjpval_{pval_adj_thresh}_logFC_{logFC_thresh}.csv"
        df = pd.read_csv(de_fname)
        df["Comparison"] = comp
        dfs.append((comp, df))

    # ---- Build a stable Cluster -> Lineage mapping across all comparisons ----
    all_df = pd.concat([d for _, d in dfs], ignore_index=True)

    if onlyMinDEGsPop:
        pass_by_comp = (
            all_df
            .assign(_pass=(all_df["DEG_up"] + all_df["DEG_down"] >= mimDEGs_thresh))
            .groupby(["Cluster", "Lineage", "Comparison"])["_pass"]
            .any())
        keep_pops = pass_by_comp.groupby(["Cluster", "Lineage"]).any()
        all_df = all_df.merge(keep_pops.rename("_keep").reset_index(), on=["Cluster", "Lineage"], how="left")
        all_df = all_df.loc[all_df["_keep"].fillna(False)].drop(columns=["_keep"]).copy()

    cluster_to_lin = (
        all_df.dropna(subset=["Cluster", "Lineage"])
              .drop_duplicates(subset=["Cluster"])
              .set_index("Cluster")["Lineage"]
              .to_dict()
    )

    def lin_rank(lin):
        return lineage_order.index(lin) if lin in lineage_order else 999

    if len(global_order) == 0:
        global_order = sorted(cluster_to_lin.keys(),key=lambda c: (lin_rank(cluster_to_lin.get(c, "Other")),natural_key(c)))

    fig, axes = plt.subplots(nrows=n_comp, ncols=1,figsize=(width, 2 * n_comp),squeeze=False, sharex=True, sharey=share_y)

    pad_template = pd.DataFrame({
        "Cluster": global_order,
        "Lineage": [cluster_to_lin.get(c, "Other") for c in global_order],
        "DEG_up": 0.0,
        "DEG_down": 0.0
    })

    # ---- Build all padded dfs first (so we can compare) ----
    padded = []
    for comp, df in dfs:
        df_small = df[["Cluster", "DEG_up", "DEG_down"]].copy()
        df_padded = pad_template.merge(df_small, on="Cluster", how="left", suffixes=("", "_orig"))
        for col in ["DEG_up", "DEG_down"]:
            df_padded[col] = df_padded[f"{col}_orig"].fillna(df_padded[col]).astype(float)
            df_padded.drop(columns=[f"{col}_orig"], inplace=True)
        padded.append((comp, df_padded))

    # ---- NEW (optional): compute which bars get borders ----
    # store boolean masks per comparison for up/down borders (length = n_clusters)
    border_masks = None
    if highlight_deg_diffs and len(padded) == 2:
        (comp0, d0), (comp1, d1) = padded

        up0 = d0["DEG_up"].to_numpy()
        up1 = d1["DEG_up"].to_numpy()
        dn0 = d0["DEG_down"].to_numpy()
        dn1 = d1["DEG_down"].to_numpy()

        # "same directionality" here means: compare UP with UP, DOWN with DOWN (not mixing)
        up_diff = up0 - up1
        dn_diff = dn0 - dn1

        up_border0 = up_diff >= min_DEGs_diff
        up_border1 = up_diff <= -min_DEGs_diff
        dn_border0 = dn_diff >= min_DEGs_diff
        dn_border1 = dn_diff <= -min_DEGs_diff

        border_masks = {
            comp0: {"up": up_border0, "down": dn_border0},
            comp1: {"up": up_border1, "down": dn_border1},
        }

    # ---- Plot ----
    for i, (comp, df_padded) in enumerate(padded):
        show_xticks = (i == n_comp - 1)

        draw_de_on_ax(
            df=df_padded,
            ax=axes[i, 0],
            x_col="Cluster",
            deg_up_col="DEG_up",
            deg_down_col="DEG_down",
            lineage_col="Lineage",
            order=global_order,
            lineage_order=lineage_order,
            palette=_DEF_PALETTE,
            DEGs_transparency_thresh=DEGs_transparency_thresh,
            title=f"{comp} – DEGs",
            show_xticklabels=show_xticks,
            show_separators=True,
            xticks_fontsize=xticks_fontsize
        )

        axes[i, 0].set_xlim(-0.5, len(global_order) - 0.5)

        # ---- NEW (optional): apply borders to the correct patches ----
        if border_masks is not None and comp in border_masks:
            n = len(global_order)
            patches = axes[i, 0].patches
            if len(patches) >= 2 * n:
                up_mask = border_masks[comp]["up"]
                dn_mask = border_masks[comp]["down"]

                # first n patches are up bars; next n are down bars
                for j in range(n):
                    if up_mask[j]:
                        patches[j].set_edgecolor("black")
                        patches[j].set_linewidth(border_lw)
                    if dn_mask[j]:
                        patches[n + j].set_edgecolor("black")
                        patches[n + j].set_linewidth(border_lw)

    # ---- One legend ----
    handles = [
        Patch(facecolor=_DEF_PALETTE[lin], edgecolor="black", label=str(lin), alpha=0.95)
        for lin in lineage_order if lin in _DEF_PALETTE
    ]
    fig.legend(handles=handles, title="Lineage", frameon=False,
               loc="center left", bbox_to_anchor=(1.02, 0.5))

    fig.tight_layout(rect=(0.0, 0.0, 0.88, 1.0))
    fig.savefig(pdf_name)
    plt.close(fig)



def multiDA_BarGrid(
    comparisons,
    pdf_name="DA_grid.pdf",
    path_to_stats=".",
    da_pval_thresh=0.05, # This threshold is used to have a transparent fill for populations with pval< <thresh>
    onlyMinDAPop=False,   # include only populations with a (p_val<0.05) AND (FC>minDA_absFC_thresh) in ANY comparison
    minDA_absFC_thresh=0.0, # threshold applied to abs(log(FC))
    share_y=False,
    xticks_fontsize=5,
    width=17,
    highlight_diff_between_comparisons=False,  #  turn border logic on/off - comparing prePE to TermPE
    border_lw_diffdir=1.2,                     #  border width when direction differs (both)
    border_lw_sigonly=1.2,                      #  border width thicker border when sig only in one
    global_order=[]   # order of cell subsets on x axis
):
    lineage_order = ['Trophoblasts', 'EVT', 'Myeloid', 'T', 'NK', 'B_Plasma', 'fPV', 'fFB', 'dS', 'Endothelial']
    n_comp = len(comparisons)
    stat_df = pd.read_csv(f"{path_to_stats}/all_DA_stats.csv")
    dfs = []
    for comp in comparisons:
        comp_for_DA = comp.replace("_vs_", "-")
        df = stat_df.loc[stat_df["Comparison"] == comp_for_DA].copy()
        df["Comparison"] = comp
        dfs.append((comp, df))

    all_df = pd.concat([d for _, d in dfs], ignore_index=True)
    if onlyMinDAPop:
        pass_by_comp = (
            all_df.assign(_pass=((all_df["P-val"] < 0.05) & (all_df["log(FC)"].abs() > minDA_absFC_thresh)))
                  .groupby(["Cluster", "Lineage", "Comparison"])["_pass"].any())
        keep_pops = pass_by_comp.groupby(["Cluster", "Lineage"]).any()
        all_df = all_df.merge(keep_pops.rename("_keep").reset_index(), on=["Cluster", "Lineage"], how="left")
        all_df = all_df.loc[all_df["_keep"].fillna(False)].drop(columns=["_keep"]).copy()

    cluster_to_lin = (all_df.dropna(subset=["Cluster", "Lineage"])
              .drop_duplicates(subset=["Cluster"])
              .set_index("Cluster")["Lineage"]
              .to_dict())

    def lin_rank(lin):
        return lineage_order.index(lin) if lin in lineage_order else 999

    if len(global_order)==0:
        global_order = sorted(cluster_to_lin.keys(), key=lambda c: (lin_rank(cluster_to_lin.get(c, "Other")), c))

    fig, axes = plt.subplots(
        nrows=n_comp, ncols=1,
        figsize=(width, 2 * n_comp),
        squeeze=False, sharex=True, sharey=share_y)

    pad_template = pd.DataFrame({
        "Cluster": global_order,
        "Lineage": [cluster_to_lin.get(c, "Other") for c in global_order],
        "log(FC)": 0.0,
        "P-val": 1.0,
    })

    # build padded dfs once
    padded = []
    for comp, df in dfs:
        df_small = df[["Cluster", "log(FC)", "P-val"]].copy()
        df_padded = pad_template.merge(df_small, on="Cluster", how="left", suffixes=("", "_orig"))
        for col in ["log(FC)", "P-val"]:
            df_padded[col] = df_padded[f"{col}_orig"].fillna(df_padded[col]).astype(float)
            df_padded.drop(columns=[f"{col}_orig"], inplace=True)
        padded.append((comp, df_padded))

    # ---- NEW (optional): compute border linewidths ----
    border_lw_by_comp = None
    if highlight_diff_between_comparisons:
        border_lw_by_comp = {comp: np.zeros(len(global_order), dtype=float) for comp, _ in padded}

        # only meaningful for the "compare two conditions" use-case
        if len(padded) == 2:
            (comp0, d0), (comp1, d1) = padded

            sig0 = (d0["P-val"].to_numpy() < da_pval_thresh) & (d0["log(FC)"].abs().to_numpy() > minDA_absFC_thresh)
            sig1 = (d1["P-val"].to_numpy() < da_pval_thresh) & (d1["log(FC)"].abs().to_numpy() > minDA_absFC_thresh)

            sgn0 = np.sign(d0["log(FC)"].to_numpy())
            sgn1 = np.sign(d1["log(FC)"].to_numpy())

            diff_dir   = sig0 & sig1 & (sgn0 != 0) & (sgn1 != 0) & (sgn0 != sgn1)
            sig_only0  = sig0 & (~sig1)
            sig_only1  = sig1 & (~sig0)

            lw0 = np.zeros(len(global_order), dtype=float)
            lw1 = np.zeros(len(global_order), dtype=float)

            lw0[diff_dir] = border_lw_diffdir
            lw1[diff_dir] = border_lw_diffdir
            lw0[sig_only0] = border_lw_sigonly
            lw1[sig_only1] = border_lw_sigonly

            border_lw_by_comp[comp0] = lw0
            border_lw_by_comp[comp1] = lw1

    # ---- Plot ----
    for i, (comp, df_padded) in enumerate(padded):
        show_xticks = (i == n_comp - 1)

        draw_da_on_ax(
            df=df_padded,
            ax=axes[i, 0],
            x_col="Cluster",
            fc_col="log(FC)",
            pval_col="P-val",
            lineage_col="Lineage",
            order=global_order,
            lineage_order=lineage_order,
            palette=_DEF_PALETTE,
            da_pval_thresh=da_pval_thresh,
            title=f"{comp} – DA",
            show_xticklabels=show_xticks,
            show_separators=True,
        )

        if show_xticks:
            axes[i, 0].tick_params(axis="x", rotation=90, labelsize=xticks_fontsize)

        axes[i, 0].set_xlim(-0.5, len(global_order) - 0.5)

        # ---- NEW (optional): apply borders ----
        if highlight_diff_between_comparisons and border_lw_by_comp is not None:
            lw_vec = border_lw_by_comp.get(comp, None)
            if lw_vec is not None and len(axes[i, 0].patches) == len(global_order):
                for patch, lw in zip(axes[i, 0].patches, lw_vec):
                    patch.set_edgecolor("black")
                    patch.set_linewidth(float(lw))

    # legend
    handles = [
        Patch(facecolor=_DEF_PALETTE[lin], edgecolor="black", label=str(lin), alpha=0.8)
        for lin in lineage_order if lin in _DEF_PALETTE
    ]
    fig.legend(handles=handles, title="Lineage", frameon=False, loc="center left", bbox_to_anchor=(1.02, 0.5))

    fig.tight_layout(rect=(0.0, 0.0, 0.88, 1.0))
    fig.savefig(pdf_name)
    plt.close(fig)


# Plots a double bar plot - on the top DA changes per subset, on the bottom number of DEGs per subset
# Only includes cell subsets that have at least
# comp - the comparison to plot
# pval_adj_thresh, logFC_thresh - the log(FC) and padj thresholds at which a DEG bar will be considered "significant" - bright colors
# pdf - whether to save the figure in a PDF file
# path - path to the files
# onlySig - leave in only significant values - either by DA or DE
# DE_file - if there is a specificDE stat file this function should use - mention it here, otherwise it will use the general file under <path> and based on the comp and thresholds
# mimDEGs_thresh - the minimum number of DEGs (summing up up and down) that a subset should have in order to be considered an important population and included in the barplot
def doubleBarPlot_DA_DE(
    comp,
    pval_adj_thresh,
    logFC_thresh,
    pdf=False,
    path=".",
    onlySig=False,
    figsize=(15, 6),
    DEGs_transparency_thresh=50,
    mimDEGs_thresh=50,
    DE_file="",
    da_pval_thresh=0.05,    # DA alpha threshold and the threshold to filter by if onlySig==TRue
    lineage_order=None,
    legend_rect=(0.0, 0.0, 0.84, 1.0),
    legend_anchor=(0.84, 0.5),
    global_order=[]
):
    if lineage_order is None:
        lineage_order = ['Trophoblasts', 'EVT', 'Myeloid', 'T', 'NK', 'B_Plasma', 'fPV', 'fFB', 'dS', 'Endothelial']
    palette = _DEF_PALETTE

    # --- load DE + DA ---
    if DE_file == "":
        de_path = f"{path}/DE_{comp}_adjpval_{pval_adj_thresh}_logFC_{logFC_thresh}.csv"
        DE_df = pd.read_csv(de_path)
    else:
        DE_df = pd.read_csv(DE_file)

    stat_df = pd.read_csv(f"{path}/all_DA_stats.csv")
    comp_for_DA = comp.replace("_vs_", "-")
    stat_curr = stat_df.loc[stat_df["Comparison"] == comp_for_DA, :].copy()
    df = stat_curr.merge(DE_df, on=["Cluster", "Lineage"], how="left")

    # fill missing DE counts if merge produced NaNs
    for c in ["DEG_up", "DEG_down"]:
        if c in df.columns:
            df[c] = df[c].fillna(0)

    # --- optional filtering ---
    if onlySig:
        df = df.loc[((df["DEG_up"] + df["DEG_down"]) >= mimDEGs_thresh) | (df["P-val"] <= da_pval_thresh), ].copy()

     # --- stacked layout: DA (top) + DE (bottom) ---
    plt.clf()
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(nrows=2, ncols=1, height_ratios=[1, 1], hspace=0.0)
    ax_top = fig.add_subplot(gs[0, 0])
    ax_bot = fig.add_subplot(gs[1, 0], sharex=ax_top)

    # top: DA
    used_order = draw_da_on_ax(
        df=df,
        ax=ax_top,
        x_col="Cluster",
        fc_col="log(FC)",
        pval_col="P-val",
        lineage_col="Lineage",
        order=None,  # let helper compute
        lineage_order=lineage_order,
        palette=palette,
        bar_width=0.85,
        fc_label="log(Fold-change)",
        title=comp,
        da_pval_thresh=da_pval_thresh,
        dim_alpha=0.25,
        bright_alpha=1.0,
        grid=True,
        show_xticklabels=False,
        show_separators=True,
    )

    # bottom: DE (same order)
    draw_de_on_ax(
        df=df,
        ax=ax_bot,
        x_col="Cluster",
        deg_up_col="DEG_up",
        deg_down_col="DEG_down",
        lineage_col="Lineage",
        order=used_order,
        lineage_order=lineage_order,
        palette=palette,
        bar_width=0.85,
        deg_label="# DE genes (up / down)",
        title=None,
        DEGs_transparency_thresh=DEGs_transparency_thresh,
        dim_alpha=0.25,
        bright_alpha=1.0,
        grid=True,
        show_xticklabels=True,
        show_separators=True,
        xticks_fontsize=8,
    )

    # cosmetics to make the two panels look “stacked”
    ax_top.spines["bottom"].set_visible(False)
    ax_bot.spines["top"].set_visible(False)
    ax_top.tick_params(axis="x", labelbottom=False)

    ax_top.set_xlim(-0.5, len(used_order) - 0.5)
    ax_bot.set_xlim(-0.5, len(used_order) - 0.5)
    ax_top.margins(x=0)
    ax_bot.margins(x=0)

    # --- one legend outside ---
    handles = [Patch(facecolor=palette[lin], edgecolor="black", label=str(lin), alpha=0.95)
                for lin in lineage_order
                if lin in palette]
    fig.legend(handles=handles, title="Lineage", frameon=False, loc="center left", bbox_to_anchor=legend_anchor)
    fig.tight_layout(rect=legend_rect)

    # --- save/show ---
    if pdf:
        out = f"{comp}_barplots.pdf"
        fig.savefig(out, bbox_inches="tight")
        plt.close(fig)
    else:
        plt.show()

    return fig, (ax_top, ax_bot)

###############################
# Grid DA/ DE plots
###############################

def natural_key(s):
    return [int(x) if x.isdigit() else x for x in re.split(r"(\d+)", str(s))]

# Helpers that draw DA onto given axes
def draw_da_on_ax(
    df,
    ax,
    x_col="Cluster",
    fc_col="log(FC)",
    pval_col="P-val",
    lineage_col="Lineage",
    order=None,
    lineage_order=None,
    palette=None,
    bar_width=0.85,
    fc_label="log(Fold-change)",
    title=None,
    da_pval_thresh=None,
    dim_alpha=0.25,
    bright_alpha=1,
    grid=True,
    show_xticklabels=True,
    show_separators=True,
):
    tit_for_file = title.replace(" ", "_")
    df.to_csv(tit_for_file+".csv")

    # sanity
    for col in [x_col, fc_col, pval_col, lineage_col]:
        if col not in df.columns:
            raise ValueError(f"Missing column: {col}")

    # --- Ordering clusters on x axis ---
    if lineage_order is None:
        lineage_order = list(df[lineage_col].dropna().unique())
    cat = pd.CategoricalDtype(categories=lineage_order, ordered=True)
    df[lineage_col] = df[lineage_col].astype(cat)

    if order is None:
        order = (df.sort_values(lineage_col).groupby(lineage_col, sort=False, observed=True)[x_col].apply(
            lambda x: sorted(x, key=natural_key)).explode().tolist())

    if palette is None:
        palette = _DEF_PALETTE.copy()
    elif isinstance(palette, str):
        cols = sns.color_palette(palette, len(lineage_order))
        palette = dict(zip(lineage_order, cols))

    plot_df = df.set_index(x_col).loc[order].reset_index()

    # DA highlight
    if da_pval_thresh is None:
        da_mask = np.zeros(len(plot_df), dtype=bool)
    else:
        da_mask = (plot_df[pval_col].to_numpy() <= da_pval_thresh)

    # main plot
    sns.barplot(data=plot_df,x=x_col, y=fc_col,ax=ax,dodge=False,edgecolor="black",linewidth=0.3,color="lightgray", width = 0.8)

    for idx, (patch, lin) in enumerate(zip(ax.patches, plot_df[lineage_col])):
        patch.set_facecolor(palette[lin])
        patch.set_alpha(bright_alpha if da_mask[idx] else dim_alpha)
        patch.set_width(bar_width)

    # separators between lineage blocks
    if show_separators:
        lineage_series = plot_df[lineage_col].tolist()
        boundaries = [
            i for i in range(1, len(lineage_series))
            if lineage_series[i] != lineage_series[i - 1]
        ]
        for b in boundaries:
            ax.axvline(b - 0.5, color="k", linewidth=0.6, alpha=0.35)

    if grid:
        ax.yaxis.grid(True, linestyle=":", alpha=0.5)

    ax.set_ylabel(fc_label)

    if show_xticklabels:
        ax.tick_params(axis="x", rotation=90, labelsize=5)
    else:
        ax.tick_params(axis="x", labelbottom=False)

    if title:
        ax.set_title(title, fontsize=10)

    return order  # we reuse this order for DE

# Helpers that draw DE onto given axes
def draw_de_on_ax(
    df,
    ax,
    x_col="Cluster",
    deg_up_col="DEG_up",
    deg_down_col="DEG_down",
    lineage_col="Lineage",
    order=None,
    lineage_order=None,
    palette=None,
    bar_width=0.85,
    deg_label="# DE genes (up / down)",
    title=None,
    DEGs_transparency_thresh=None,
    dim_alpha=0.25,
    bright_alpha=1,
    grid=True,
    show_xticklabels=True,
    show_separators=True,
    xticks_fontsize=7
):
    tit_for_file = title.replace(" ", "_")
    df.to_csv(tit_for_file+".csv")

    # sanity
    for col in [x_col, deg_up_col, deg_down_col, lineage_col]:
        if col not in df.columns:
            raise ValueError(f"Missing column: {col}")

    if lineage_order is None:
        lineage_order = list(df[lineage_col].dropna().unique())
    cat = pd.CategoricalDtype(categories=lineage_order, ordered=True)
    df[lineage_col] = df[lineage_col].astype(cat)

    if order is None:
        #order = (df.sort_values([lineage_col, x_col])[x_col].drop_duplicates().tolist())
        order = (df.sort_values(lineage_col).groupby(lineage_col, sort=False, observed=True)[x_col].apply(
                lambda x: sorted(x, key=natural_key)).explode().tolist())

    if palette is None:
        palette = _DEF_PALETTE.copy()
    elif isinstance(palette, str):
        cols = sns.color_palette(palette, len(lineage_order))
        palette = dict(zip(lineage_order, cols))

    plot_df = df.set_index(x_col).loc[order].reset_index()

    up_vals   = plot_df[deg_up_col].fillna(0).to_numpy().astype(float)
    down_vals = plot_df[deg_down_col].fillna(0).to_numpy().astype(float)

    if DEGs_transparency_thresh is None:
        deg_mask_up = np.zeros(len(plot_df), dtype=bool)
        deg_mask_down = np.zeros(len(plot_df), dtype=bool)
    else:
        deg_mask_up   = (up_vals + down_vals) >= DEGs_transparency_thresh
        deg_mask_down = (up_vals + down_vals) >= DEGs_transparency_thresh

    xpos = np.arange(len(plot_df))
    colors = [palette[lin] for lin in plot_df[lineage_col]]

    # up bars
    bars_up = ax.bar(xpos, up_vals, width=bar_width, edgecolor="black", linewidth=0.3, color=colors)
    for b, m in zip(bars_up, deg_mask_up):
        b.set_alpha(bright_alpha if m else dim_alpha)

    # down bars
    bars_down = ax.bar(xpos, -down_vals, width=bar_width, edgecolor="black", linewidth=0.3, color=colors)
    for b, m in zip(bars_down, deg_mask_down):
        b.set_alpha(bright_alpha if m else dim_alpha)

    # separators between lineage blocks
    if show_separators:
        lineage_series = plot_df[lineage_col].tolist()
        boundaries = [
            i for i in range(1, len(lineage_series))
            if lineage_series[i] != lineage_series[i - 1]
        ]
        for b in boundaries:
            ax.axvline(b - 0.5, color="k", linewidth=0.6, alpha=0.35)

    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_ylabel(deg_label)

    if show_xticklabels:
        ax.set_xticks(xpos)
        ax.set_xticklabels(plot_df[x_col], rotation=90, fontsize=xticks_fontsize)
    else:
        ax.set_xticks(xpos)
        ax.set_xticklabels([""] * len(xpos))

    if grid:
        ax.yaxis.grid(True, linestyle=":", alpha=0.3)

    if title:
        ax.set_title(title, fontsize=10)

    ax.xaxis.grid(False)





# This function plots 3 panels - on the left a UMAP colored by cluster_col, on the right the same UMAP colored by fetal/ maternal cell annotations and in the middle a barplot quantifying the frequecnies of fetal/ maternal cells per cluster
def plot_UMAP_MatFet_trio(
    lin_obj_path,  # either the object itself or the path to the object
    title,
    cluster_col="sub_lin2",
    cluster_order=None,
    xlabel="% cells",
    xlabel_fontsize=12,
    legend_fontsize=10,
    legend_title_fontsize=10,
    figsize=(10, 2.6),
    cluster_palette=None, # colors for cluster UMAP
    add_bar_count= False, # Whether to add a fourth panel
    count_ticks=None, # list of ticks to have for x axis in the count barplot (the third)
    dotsize=None      #
):
    if (type(lin_obj_path)==str):
        lin_data = mu.read(lin_obj_path)
    else:
        lin_data = lin_obj_path
    ad = lin_data["GEX"]
    obs = ad.obs

    # fallback if user passed sub_lin2 but only sub_lin exists
    if cluster_col not in obs.columns:
        if cluster_col == "sub_lin2" and "sub_lin" in obs.columns:
            cluster_col = "sub_lin"
        else:
            raise ValueError(f"{cluster_col} not in adata.obs")

    # ---- choose ONE canonical cluster order ----
    if cluster_order is None:
        if pd.api.types.is_categorical_dtype(obs[cluster_col]):
            cluster_order = list(obs[cluster_col].cat.categories)
        else:
            cluster_order = list(pd.unique(obs[cluster_col].astype(str)))

    obs[cluster_col] = pd.Categorical(obs[cluster_col].astype(str),
                                     categories=cluster_order, ordered=True)

    # Fet/Mat colors
    obs["FetMat_Annotation"] = obs["FetMat_Annotation"].astype(str).replace({"nan": "Unassigned"})
    FetMat_col = {"Fetal": "black", "Maternal": "#B58B47", "Unassigned": "lightgrey"}
    FetMat_col = {k: FetMat_col[k] for k in pd.unique(obs["FetMat_Annotation"]) if k in FetMat_col}

    # ---- layout ----
    fig = plt.figure(figsize=figsize)
    if (add_bar_count):
        gs = gridspec.GridSpec(1, 4, width_ratios=[1, 0.7, 0.7, 1], wspace=0.65)
        ax0 = fig.add_subplot(gs[0, 0])
        ax1 = fig.add_subplot(gs[0, 1])
        ax3 = fig.add_subplot(gs[0, 2])  # counts
        ax2 = fig.add_subplot(gs[0, 3])  # right UMAP
    else:
        gs = gridspec.GridSpec(1, 3, width_ratios=[1, 0.7, 1], wspace=0.65)
        ax0 = fig.add_subplot(gs[0, 0])
        ax1 = fig.add_subplot(gs[0, 1])
        ax2 = fig.add_subplot(gs[0, 2])

    # Left UMAP (clusters)
    sc.pl.umap(ad, color=cluster_col, title="", ax=ax0, palette=(cluster_palette if isinstance(cluster_palette, dict) else None))

    # --- dots strip from cluster legend colors ---
    handles, labels = ax0.get_legend_handles_labels()
    if ax0.get_legend() is not None:
        ax0.get_legend().remove()

    if (cluster_palette==None):
        dot_colors = []
        for h in handles:
            c = None
            if hasattr(h, "get_facecolor"):
                fc = h.get_facecolor()
                if fc is not None and len(fc) > 0:
                    c = fc[0]
            if c is None and hasattr(h, "get_color"):
                c = h.get_color()
            dot_colors.append(c)
    else:
        # if dict: map cluster_order -> colors; if list: assume aligned to cluster_order
        if isinstance(cluster_palette, dict):
            dot_colors = [cluster_palette.get(c, "lightgrey") for c in cluster_order]
        else:
            dot_colors = list(cluster_palette)

    ax0dots = ax0.inset_axes([1.01, 0.0, 0.03, 1.0], transform=ax0.transAxes)
    n = len(dot_colors)
    if n > 0:
        ys = np.linspace(0.95, 0.05, n)
        if (cluster_palette==None):
            ax0dots.scatter(np.full(n, 0.5), ys, s=25, c=dot_colors, marker="o")
        else:
            ax0dots.scatter(np.full(n, 0.5), ys, s=25, c=dot_colors, marker="o")
    ax0dots.set_xlim(0, 1)
    ax0dots.set_ylim(0, 1)
    ax0dots.axis("off")

    # --- Fet/Mat legend (text) aligned left + larger font ---
    # Build barplot
    df = pd.crosstab(obs[cluster_col], obs["FetMat_Annotation"]).reindex(cluster_order).dropna(how="all")
    df = df.div(df.sum(axis=1), axis=0) * 100.0
    df.plot(kind="barh", stacked=True, grid=False, color=FetMat_col, ax=ax1, width=0.8)
    ax1.invert_yaxis()
    ax1.set_xlabel(xlabel, fontsize=xlabel_fontsize)
    ax1.set_ylabel("")

    leg1 = ax1.legend(title="",frameon=False,loc="upper left",bbox_to_anchor=(1.02, 1.0),fontsize=legend_fontsize,handletextpad=0.6,borderaxespad=0.0,labelspacing=0.4,alignment="left")   # <- left-align legend entries
    if leg1.get_title() is not None:
        leg1.get_title().set_fontsize(legend_title_fontsize)

    # Right UMAP (Fet/Mat) — remove legend to avoid duplication
    if (dotsize!=None):
        sc.pl.umap(ad, color="FetMat_Annotation", title="", palette=FetMat_col, ax=ax2, size=dotsize)
    else:
        sc.pl.umap(ad, color="FetMat_Annotation", title="", palette=FetMat_col, ax=ax2)
    if ax2.get_legend() is not None:
        ax2.get_legend().remove()

    # ---- OPTIONAL: count barplot (same cluster order as legend/barh) ----
    if add_bar_count:
        counts = obs[cluster_col].value_counts().reindex(cluster_order).fillna(0).astype(int)
        df_counts = counts.reset_index()
        df_counts.columns = ["Cluster", "Count"]

        sns.barplot(data=df_counts, y="Cluster", x="Count", hue="Cluster", palette=cluster_palette, dodge=False, ax=ax3, edgecolor="black", linewidth=0.3, width = 0.8)
        if ax3.get_legend() is not None:
            ax3.get_legend().remove()
       # ax3.invert_yaxis()  # match ax1/cluster order
        ax3.set_xlabel("Cells", fontsize=xlabel_fontsize)
        ax3.set_ylabel("")
        if (count_ticks!=None):
            ax3.set_xticks(count_ticks)
            ax3.set_xticklabels([f"{int(round(t))}" for t in count_ticks])
        ax3.tick_params(axis="y", labelleft=False)  # cluster labels already shown on ax1
        ax3.grid(False)

    # ---- thin borders + ticks everywhere ----
    axes_all = (ax0, ax1, ax2) if not add_bar_count else (ax0, ax1, ax3, ax2)
    for ax in axes_all:
        for s in ["left", "right", "top", "bottom"]:
            ax.spines[s].set_linewidth(0.5)
            ax.tick_params(width=0.5, length=3, labelsize=8)

    fig.savefig(f"{title}_umap_trio.pdf", bbox_inches="tight", dpi=600)



def composition_barplot_h(adata, cluster_col, yattr, color_dict=None, ax=None):
    obs = adata.obs
    # desired order = categorical order if present, else sorted unique
    if pd.api.types.is_categorical_dtype(obs[cluster_col]):
        order = list(obs[cluster_col].cat.categories)
    else:
        order = sorted(obs[cluster_col].astype(str).unique())

    df = pd.crosstab(obs[cluster_col].astype(str), obs[yattr].astype(str))
    df = df.reindex(order).dropna(how="all")               # enforce order
    df = df.div(df.sum(axis=1), axis=0) * 100.0

    ax = df.plot(kind="barh", stacked=True, grid=False, color=color_dict, ax=ax)
    ax.invert_yaxis()  # optional: keep top-to-bottom same as factor order
    return ax





# lin_freq_file created by barplot_lin_per_pt
# withSampAnnotation - boolean variable that indicates if the sample annotation should be included in the plot
def boxplot_DA_lin(lin_freq_file,
                   pairs,
                   lin,
                   pdf=None,
                   condition_col = "Condition",
                   figsize=(4,4),
                   ):
    # samp_cond - is a csv file that includes the association of every sample to its condition. A sample can be associated with more than one condition
    condition_colors = dict(zip(["Cont_labored", "Cont_unlabored", 'Cont_PTB', "DM1", "FGR", "PrePE", "TermPE", "PTB"],
     # ["#1e4923", "#489142", "#99cf91", "#416eb6", "#533678", "#a06da5", "#a7a3d1", "#98989a"]))
        ["#489142", "#489142", "#489142", "#BF0A3A", "#BF0A3A", "#BF0A3A", "#BF0A3A", "#BF0A3A"]))
    samp_cond = pd.read_csv("/projects/Pregnancy_placenta/Results/13_DE/sample_condition_list.csv")
    lin_freq = pd.read_csv(lin_freq_file)

    sub_freq = lin_freq[["Sample",lin]]
    sub_freq = samp_cond.merge(sub_freq, on="Sample", how="left")
    sub_freq["Condition"] = pd.Categorical(sub_freq[condition_col], ordered=True, categories=["Cont_labored","Cont_unlabored","Cont_PTB","PTB","DM1","FGR","PrePE","TermPE"])

    # All conditions
    plt.clf()
    fig, ax = plt.subplots(1, 1, figsize=figsize)
    ax = sns.boxplot(data = sub_freq, x="Condition", y=lin, hue="Condition", width=0.75, fliersize=0, ax=ax, palette=condition_colors)
    sns.stripplot(data=sub_freq, x="Condition", y=lin, dodge=True, size=4, linewidth=0.25, ax=ax, color="black")
    fig.suptitle("Abundnace of " + lin + " (Mann-Whitney p-value)")
    plt.ylabel("Abundnace of "+lin)
    plt.tick_params(axis='x', labelrotation=90, size=6)

    annotator = Annotator(data=sub_freq, x="Condition", y=lin, ax=ax, pairs=pairs)
    annotator.configure(test='Mann-Whitney',
                        hide_non_significant=True,
                        line_width=1,
                        fontsize=8,
                        line_offset=0.01,  # <- distance from data to first bracket (smaller)
                        line_offset_to_group=0.005,  # <- distance between stacked brackets (smaller)
                        line_height=0  # <- bracket "height" (smaller)
                        ).apply_and_annotate()

    plt.tight_layout()
    if (pdf == None):
        plt.show()
    else:
        pdf.savefig(fig)






def dot_plot_with_umap_stack(
    rna_anndata,
    prot_anndata=None,
    groupby="Final_cluster",
    strip_categories=None,
    horizontal_strip=True,
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
    text_scale=0.7,
    dot_size_exponent=1.3,
    stroke_width=0.3
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

    # --- global styling scalers ---
    ytick_fontsize = int(ytick_fontsize * text_scale)
    xtick_fontsize = int(xtick_fontsize * text_scale)
    category_size = int(category_size * text_scale)
    legend_fontsize = int(legend_fontsize * text_scale)

    legend_fontoutline = max(1, int(legend_fontoutline * text_scale))
    panel_width = stroke_width

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
                raise ValueError(f"These protein markers are missing from prot_anndata.var_names: {list(missing_prot)}")

            prot_pb = get_mean_pct(
                prot_anndata[:, protein.tolist()].copy(),
                groupby,
                standard_scale=standard_scale
            )
            prot_pb["feature"] = "protein"
            prot_pb["block"] = prot_pb["Gene"].astype(str).replace(
                dict(zip(strip_categories.loc[strip_categories["types"] == "Protein", "markers"],
                         strip_categories.loc[strip_categories["types"] == "Protein", "blocks"])))
            pb_list.append(prot_pb)

    if rna_anndata is not None:
        genes = strip_categories.loc[strip_categories["types"] == "RNA", "markers"]
        if len(genes) > 0:
            missing_rna = genes[~genes.isin(rna_anndata.var_names)]
            if len(missing_rna) > 0:
                raise ValueError(f"These RNA markers are missing from rna_anndata.var_names: {list(missing_rna)}")

            rna_pb = get_mean_pct(rna_anndata[:, genes.tolist()].copy(),groupby,standard_scale=standard_scale)
            rna_pb["feature"] = "RNA"
            rna_pb["block"] = rna_pb["Gene"].astype(str).replace(
                dict(zip(
                    strip_categories.loc[strip_categories["types"] == "RNA", "markers"],
                    strip_categories.loc[strip_categories["types"] == "RNA", "blocks"])))
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
        size=umap_dotsize,
        return_fig=False,
        show=False,
        ax=umap_ax,
        legend_loc=umap_legend_loc,
        legend_fontsize=legend_fontsize,
        legend_fontoutline=legend_fontoutline
    )
    umap_ax.set_aspect(aspect=umap_aspect, adjustable="box")

    legend_colors = umap_anndata.uns[color_key_for_legend]
    if len(legend_colors) < len(full_categories):
        legend_colors = (legend_colors * (len(full_categories) // len(legend_colors) + 1))[:len(full_categories)]

    cat_to_color = dict(zip(full_categories, legend_colors))
    plot_colors = [cat_to_color[c] for c in categories_to_plot[::-1]]

    sns.scatterplot(
        x=np.repeat(1, len(categories_to_plot)),
        y=pd.Series(np.arange(0, len(categories_to_plot)) + margin),
        c=plot_colors,
        s=legend_dotsize,
        linewidth=stroke_width,
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

    legend_ax.scatter(legend["x"], legend["y"], s=legend["size"] ** dot_size_exponent, c="grey", edgecolors="black",linewidth=stroke_width)
    legend_ax.tick_params(which='major', length=7)
    legend_ax.set_xticks(np.arange(1, 6, 1) ** legend_scale)
    legend_ax.set_xlim([-1 * margin, 4])
    legend_ax.set_yticklabels([])
    legend_ax.set_xticklabels(legend["size"], size=int(8 * text_scale))
    legend_ax.spines['top'].set_visible(False)
    legend_ax.spines['right'].set_visible(False)
    legend_ax.spines['left'].set_visible(False)
    legend_ax.spines['bottom'].set_visible(False)
    legend_ax.set_title("Fraction in group\nexpress genes (%)", loc="center", size=int(10 * text_scale))

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
            im = ax[k + 1].scatter(block_pb["Gene"], block_pb[groupby], c=block_pb["mean"], s=block_pb["pct"].astype(float) ** dot_size_exponent, cmap="Reds", vmin=rna_vmin, vmax=rna_vmax, edgecolors="black", linewidth=stroke_width)
            if first_rna_im is None:
                first_rna_im = im
        else:
            im = ax[k + 1].scatter(block_pb["Gene"], block_pb[groupby], c=block_pb["mean"], s=block_pb["pct"].astype(float) ** dot_size_exponent, cmap="Blues",
                vmin=prot_vmin, vmax=prot_vmax, edgecolors="black", linewidth=stroke_width)
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

        # --- attach colorbars ---
        rna_cbar = fig.colorbar(first_rna_im, cax=cbar_ax2, orientation=cbar_orientation) if (
                    cbar_ax2 is not None and first_rna_im is not None) else None
        prot_cbar = fig.colorbar(first_prot_im, cax=cbar_ax, orientation=cbar_orientation) if (
                    cbar_ax is not None and first_prot_im is not None) else None

        tick_fs = int(10 * text_scale)
        label_fs = int(10 * text_scale)

        # --- format (works for RNA-only, CITE-only, or both) ---
        for cb, ax_cb, lab in [(prot_cbar, cbar_ax, "CITE"), (rna_cbar, cbar_ax2, "RNA")]:
            if cb is None:
                continue
            cb.ax.tick_params(labelsize=tick_fs, width=stroke_width, length=2)
            if cbar_orientation == "vertical":
                ax_cb.set_xlabel(lab, fontsize=label_fs, labelpad=4)
            else:
                cb.set_label(lab, fontsize=label_fs, labelpad=2)
                cb.ax.xaxis.set_label_position("bottom")
            cb.update_ticks()

        # title: put it once (on CITE if present, else on RNA)
        if cbar_orientation == "vertical":
            if cbar_ax is not None and prot_cbar is not None:
                cbar_ax.set_title(title, pad=12, fontsize=label_fs)
            elif cbar_ax2 is not None and rna_cbar is not None:
                cbar_ax2.set_title(title, pad=12, fontsize=label_fs)
        else:
            if cbar_ax is not None and prot_cbar is not None:
                cbar_ax.text(0.5, 1.65, title, transform=cbar_ax.transAxes, ha="center", va="bottom",
                             fontsize=int(8 * text_scale))
            elif cbar_ax2 is not None and rna_cbar is not None:
                cbar_ax2.text(0.5, 1.65, title, transform=cbar_ax2.transAxes, ha="center", va="bottom",
                              fontsize=int(8 * text_scale))
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




def feature_plots_pdf(adata, gene_list, title, use_raw=False, hex=False, gridsize=100, umap="umap", pdf=True, figsize=(10, 7), dot_size=5, axes_labelsize=25, axes_titlesize= 35, legend_fontsize= 20, nrows=4, ncols=4, out_dir=".", ticks = [0, 1, 2, 3, 4, 5,6], tick_labels = ["0", "1", "2", "3", "4", "5", "6+"], vmax=6, cmap='YlOrRd'):
    sc.set_figure_params(vector_friendly=False, dpi_save=400)
    genes = [g for g in gene_list if g in adata.var_names]
    print(str(len(genes)) + " out of " + str(len(gene_list)) + " were measured")
    gene_per_page = ncols*nrows
    print("gene per page: "+str(gene_per_page))
    page_no = math.floor(len(genes)/gene_per_page)
    if(len(genes)%gene_per_page != 0):
        page_no+=1

    if (umap != "umap"):
        adata = adata.copy()
        adata.obsm["X_umap"] = adata.obsm["X_" + umap]  # sc.pl.scatter doesn't have all the feature os sc.pl.umap to plot correctly

    def plot(dot_size):
        fig_list = []
        count = 0
        for page in range(page_no):
            print("page: " + str(page))
            if (hex):
                x = adata.obsm["X_umap"][:, 0]
                y = adata.obsm["X_umap"][:, 1]
                fig, axs = plt.subplots(ncols=ncols, nrows=nrows, figsize=figsize)
                for i in range(nrows):
                    for j in range(ncols):
                        if (count > (len(genes) - 1)):
                            break
                        gene = genes[count]
                        count += 1
                        ax = axs[i][j]
                        hb = ax.hexbin(x, y, C=adata[:, gene].to_df()[gene], cmap=cmap, gridsize=gridsize, linewidths=0, rasterized=False, vmin=0, vmax=vmax)
                        ax.set_title(gene, fontsize=8)
                        ax.set_xticks([])
                        ax.set_yticks([])
                        ax.set_ylabel("")
                        ax.set_xlabel("")
                        #ax.set_rasterization_zorder(10)
                        ax.grid(False)
                        ax.tick_params(labelsize=2)
                        cb= fig.colorbar(hb, ax=ax)
                        #cb.set_label(f"{gene} mean expression (capped)")
                        cb.set_ticks(ticks)
                        cb.set_ticklabels(tick_labels)
                        cb.ax.tick_params(labelsize=6)
            else:
                parameters = {'axes.labelsize': 0, 'axes.titlesize': axes_titlesize, 'axes.grid': False, 'legend.fontsize': "medium", 'figure.figsize':figsize}
                plt.rcParams.update(parameters)
                fig = sc.pl.umap(adata, color=genes[gene_per_page * page:gene_per_page * (page + 1)], size=dot_size, return_fig=True, show=False,cmap=cmap, use_raw=use_raw, ncols= ncols)
                fig.set_size_inches(figsize[0], figsize[1])

                # rasterize only plotted collections
                def is_umap_axis(ax, n_cells):   # Keep only axes that contain the actual UMAP scatter points. Exclude colorbar axes.
                    for coll in ax.collections:
                        if hasattr(coll, "get_offsets"):
                            offsets = coll.get_offsets()
                            if hasattr(offsets, "shape") and offsets.shape[0] == n_cells:
                                return True
                    return False

                plot_axes = [ax for ax in fig.axes if is_umap_axis(ax, adata.n_obs)]
                for ax in plot_axes:
                    for coll in ax.collections:
                        coll.set_rasterized(True)

                # remove axes titles so they don't get flattened with the panel
                for ax in plot_axes:
                    ax.set_title("")
                # force layout before placing figure-level text
                fig.canvas.draw()
                # add titles back as figure text (editable in Illustrator)
                page_genes = genes[gene_per_page * page:gene_per_page * (page + 1)]
                for ax, gene in zip(plot_axes, page_genes):
                    pos = ax.get_position()
                    fig.text(pos.x0 + pos.width / 2, pos.y1 + 0.005, gene.replace("G:", ""), ha="center", va="bottom", fontsize=8)

            for ax in fig.axes:
                ax.set_title(ax.get_title().replace("G:", ""))
                ax.set_title(ax.get_title().replace("P:", ""))

            fig_list.append(fig)

        return(fig_list)

    if (out_dir != ".") & (~os.path.exists(out_dir)):
        os.mkdir(out_dir)
    if (pdf):
        if (dot_size==""):
            dot_size=25
        with PdfPages(out_dir+"/"+title.replace(" ","_")+".pdf") as pdf:
            fig_list = plot(dot_size=dot_size)
            for fig in fig_list:
                #fig.suptitle(title, fontsize=axes_titlesize, y=0.995)
                fig.subplots_adjust(top=0.90)
                pdf.savefig(fig, bbox_inches="tight", dpi=300)
    else:
        if (dot_size==""):
            dot_size=20
        fig_list = plot(dot_size=dot_size)
        for fig in fig_list:
            plt.suptitle(title, fontsize=nrows*15, y=0.995)
            plt.subplots_adjust(hspace=0.8)
            fig.show()






# Computes mean expression of a gene by sample
def mean_exp_by_sample(adata, gene, subset_col="short_sublin2", subset=None, sample_col="Sample"):
    gene = "G:" + gene
    if gene not in adata.var_names:
        return None
    x = adata[adata.obs[subset_col] == subset, gene]
    if x.n_obs == 0:
        return None
    vals = x.X.toarray().ravel() if hasattr(x.X, "toarray") else np.ravel(x.X)
    out = pd.DataFrame({sample_col: x.obs[sample_col].values, "exp": vals})
    return out.groupby(sample_col).agg(exp=("exp", "mean"), n_cells=("exp", "size")).reset_index()



# Plot correlation scatter plots between myeloid cells and T cells
# including only samples with at least 30 cells per subset
# This functions assumes there are 2 T cell subsets and 1 myeloid subset
def scatter_co_expression(t_data, mMac_data, subset_pairs, gene_list, pdf_title):
    # sample-condition table from all relevant cells
    samp_cond = pd.concat([t_data.obs[["Sample", "Condition"]], mMac_data.obs[["Sample", "Condition"]]]).drop_duplicates()

    # Colelcting L-R gene expression and computing correlation
    exp_df = pd.DataFrame()
    stat_df = pd.DataFrame()
    for _, row in gene_list.iterrows():
        T_gene = row["T_gene"]                # T_gene = "CD80"
        myeloid_gene = row["myeloid_gene"]    # myeloid_gene = "CTLA4"
        print(T_gene+" - "+myeloid_gene)

        # first pass: compute all correlations/p-values
        for (t_sub, my_sub) in subset_pairs:         # t_sub= 'mT_4_CD8_eff' ; my_sub = 'mMac_5'
            T_samp_mean = mean_exp_by_sample(t_data, T_gene, subset=t_sub)
            M_samp_mean = mean_exp_by_sample(mMac_data, myeloid_gene, subset=my_sub)
            if (T_samp_mean is None) or (M_samp_mean is None):
                results.append({"t_sub": t_sub, "my_sub": my_sub, "df": None, "corr": np.nan, "pval": np.nan})
                continue
            df = pd.merge(T_samp_mean, M_samp_mean, on="Sample", suffixes=("_T", "_Myeloid"))
            df = df.rename(columns={"exp_T": "T_exp", "n_cells_T": "T_ncells", "exp_Myeloid": "Myeloid_exp", "n_cells_Myeloid": "Myeloid_ncells"})

            # filter for samples with at least 30 cells in both subsets
            df = df[(df["T_ncells"] >= 30) & (df["Myeloid_ncells"] >= 30)]
            df = pd.merge(df, samp_cond, on="Sample", how="left").drop_duplicates()

            if len(df) < 10:
                stat_df.append({"t_sub": t_sub, "my_sub": my_sub, "df": None, "corr": np.nan, "pval": np.nan})
                continue

            corr, pval = spearmanr(df["Myeloid_exp"], df["T_exp"])
            curr_stat = pd.DataFrame({"T_sub": t_sub, "T_gene":T_gene, "Myeloid_sub": my_sub, "Myeloid_gene":myeloid_gene, "No_samp_wOver30cells": df.shape[0], "corr": corr, "pval": pval}, index=[0])
            stat_df = pd.concat([stat_df, curr_stat])

            # collecting info for exp df
            df["t_sub"] = t_sub
            df["my_sub"] = my_sub
            df["T_gene"] = T_gene
            df["Myeloid_gene"] = myeloid_gene
            exp_df = pd.concat([exp_df, df])

    # FDR correction across all ligand-receptor pairs of one cell tpye
    stat_df.reset_index(drop=True, inplace=True)
    stat_df["adj_pval"] = np.nan
    for cell in stat_df["T_sub"].unique():
        mask = stat_df["T_sub"] == cell
        stat_df.loc[mask, "adj_pval"] = multipletests(stat_df.loc[mask, "pval"],method="fdr_bh")[1]

    # Saving files
    exp_df.to_csv("Fig6G_expression_coinhibitory.csv", index=False)
    stat_df.to_csv("Fig6G_stats_coinhibitory.csv", index=False)

    # Plotting
    t_pops = list(set(stat_df["T_sub"]))
    myeloid_pop = list(set(stat_df["Myeloid_sub"]))[0]
    sns.set_theme(rc={"figure.dpi": 300})
    condition_colors = dict(zip(["Control_lbr", "Control_unlbr", 'PTB (GA Controls)', "DM1", "FGR", "PrePET", "TermPET", "PTB",'PrePET AND FGR'],
                                ["#274001", "#03A63C", "#B5FF63", "#7365A6", "#FF4F00", "#F05C8F", "#D9A711", "#BF0A3A", "#FF4F00"]))
    sns.set_style("ticks")
    with PdfPages(pdf_title) as pdf:
        for _, row in gene_list.iterrows():
            fig, axes = plt.subplots(1, 2, figsize=(5, 3))
            T_gene = row["T_gene"]
            myeloid_gene = row["myeloid_gene"]

            for i, (ax, t_sub) in enumerate(zip(axes, t_pops)):
                curr_exp = exp_df.loc[
                    (exp_df["t_sub"] == t_sub) &
                    (exp_df["T_gene"] == T_gene) &
                    (exp_df["Myeloid_gene"] == myeloid_gene), :
                ]
                curr_stat = stat_df.loc[
                    (stat_df["T_sub"] == t_sub) &
                    (stat_df["T_gene"] == T_gene) &
                    (stat_df["Myeloid_gene"] == myeloid_gene), :
                ]

                if curr_exp.empty:
                    ax.set_axis_off()
                    ax.set_title(f"{my_sub} vs {t_sub}\nmissing gene/data", fontsize=10, color="black")
                    continue

                sns.scatterplot(
                    data=curr_exp,
                    x="T_exp",
                    y="Myeloid_exp",
                    hue="Condition",
                    palette=condition_colors,
                    ax=ax,
                    legend=(i == 0)  # only first subplot creates legend
                )

                adj_pval = curr_stat["adj_pval"].iloc[0]
                pval = curr_stat["pval"].iloc[0]
                corr = curr_stat["corr"].iloc[0]

                title_color = "red" if adj_pval < 0.05 else "black"

                ax.set_title(
                    f"{myeloid_pop} vs {t_sub}\nR={corr:.2f}, p={pval:.3g}, FDR={adj_pval:.3g}",
                    fontsize=10,
                    color=title_color
                )
                ax.set_xlabel(f"{T_gene} in {t_sub}", fontsize=9)
                ax.set_ylabel(f"{myeloid_gene} in {myeloid_pop}", fontsize=9)

            handles, labels = axes[0].get_legend_handles_labels()
            leg = axes[0].get_legend()
            if leg is not None:
                leg.remove()

            if handles:
                fig.legend(handles, labels, loc="center left", bbox_to_anchor=(1.02, 0.5))

            fig.suptitle(f"{myeloid_gene} (myeloid) vs {T_gene} (T)", fontsize=14)
            plt.tight_layout(rect=[0, 0, 0.88, 0.96])
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)




