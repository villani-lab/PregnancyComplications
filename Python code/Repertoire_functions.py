import scanpy as sc
import anndata as an
import seaborn as sns
import numpy as np
import pandas as pd
from fontTools.ttLib.tables.F__e_a_t import Feature
from matplotlib.colors import TwoSlopeNorm
from collections import Counter
import matplotlib.pyplot as plt
from statannotations.Annotator import Annotator
from matplotlib.backends.backend_pdf import PdfPages
import re
from statsmodels.stats.multitest import multipletests
import math

code_path="/projects/Pregnancy_placenta/Code/Publication_code"


# Compute clone proportion per sample
# A clone is considered expanded if it has at least  min_cell_expanded cells and has a proportion of above min_prop_expanded per sample
def calc_clone_percent(adata, sample_col, min_cell_expanded = 5, min_prop_expanded = 0.01) :
    clone_perc_dict = {}
    clone_exp_dict = {}
    for cat in set(adata.obs[sample_col]):
        subj_data = adata[adata.obs[sample_col] == cat]
        clone_df = pd.DataFrame(Counter(subj_data.obs["clone"]), index = ["n"]).T
        # Taking into account only cells with legit clones
        if (any(clone_df.index == "")):
            clone_df = clone_df.drop(index = "")
        clone_df["perc"] = (clone_df["n"]*100) / np.sum(clone_df["n"])
        clone_dict = dict(zip(clone_df.index.values, clone_df["perc"]))
        clone_perc_dict[cat] = clone_dict
        # Lets define expanded clones
        clone_df["expanded"] = [True if (n > min_cell_expanded) & (perc > min_prop_expanded) else False for n, perc in zip(clone_df["n"], clone_df["perc"])]
        exp_dict = dict(zip(clone_df.index.values, clone_df["expanded"]))
        clone_exp_dict[cat] = exp_dict

    # Add the prop info to the data
    adata.obs["clone_percent"] = [clone_perc_dict[subj][clone] if clone in clone_perc_dict[subj].keys()
                                      else "" for subj, clone in zip(adata.obs[sample_col], adata.obs["clone"])]
    adata.obs["clone_percent"] = pd.to_numeric(adata.obs["clone_percent"])

    # Add the expanded info
    adata.obs["expanded"] = [clone_exp_dict[subj][clone] if clone in clone_exp_dict[subj].keys()
                           else False for subj, clone in zip(adata.obs[sample_col], adata.obs["clone"])]




# This function creates bar plots for the frequencies of cells that have chain data for every cluster
# adata - an anndata object with repertoire data in obs
# cluster_col - the name os the cluster column in obs
def plot_chains(adata, cluster_col):
    def get_perc_out_of_cluster(adata, cluster_col, feature):
        df_feature = adata.obs[[cluster_col, feature]].value_counts().reset_index()
        df_feature = df_feature.loc[df_feature[feature]==True,[cluster_col,"count"]].copy().rename(columns={"count":"total_"+feature})
        df_cluster = adata.obs[cluster_col].value_counts().reset_index().rename(columns={"count":"total"})
        df = df_feature.merge(df_cluster, on = cluster_col)
        df["prop_"+feature] = (df["total_"+feature]*100)/df["total"]
        return(df)

    df_alpha = get_perc_out_of_cluster(adata, cluster_col=cluster_col, feature="has_alpha")
    df_beta = get_perc_out_of_cluster(adata, cluster_col=cluster_col, feature="has_beta")
    df_tcr = get_perc_out_of_cluster(adata, cluster_col=cluster_col, feature="has_tcr")
    df = df_alpha.merge(df_beta, on=cluster_col).merge(df_tcr, on=cluster_col)

    df = df[[cluster_col, "prop_has_alpha", "prop_has_beta", "prop_has_tcr"]].melt(id_vars=cluster_col,value_name="Percent", var_name="Feature")
    df = df.replace({"prop_has_beta": "Beta", "prop_has_alpha": "Alpha", "prop_has_tcr": "TCR"})
    df[cluster_col] = df[cluster_col].astype(str)

    # Setting cluster order
    clust_order = adata.obs[cluster_col].cat.categories
    df[cluster_col] = pd.Categorical(df[cluster_col], ordered=True, categories=clust_order)

    # Plotting barplots - per  cluster
    with PdfPages("TCR_stats_per_cluster.pdf") as pdf:
        plt.clf()
        plt.figure(figsize=(5, 5))
        sns.set(font_scale=1)
        sns.set_style("white")
        g = sns.barplot(data=df, y=cluster_col, x="Percent", hue="Feature", orient='h', linewidth=0, edgecolor="white")#, width=0.9
        g.set_xticks([0,25,50,75,100])
        plt.legend(ncol=1, loc='center left', bbox_to_anchor=(1, 0.5))
        plt.tight_layout()
        pdf.savefig(g.figure)




# Running permutation test and ploting for all conditions together and for every condition separately
# by_cond_or_samp - should the enrichment be tested byu condition ("Condition") or sample ("Sample")
def test_enrichment_expanded_per_cluster(T_data, n_perm, title, onlyPlot=False, by_cond_or_samp="Condition", cluster_col = "sub_lin_2"):
    if (onlyPlot==False):
        # Filtering to cells that have TCR info
        T_data_filt = T_data["GEX"][T_data["GEX"].obs["has_tcr"] == True, :].copy()
        # df includes one row per clone
        df = T_data_filt.obs[["Sample", cluster_col, "clone", "expanded"]].copy()

        results = {}
        if (by_cond_or_samp=="Condition"):
            print("** Run per condition **")
            # run "All" once, then per-condition in a loop
            results["All"] = permutation_test_cluster_enrichment_clone_weighted(df, sample_col="Sample",cluster_col=cluster_col, clone_col="clone", expanded_col="expanded", n_perm=n_perm, seed=1)

            samp_cond = pd.read_csv(code_path+"/AuxFiles/sample_condition_list.csv")
            for cond in samp_cond["Condition_group"].dropna().unique():
                print(cond)
                curr_samp = samp_cond.loc[samp_cond.Condition_group == cond, "Sample"]
                df_sub = df.loc[df["Sample"].isin(curr_samp), :].copy()
                results[str(cond)] = permutation_test_cluster_enrichment_clone_weighted(
                    df_sub,
                    sample_col="Sample",
                    cluster_col=cluster_col,
                    clone_col="clone",
                    expanded_col="expanded",
                    n_perm=n_perm,
                    seed=1
                )
        else:
            print("** Run per Sample **")
            samples = list(set(df["Sample"]))
            for curr_samp in samples:   # curr_samp = 'BMP-0135'
                print(curr_samp)
                df_sub = df.loc[df["Sample"]==curr_samp, :].copy()
                results[curr_samp] = permutation_test_cluster_enrichment_clone_weighted(
                    df_sub,
                    sample_col="Sample",
                    cluster_col=cluster_col,
                    clone_col="clone",
                    expanded_col="expanded",
                    n_perm=n_perm,
                    seed=1
                )

        # save all results in one long CSV
        all_long = []
        if (by_cond_or_samp=="Condition"):
            feature_title = "Condition"
            cvs_name = "perm_test_" + str(n_perm) + "_all_and_by_condition_" + title + ".csv"
        else:
            feature_title="Sample"
            cvs_name = "perm_test_" + str(n_perm) + "_all_and_by_Sample" + title + ".csv"

        for feature, res in results.items():
            tmp = res.reset_index()  # has columns: cluster, stat_z, p_perm_two_sided
            tmp[feature_title] = feature
            all_long.append(tmp)

        results_long = pd.concat(all_long, ignore_index=True)

        # multiple hypothesis correction - by cell subset
        for clust in results_long["cluster"].unique():
            mask = results_long["cluster"] == clust
            curr_pvals = results_long.loc[mask, "p_perm_two_sided"]
            results_long.loc[mask, "p_adj_by_cluster"] = multipletests(curr_pvals, method="fdr_bh")[1]

        results_long.to_csv(cvs_name, index=False)
    else:
        if(feature_title == "Sample"):
            results_long = pd.read_csv("perm_test_" + str(n_perm) + "_all_and_by_Sample" + title + ".csv")
        else:
            results_long = pd.read_csv("perm_test_" + str(n_perm) + "_all_and_by_condition_" + title + ".csv")

    plot_expansion_enrichment(results_long, title, feature=by_cond_or_samp)





# Feature can be either condition or sample
def plot_expansion_enrichment(plot_df, title, feature="condition", cluster_order=None, feature_order=None, figsize=None, fontsize=20, scale_variable=1.0):
    if (type(plot_df)==str):
        plot_df = pd.read_csv(plot_df)

    if (cluster_order != None):
        plot_df["cluster"] = plot_df["cluster"].astype(str)
        plot_df["cluster"] = pd.Categorical(plot_df["cluster"], categories=cluster_order, ordered=True)
    if (feature_order != None):
        plot_df[feature] = plot_df[feature].astype(str)
        plot_df[feature] = pd.Categorical(plot_df[feature], categories=feature_order, ordered=True)

    # assign each feature to an x-position (one vertical column per feature)
    if (feature_order==None):
        feature_order = list(set(plot_df[feature]))  # ["All", cond1, cond2, ...]
        if (feature=="Sample"):
            feature_order = sorted(feature_order, key=lambda s: int(re.search(r"\d+", s).group()))
    x_map = {c: i for i, c in enumerate(feature_order)}
    plot_df["x"] = plot_df[feature].map(x_map)

    # dot size ~ -log10(p)
    plot_df["size"] = -np.log10(plot_df["p_adj_by_cluster"].clip(lower=1e-300))
    sizes = (30 + 200 * (plot_df["size"] / plot_df["size"].max())) * scale_variable

    plt.close("all")
    if (figsize==None):
        fig, ax = plt.subplots(figsize=(max(6, 1.2 * len(feature_order) + 2), 0.25 * len(set(plot_df["cluster"])) + 2))
    else:
        if (feature=="Sample"): # add bottom axis with annotations of samples
            fig, axes = plt.subplots(nrows=2, ncols=1 ,figsize=figsize, sharex=True,   gridspec_kw={"height_ratios": [4, 1.5], "hspace": 0.05})
            ax = axes[0]
            ax_bottom = axes[1]
            ax_bottom = create_annotation_bottom_plot(ax_bottom, samples=feature_order)
        else:
            fig, ax = plt.subplots(figsize=figsize)
    y = plot_df["cluster"].cat.codes.to_numpy()
    z = plot_df["stat_z"].to_numpy()
    norm = TwoSlopeNorm(vmin=z.min(), vcenter=0.0, vmax=z.max())  # force white at zero (global across all columns so colors comparable)
    lw = np.where(plot_df["p_adj_by_cluster"].to_numpy() < 0.05, 1.8, 0.3)
    sc = ax.scatter(plot_df["x"].to_numpy(), y, c=z, s=sizes, cmap="RdBu_r", norm=norm, edgecolors="k", linewidths=lw)

    ax.set_yticks(np.arange(len(cluster_order)))
    ax.set_yticklabels(cluster_order, fontsize=fontsize * scale_variable)
    ax.invert_yaxis()
    ax.set_xticks(np.arange(len(feature_order)))
    ax.set_xticklabels(feature_order, rotation=90, ha="center", fontsize=fontsize * scale_variable)
    pad = 0.2
    ax.set_ylim(len(cluster_order) - 0.5 + pad, -0.5 - pad)
    ax.set_xlim(-0.5 - pad, len(feature_order) - 0.5 + pad)

    ax.grid(False)
    ax.set_xlabel("")
    ax.set_ylabel("")
    fig.suptitle(title, fontsize=fontsize * scale_variable)

    # Colorbar
    if (feature=="Sample"): # make space for the bottom axis
        cbar = fig.colorbar(sc, ax=[ax, ax_bottom], label="Z score")
    else:
        cbar = fig.colorbar(sc, ax=ax, label="Z score")
    cbar.ax.tick_params(labelsize=fontsize * scale_variable)

    # size legend using representative p-values
    p_examples = [0.05, 0.01, 0.001]
    handles = []
    labels = []
    size_raw = plot_df["size"].to_numpy()
    size_max = size_raw.max()

    for p_ex in p_examples:
        # same mapping as the plot: size_raw = -log10(p), then scaled to 30..230
        s_raw = -np.log10(max(p_ex, 1e-300))
        s_ex = (30 + 200 * (s_raw / size_max)) * scale_variable
        h = ax.scatter([], [], s=s_ex, c="lightgray", edgecolors="k", linewidths=0.8)
        handles.append(h)
        labels.append(f"p_adj = {p_ex:g}")

    ax.legend(handles, labels, title="Significance (adj_P)", frameon=True, loc="upper left",
              bbox_to_anchor=(1.5, 1.0), borderaxespad=0.,
              fontsize=fontsize * scale_variable, title_fontsize=fontsize * scale_variable)
    plt.tight_layout()
    fig.savefig("TCR_expansion_enrich_" + title + ".pdf")




def observed_stats_clone_weighted_from_cells(df, sample_col="Sample", cluster_col="cluster", clone_col="clone",    expanded_col="expanded"):
    """
    Compute a clone-level enrichment statistic for each cluster.

    Input
    -----
    df : One row per cell. Must contain: sample_col, cluster_col, clone_col,  expanded_col (0/1 or bool)

    Rationale
    ---------
    For each clone c in sample s:
      - n_c       = total number of cells in that clone
      - n_ck      = number of cells from that clone found in cluster k
      - frac_ck   = n_ck / n_c  (fraction of the clone that lies in cluster k)
    Each clone contributes to cluster k with weight:
      w_ck = frac_ck * sqrt(n_c)
    For each sample s:
      p_s = fraction of clones in that sample that are expanded
    For each sample s and cluster k:
      obs_sk = sum of w_ck over expanded clones only
      n_eff_sk = sum of w_ck over all clones
      exp_sk = n_eff_sk * p_s
      var_sk = sum(w_ck^2) * p_s * (1 - p_s)

    Final cluster statistic:
      stat_k = sum_s (obs_sk - exp_sk) / sqrt(sum_s var_sk)

    Output
    ------
        Index = cluster names
        Values = z-like enrichment statistic per cluster
    """

    # Keep only needed columns and convert expanded flag to integer 0/1
    x = df[[sample_col, clone_col, cluster_col, expanded_col]].copy()
    x[expanded_col] = x[expanded_col].astype(int)

    # ------------------------------------------------------------------
    # clone_level: one row per (sample, clone)
    # ------------------------------------------------------------------
    # expanded_clone = whether this clone is expanded
    # If expanded is already clone-level, max() is just a safe way to collapse
    clone_level = (x.groupby([sample_col, clone_col], observed=True)[expanded_col].max().rename("expanded_clone").reset_index())

    # ------------------------------------------------------------------
    # cc: one row per (sample, clone, cluster)
    # ------------------------------------------------------------------
    # n_ck = number of cells from clone c found in cluster k
    cc = (x.groupby([sample_col, clone_col, cluster_col], observed=True).size().reset_index(name="n_ck"))

    # clone_sizes: one row per (sample, clone)
    # clone_size = total number of cells in that clone
    clone_sizes = (cc.groupby([sample_col, clone_col], observed=True)["n_ck"].sum().rename("clone_size").reset_index())

    # Add clone size and expanded status onto the clone-cluster table
    cc = cc.merge(clone_sizes, on=[sample_col, clone_col], how="left")
    cc = cc.merge(clone_level, on=[sample_col, clone_col], how="left")

    # ------------------------------------------------------------------
    # Per-clone contribution to each cluster
    # ------------------------------------------------------------------
    # frac = fraction of clone cells in this cluster
    cc["frac"] = cc["n_ck"] / cc["clone_size"]

    # w = clone contribution to this cluster
    #   = fraction of clone in cluster * dampened clone size
    cc["w"] = cc["frac"] * np.sqrt(cc["clone_size"])

    # sample-level fraction of expanded clones
    # p_s[s] = fraction of clones in sample s that are expanded
    p_s = clone_level.groupby(sample_col, observed=True)["expanded_clone"].mean()

    # ------------------------------------------------------------------
    # Terms used for the observed statistic and variance
    # ------------------------------------------------------------------
    # obs_part = contribution from expanded clones only
    cc["obs_part"] = cc["w"] * cc["expanded_clone"]

    # w2 = squared weights, used in the variance
    cc["w2"] = cc["w"] ** 2

    # ------------------------------------------------------------------
    # Aggregate within each sample x cluster
    # ------------------------------------------------------------------
    # obs    = observed total contribution of expanded clones
    # n_eff  = total contribution of all clones
    # sum_w2 = sum of squared clone weights
    agg_sk = (cc.groupby([sample_col, cluster_col], observed=True).agg(obs=("obs_part", "sum"),n_eff=("w", "sum"),sum_w2=("w2", "sum")).reset_index())

    # Add sample-level expanded-clone frequency
    agg_sk["p_s"] = agg_sk[sample_col].map(p_s)

    # Expected expanded-clone contribution under the null
    agg_sk["exp"] = agg_sk["n_eff"] * agg_sk["p_s"]

    # Variance under the null
    agg_sk["var"] = agg_sk["sum_w2"] * agg_sk["p_s"] * (1 - agg_sk["p_s"])

    # ------------------------------------------------------------------
    # Aggregate across samples to get one statistic per cluster
    # ------------------------------------------------------------------
    agg = agg_sk.groupby(cluster_col)[["obs", "exp", "var"]].sum()

    # z-like statistic
    agg["stat"] = (agg["obs"] - agg["exp"]) / np.sqrt(agg["var"] + 1e-12)

    return agg["stat"]


def permutation_test_cluster_enrichment_clone_weighted(
    df,
    n_perm=2000,
    seed=0,
    sample_col="sample",
    cluster_col="cluster",
    clone_col="clone",
    expanded_col="expanded"
):
    """
    Permute cell-to-cluster assignment within sample,
    but compute enrichment statistic at the clone level.
    """
    df = df[[sample_col, cluster_col, clone_col, expanded_col]].copy()
    df[expanded_col] = df[expanded_col].map({True: 1, False: 0}).fillna(df[expanded_col]).astype(int)
    rng = np.random.default_rng(seed)

    obs = observed_stats_clone_weighted_from_cells(
        df,
        sample_col=sample_col,
        cluster_col=cluster_col,
        clone_col=clone_col,
        expanded_col=expanded_col
    )
    clusters = obs.index.to_numpy()
    obs_vals = obs.to_numpy()

    perm_stats = np.zeros((n_perm, len(clusters)), dtype=float)

    for b in range(n_perm):
        dfp = permute_cluster_labels_within_sample(
            df, rng,
            sample_col=sample_col,
            cluster_col=cluster_col
        )
        st = observed_stats_clone_weighted_from_cells(
            dfp,
            sample_col=sample_col,
            cluster_col=cluster_col,
            clone_col=clone_col,
            expanded_col=expanded_col
        )
        perm_stats[b, :] = st.reindex(clusters).to_numpy()

    abs_obs = np.abs(obs_vals)
    abs_perm = np.abs(perm_stats)
    pvals = (1.0 + (abs_perm >= abs_obs[None, :]).sum(axis=0)) / (n_perm + 1.0)

    out = pd.DataFrame({
        "cluster": clusters,
        "stat_z": obs_vals,
        "p_perm_two_sided": pvals,
    }).set_index("cluster").sort_values("p_perm_two_sided")

    return out



def permute_cluster_labels_within_sample(df, rng, sample_col="sample", cluster_col="cluster"):
    """
    Returns a copy of df with cluster labels permuted within each sample.
    Preserves sample-level cluster counts exactly.
    """
    dfp = df.copy()
    for s, idx in dfp.groupby(sample_col).groups.items():
        labels = dfp.loc[idx, cluster_col].to_numpy()
        rng.shuffle(labels)
        dfp.loc[idx, cluster_col] = labels
    return dfp




# This function plots a given feature in a bosplot split by condition as defined int he condition-sampl convertion file
def boxplot_feature(df, feature, title):
    # Collecting daat
    samp_cond = pd.read_csv(code_path+"/AuxFiles/sample_condition_list.csv")
    cond_grp = list(set(samp_cond["Condition_group"]))
    df_to_plot = pd.DataFrame()
    for cond in cond_grp:
        curr_samp = list(samp_cond.loc[samp_cond["Condition_group"]==cond, "Sample"])
        curr = df.loc[df["Sample"].isin(curr_samp), ["Sample", feature]]
        curr["Cond_group"] = cond
        df_to_plot = pd.concat([df_to_plot, curr])

    # Plotting
    condition_colors = dict(zip(["Cont_labored", "Cont_unlabored", 'Cont_PTB', "DM1", "FGR", "PrePE", "TermPE", "PTB"],
                                ["#489142", "#489142", "#489142", "#BF0A3A", "#BF0A3A", "#BF0A3A", "#BF0A3A", "#BF0A3A"]))
    df_to_plot["Cond_group"] = pd.Categorical(df_to_plot["Cond_group"], ordered=True, categories=["Cont_labored", "Cont_unlabored", 'Cont_PTB', "PTB", "PrePE", "TermPE", "DM1", "FGR"]) # order of conditions
    plt.close("all")
    fig, ax = plt.subplots(1, 1, figsize=(6, 5))
    sns.boxplot(data = df_to_plot, x = "Cond_group", y= feature, hue="Cond_group", palette=condition_colors, ax = ax)
    sns.stripplot(data = df_to_plot, x = "Cond_group", y= feature, color="black", ax = ax)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    ax.set_title(feature, fontweight='bold')
    plt.tight_layout()
    fig.savefig(title+".pdf")













def compute_barplot_summary(T_data):
    df = T_data.obs[["short_sublin2", "expanded"]].value_counts().reset_index()
    expanded = df.loc[df["expanded"], ["short_sublin2", "count"]].set_index("short_sublin2")
    total = pd.DataFrame(T_data.obs["short_sublin2"].value_counts())
    df_perc_expand = expanded * 100 / total

    df_summ = pd.merge(total, df_perc_expand, left_index=True, right_index=True)
    df_summ.columns = ["Total_cells_with_TCRs", "Percent_expanded"]

    df_summ = df_summ.loc[~df_summ.index.isin(["mT_13_gd_DN","mT_14_gd_CD8","mT_10_gd_eff","T_11_gd_naive"]),:].copy()

    df_summ["CD8_CD4"] = df_summ.index.map({
        "mT_4_CD8_eff": "CD8",
        "mT_5_cycling": "CD8",
        "mT_9_CD4_eff": "CD4",
        "mT_2_CD8_eff": "CD8",
        "mT_6_Treg": "CD4",
        "mT_7_CD8_eff": "CD8",
        "mT_15_CD4": "CD4",
        "mT_12_MAIT": "CD8",
        "mT_1_CD4_eff": "CD4",
        "fT_8_CD8": "CD8",
        "fT_3_CD4": "CD4"
    })

    df_summ.sort_values(["CD8_CD4", "Percent_expanded"], inplace=True, ascending=False)
    return df_summ


def barplot_expanded_on_axes(T_data, title, ax_top, ax_bottom, show_xlabels=True):
    df_summ = compute_barplot_summary(T_data)
    split_idx = np.where(df_summ["CD8_CD4"].values[:-1] != df_summ["CD8_CD4"].values[1:])[0][0] + 0.5
    x = np.arange(len(df_summ))

    ax_top.bar(x, df_summ["Total_cells_with_TCRs"])
    ax_top.set_ylabel("Total cells", fontsize=7)
    ax_top.grid(False)
    ax_top.axvline(split_idx, color="black", lw=0.8)
    ax_top.set_title(title, fontsize=9, fontweight="bold", pad=2)
    ax_top.spines[['top', 'right']].set_visible(False)
    ax_top.tick_params(axis="x", bottom=False, labelbottom=False)
    ax_top.tick_params(axis="y", labelsize=6)

    top_max = nice_top_tick(df_summ["Total_cells_with_TCRs"].max())
    ax_top.set_yticks([0, top_max // 2, top_max])

    ax_bottom.bar(x, df_summ["Percent_expanded"])
    ax_bottom.set_ylabel("% Expanded", fontsize=7)
    ax_bottom.grid(False)
    ax_bottom.axvline(split_idx, color="black", lw=0.8)
    ax_bottom.spines[['top', 'right']].set_visible(False)
    ax_bottom.tick_params(axis="y", labelsize=6)

    bottom_max = nice_top_tick(df_summ["Percent_expanded"].max())
    ax_bottom.set_yticks([0, bottom_max // 2, bottom_max])

    ax_bottom.set_xticks(x)
    if show_xlabels:
        ax_bottom.set_xticklabels(df_summ.index, rotation=90, fontsize=5)
    else:
        ax_bottom.set_xticklabels([])

    return df_summ


def nice_top_tick(x):
    if x <= 10:
        step = 1
    elif x <= 100:
        step = 10
    elif x <= 1000:
        step = 100
    else:
        step = 1000
    return int(np.ceil(x / step) * step)


def plot_expanded_grid(T_data, samp_cond_conv, out_pdf, nrows=3, ncols=3):
    conditions = sorted(samp_cond_conv["Condition_group"].dropna().unique())

    fig = plt.figure(figsize=(6, 5), constrained_layout=True)   # half-page-ish
    outer = fig.add_gridspec(nrows, ncols, wspace=0.45, hspace=0.7)

    # first panel = all conditions
    panels = [("All conditions", T_data)]

    # add each condition
    for cond in conditions:
        curr_samp = list(samp_cond_conv.loc[samp_cond_conv["Condition_group"] == cond, "Sample"])
        curr_data = T_data[T_data.obs["Sample"].isin(curr_samp), :].copy()
        panels.append((cond, curr_data))

    # keep only up to nrows*ncols panels
    panels = panels[:nrows * ncols]

    for i, (title, curr_data) in enumerate(panels):
        r, c = divmod(i, ncols)
        inner = outer[r, c].subgridspec(2, 1, hspace=0.05, height_ratios=[1, 1])
        ax_top = fig.add_subplot(inner[0])
        ax_bottom = fig.add_subplot(inner[1], sharex=ax_top)

        barplot_expanded_on_axes(curr_data, title=title,ax_top=ax_top, ax_bottom=ax_bottom, show_xlabels=True)

    fig.savefig(out_pdf)
    plt.close(fig)










# bottom annotation strip
def create_annotation_bottom_plot(ax_ann, samples):
    # condition matrix for multi-label annotation
    samp_cond = pd.read_csv(code_path + "/AuxFiles/sample_condition_list.csv")
    condition_colors = dict(zip(["Cont_labored", "Cont_unlabored", 'Cont_PTB', "DM1", "FGR", "PrePE", "TermPE", "PTB"],
                                ["#CCEA8D", "#A6BC09", "green", "#5DB3D8", "#FFC817", "#EF7E29", "#E53938", "#8F2F82"]))

    mat = (samp_cond.assign(val=1).pivot_table(index="Condition_group", columns='Patient_ID', values="val", fill_value=0).reindex(columns=samples).reindex(condition_colors.keys()))

    for r, grp in enumerate(mat.index):
        color = condition_colors[grp]
        for c, s in enumerate(samples):
            if mat.loc[grp, s] == 1:
                ax_ann.add_patch(plt.Rectangle((c - 0.5, r - 0.5), 1, 1, color=color))
            else:
                ax_ann.add_patch(plt.Rectangle((c - 0.5, r - 0.5), 1, 1, color="white", ec="lightgrey", lw=0.5))
    ax_ann.set_xlim(-0.5, len(samples) - 0.5)
    ax_ann.set_ylim(len(mat) - 0.5, -0.5)
    ax_ann.set_yticks(range(len(mat.index)))
    ax_ann.set_yticklabels(mat.index)
    ax_ann.set_xticks(np.arange(len(samples)))
    ax_ann.set_xticklabels(samples, rotation=90)
    ax_ann.spines[['top','right']].set_visible(False)
    ax_ann.grid(False)
    return(ax_ann)
