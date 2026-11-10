import muon as mu
import os
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import sys
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib
import scipy
sc.set_figure_params(vector_friendly=True, dpi_save=300)
matplotlib.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']

data_dir="/projects/Pregnancy_placenta/Results/21_Final_objects/"
code_dir = "/projects/Pregnancy_placenta/Code/Publication_code/"
sys.path.append(code_dir)
import FigFunctions as func
import Repertoire_functions as rep_func

os.chdir("/projects/Pregnancy_placenta/Results/22_Final_figures/SuppFigure_7/")


#===============================================
# Sup Fig 10A) TCR chains per cluster
#===============================================
t_data = mu.read(data_dir + "T_cells.h5mu")
t_data["GEX"].obs["sub_lin2"] =  pd.Categorical(t_data["GEX"].obs["sub_lin2"], ordered=True, categories=[ 'mT_1_CD4_eff', 'mT_2_CD8_eff_GZMK',  'fT_3_CD4_naive','mT_4_CD8_eff_LAG3', 'mT_5_cycling', 'mT_6_Treg', 'mT_7_CD8_eff_GZMB','fT_8_CD8_naive', 'mT_9_CD4_eff_CTLA4', 'mT_10_gd_TRDV1_eff','T_11_gd_TRDV1_naive',  'mT_12_MAIT', 'mT_13_gd_DN_TRDV2-9',
 'mT_14_gd_CD8_TRDV2-9', 'mT_15_CD4_CXCL13'])
rep_func.plot_chains(adata=t_data["GEX"], cluster_col="sub_lin2")


#=======================================================================================================================
# Sup Fig 10B) General TCR stats
#  Total T cells, Total cells with TCRs, Plot proportions of fetal / maternal clones per sample
#=======================================================================================================================
# Count maternal and fetal T cells
#-----------------------------------
t_data = mu.read(data_dir + "T_cells.h5mu")
t_gex = t_data["GEX"]

df = t_gex.obs[['Patient_ID', 'FetMat_Annotation']].value_counts().reset_index()
df2 = t_gex.obs[['Patient_ID', 'FetMat_Annotation', 'has_tcr']].value_counts().reset_index()
df2 = df2.loc[df2["has_tcr"]==True,:].copy()
df2.rename(columns={'count':'cells_wTCR'}, inplace=True)

df3 = t_gex.obs[['Patient_ID', 'FetMat_Annotation', 'expanded']].value_counts().reset_index()
df3 = df3.loc[df3["expanded"]==True,:].copy()
df3.rename(columns={'count':'cells_expanded'}, inplace=True)

df = df.merge(df2[['Patient_ID', "FetMat_Annotation", "cells_wTCR"]], how="left", on=['Patient_ID', "FetMat_Annotation"])
df = df.merge(df3[['Patient_ID', "FetMat_Annotation", "cells_expanded"]], how="left", on=['Patient_ID', "FetMat_Annotation"])

df = df.pivot_table(index='Patient_ID', columns = ["FetMat_Annotation"], values=["cells_wTCR", "count","cells_expanded"])
df.columns = ['_'.join([str(x) for x in col if x not in (None, '')]).strip() for col in df.columns.to_flat_index()]
samp_cond = t_gex.obs[['Patient_ID',"Condition"]].drop_duplicates().reset_index(drop=True)
df = df.merge(samp_cond, how="left", on='Patient_ID')
df = df[['Patient_ID',"Condition", "count_Maternal","count_Fetal","count_nan",
                               "cells_wTCR_Maternal", "cells_wTCR_Fetal","cells_wTCR_nan",
                               "cells_expanded_Maternal","cells_expanded_Fetal","cells_expanded_nan"]]

num_cols = df.select_dtypes(exclude="category").columns
df[num_cols] = df[num_cols].fillna(0)
df["total"] = df["count_Maternal"] + df["count_Fetal"]+ df["count_nan"]
df["total_withTCR"] = df["cells_wTCR_Maternal"] + df["cells_wTCR_Fetal"]+ df["cells_wTCR_nan"]
df["per_fetal"] = df["cells_wTCR_Fetal"]*100/ df["total_withTCR"]
df["per_maternal"] = df["cells_wTCR_Maternal"]*100/ df["total_withTCR"]
df["per_unassigned"] = df["cells_wTCR_nan"]*100/ df["total_withTCR"]

df.to_csv("FetMat_TCR_counts_v3.csv", index=False)



# Plot sample level T cell data
#-----------------------------------
import importlib
rep_func = importlib.reload(rep_func)

from matplotlib.gridspec import GridSpec
df = pd.read_csv("FetMat_TCR_counts_v3.csv")

samp_cond = pd.read_csv(code_dir+"AuxFiles/sample_condition_list.csv")

# sample order by original Condition, then Sample
cond_order = ["Cont_labored", "Cont_unlabored", "Cont_PTB", "DM1", "FGR", "PrePE", "TermPE", "PTB"]
sample_order = (
    samp_cond.assign(Condition_group=pd.Categorical(samp_cond["Condition_group"], categories=cond_order, ordered=True))
    .sort_values(["Condition_group", 'Patient_ID'])
    ['Patient_ID']
    .drop_duplicates())

df2 = df.set_index('Patient_ID').loc[sample_order].reset_index()
df2.to_csv("Supp_fig6B_FetMat_TCR_counts_per_samp.csv", index=False)
samples = df2['Patient_ID'].tolist()
x = np.arange(len(samples))

fig = plt.figure(figsize=(16, 7))
gs = GridSpec(4, 1, height_ratios=[1, 1, 2.2, 1.2], hspace=0.08)
ax_top = fig.add_subplot(gs[0])
ax_mid = fig.add_subplot(gs[1], sharex=ax_top)
ax_main = fig.add_subplot(gs[2], sharex=ax_top)
ax_ann = fig.add_subplot(gs[3], sharex=ax_top)

# top: total cells
ax_top.bar(x, df2["total"], color="grey")
ax_top.set_ylabel("Total T cells")
ax_top.spines[['top','right']].set_visible(False)
ax_top.tick_params(axis="x", bottom=False, labelbottom=False)

# middle: cells with TCR
ax_mid.bar(x, df2["total_withTCR"], color="dimgray")
ax_mid.set_ylabel("Cells\nwith TCR")
ax_mid.spines[['top','right']].set_visible(False)
ax_mid.tick_params(axis="x", bottom=False, labelbottom=False)

# main stacked percentages
ax_main.bar(x, df2["per_maternal"], color="#B58B47", label="Maternal")
ax_main.bar(x, df2["per_fetal"], bottom=df2["per_maternal"], color="black", label="Fetal")
ax_main.bar(x, df2["per_unassigned"], bottom=df2["per_maternal"]+df2["per_fetal"], color="lightgrey", label="Unassigned")
ax_main.set_ylabel("Percent")
ax_main.legend(frameon=False)
ax_main.spines[['top','right']].set_visible(False)
ax_main.tick_params(axis="x", bottom=False, labelbottom=False)
ax_main.legend(frameon=False, loc="center right", bbox_to_anchor=(-0.08, 0.5))

# separators by original Condition
for i in range(1, len(df2)):
    if df2.loc[i, "Condition"] != df2.loc[i-1, "Condition"]:
        for ax in [ax_top, ax_mid, ax_main, ax_ann]:
            ax.axvline(i - 0.5, color="k", lw=0.8)

# Bottom: annotation of conditions
ax_ann = rep_func.create_annotation_bottom_plot(ax_ann, samples)

ax_top.grid(False)
ax_mid.grid(False)
ax_main.grid(False)

plt.tight_layout()
plt.savefig("SuppFig6B_Tcell_data_per_samp.pdf", dpi=300)



#=============================================================================================
# Sup Fig 10C) Barplot of clone proportions per sample
#=============================================================================================
t_data= mu.read(data_dir+"/T_cells.h5mu")["GEX"]
gex_data=t_data[t_data.obs["has_tcr"],].copy()
samp_cond = pd.read_csv(code_dir+"AuxFiles/sample_condition_list.csv")

# sample order
cond_order = ["Cont_labored", "Cont_unlabored", "Cont_PTB", "DM1", "FGR", "PrePE", "TermPE", "PTB"]
sample_order = (samp_cond.assign(Condition_group=pd.Categorical(samp_cond["Condition_group"],categories=cond_order,ordered=True))
    .sort_values(["Condition_group", 'Patient_ID'])['Patient_ID']
    .drop_duplicates())

# count cells per clone within each sample
clone_sizes = (gex_data.obs.groupby(['Patient_ID', "clone"]).size().reset_index(name="count"))

# keep only wanted sample order
sample_order = [s for s in sample_order if s in clone_sizes['Patient_ID'].unique()]

# build per-sample clone-size lists
clone_lists = {}
singleton_sum = {}

for s in sample_order:
    counts = clone_sizes.loc[clone_sizes['Patient_ID'] == s, "count"]
    clone_lists[s] = sorted(counts[counts > 1].tolist(), reverse=True)
    singleton_sum[s] = int((counts == 1).sum())

# plot
fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(14, 5), gridspec_kw={"height_ratios": [4, 1.5], "hspace": 0.05}, sharex=True)
ax=axes[0]
ax1=axes[1]
# meta
ax1 = rep_func.create_annotation_bottom_plot(ax1, samples=sample_order)

# clones
bottom = [0] * len(sample_order)

max_n = max(len(v) for v in clone_lists.values()) if clone_lists else 0
for i in range(max_n):
    heights = [clone_lists[s][i] if i < len(clone_lists[s]) else 0 for s in sample_order]
    ax.bar(range(len(sample_order)), heights, bottom=bottom, width=0.9)
    bottom = [b + h for b, h in zip(bottom, heights)]

# add singleton segment in light grey
singletons = [singleton_sum[s] for s in sample_order]
ax.bar(range(len(sample_order)), singletons, bottom=bottom, width=0.9, color="lightgrey")

ax.set_xticks(range(len(sample_order)))
ax.set_xticklabels(sample_order, rotation=90)
ax.set_ylabel("Cell count")
ax.set_xlabel('Patient_ID')
ax.grid(False)
plt.tight_layout()
fig.savefig("SuppFig6C_clone_freq_per_samp.pdf")




#=============================================================================================
# Sup Fig 10D) Shannon Enrichemnet
#=============================================================================================
t_data= mu.read(data_dir+"/T_cells.h5mu")["GEX"]

# Computing diversity only for samples that have at least 100/ 250 different clones
df = t_data.obs[["Sample","clone"]].value_counts().reset_index()["Sample"].value_counts().reset_index()
df = pd.merge(df, t_data.obs[["Sample"]].value_counts(), on="Sample", suffixes=("_clones", "_cells"))
df_min100 = df.loc[df["count_clones"]>=100,:]  # 50 samples out of 65
samp_to_analyze = list(df_min100["Sample"])
t_data = t_data[t_data.obs["Sample"].isin(samp_to_analyze),].copy() # 43,620 cells for 100 min; 41,651 cells for min250
freq_clone = t_data.obs[["Sample","clone","clone_percent"]].copy()
freq_clone.reset_index(inplace=True, drop=True)
freq_clone = freq_clone.drop_duplicates() # 30,588 clones for 100 min; 28,811 cells for min250
freq_clone["freq"] = freq_clone["clone_percent"]/100


# Computing and plotting Shannon entropy per condition
freq_clone['shannon_index'] = freq_clone['freq']*np.log(freq_clone['freq'])
df_shannon = pd.DataFrame()
for samp in samp_to_analyze: # samp="BMP-0039"
    curr = freq_clone.loc[freq_clone["Sample"]==samp, :].copy()
    shannon = -np.sum(curr['shannon_index'])
    df_shannon = pd.concat([df_shannon, pd.DataFrame([[samp, shannon]], columns = ["Sample","Shannon"])])
rep_func.boxplot_feature(df_shannon, feature="Shannon", title = "SuppFig6D_shannon")
df_shannon.to_csv("SuppFig6D_shannon.csv")

# Checking differences in Shnnon diversity across conditions
comparisons = [["Cont_labored", "Cont_unlabored"],
               ["PTB","Cont_PTB"],
               ["DM1","Cont_labored"],
               ["FGR","Cont_labored"],
               ["PrePE","Cont_labored"],
               ["TermPE","Cont_labored"],
               ["DM1","Cont_unlabored"],
               ["FGR","Cont_unlabored"],
               ["PrePE","Cont_unlabored"],
               ["TermPE","Cont_unlabored"],
               ["PTB","Cont_labored"],
               ["PrePE", "Cont_PTB"],
               ["PTB","Cont_unlabored"]]
samp_cond = pd.read_csv(code_dir+"/AuxFiles/sample_condition_list.csv")
shannon_withCond = samp_cond.merge(df_shannon, on="Sample", how="left")
# removing samples that don't have shannon entropy due to too few clones
shannon_withCond = shannon_withCond.loc[~shannon_withCond["Shannon"].isna(),:].copy()

res_shannon = pd.DataFrame()
for comp in comparisons:
    pval = scipy.stats.mannwhitneyu(x=shannon_withCond.loc[shannon_withCond["Condition_group"] == comp[0], "Shannon"],
                                    y=shannon_withCond.loc[shannon_withCond["Condition_group"] == comp[1], "Shannon"],).pvalue
    res_shannon = pd.concat([res_shannon, pd.DataFrame([["-".join(comp), pval]], index=[0])])
res_shannon.columns = ["Comparison","Pval Shannon diversity"]
res_shannon.to_csv("SuppFig6D_shannon_pvals.csv")





#=======================================================================================================================
# Sup Fig 10E) Barplots of expanded clones per cluster - total cells and %expanded per cluster
#=======================================================================================================================
os.chdir("/projects/Pregnancy_placenta/Results/22_Final_figures/SuppFigure_6/")
T_data = mu.read(data_dir + "/T_cells.h5mu")
T_data = T_data["GEX"][T_data["GEX"].obs["has_tcr"], :].copy()

samp_cond_conv = pd.read_csv("/projects/Pregnancy_placenta/Results/13_DE/sample_condition_list.csv")
rep_func.plot_expanded_grid(T_data, samp_cond_conv, out_pdf="SuppFig6E_expanded_barplot_grid.pdf", nrows=3, ncols=3)



# Sanity check - How many Tregs per condition?
T_data = mu.read(data_dir+"/T_cells.h5mu")
Treg = T_data["GEX"][T_data["GEX"].obs["sub_lin2"]=="mT_6_Treg",:].copy()
Treg_counts = pd.DataFrame()
for cond in conditions:
    curr_samp = list(samp_cond_conv.loc[samp_cond_conv["Condition_group"] == cond, "Sample"])
    curr_data = Treg[Treg.obs["Sample"].isin(curr_samp),:].copy()
    tcr_data = curr_data[curr_data.obs["has_tcr"],:].copy().copy()
    expanded = tcr_data[tcr_data.obs["expanded"],:].copy().copy()
    Treg_counts = pd.concat([Treg_counts, pd.DataFrame({"Condition":cond,
                                                        "Treg_count": curr_data.shape[0],
                                                        "Treg_wTCR": tcr_data.shape[0],
                                                        "Treg_expanded":expanded.shape[0]}, index=[0])])

Treg_counts["expanded_percent"] = Treg_counts["Treg_expanded"]*100 /  Treg_counts["Treg_wTCR"]
Treg_counts.set_index("Condition", inplace=True)
T_data_wTCR = T_data[T_data["GEX"].obs["has_tcr"],:].copy()

Treg_counts["Total_T_cells_withTCR"] = 0
for cond in conditions:
    curr_samp = list(samp_cond_conv.loc[samp_cond_conv["Condition_group"] == cond, "Sample"])
    curr_data = T_data_wTCR["GEX"][T_data_wTCR["GEX"].obs["Sample"].isin(curr_samp),:].copy()
    Treg_counts.loc[cond, "Total_T_cells_withTCR"] = curr_data.shape[0]
Treg_counts["Treg_perc_out_of_total"] = Treg_counts["Treg_count"] *100 / Treg_counts["Total_T_cells_withTCR"]
Treg_counts.loc[["PrePE","PTB","Cont_PTB","Cont_unlabored","DM1","TermPE","FGR","Cont_labored"],:]






#=======================================================================================================================
# Sup Fig 10F) Enrichment of expanded clones per sample
#=======================================================================================================================
import importlib
rep_func = importlib.reload(rep_func)

T_data = mu.read(data_dir+"/T_cells.h5mu")
# removing gd T cells for this analysis and splitting CD8 and CD4 enrichment
cd4T_data = T_data[T_data["GEX"].obs["sub_lin2"].isin(["mT_1_CD4_eff","fT_3_CD4_naive","mT_6_Treg","mT_15_CD4_CXCL13","mT_9_CD4_eff_CTLA4"])].copy()
cd8T_data = T_data[T_data["GEX"].obs["sub_lin2"].isin([ "fT_8_CD8_naive", 'mT_5_cycling',"mT_7_CD8_eff_GZMB","mT_12_MAIT" ,"mT_2_CD8_eff_GZMK","mT_4_CD8_eff_LAG3"])].copy()

# Compute enrichment
rep_func.test_enrichment_expanded_per_cluster(cd4T_data, n_perm=2000, title="CD4_v3_clones_bySample", by_cond_or_samp="Patient_ID", cluster_col ="sub_lin2")
rep_func.test_enrichment_expanded_per_cluster(cd8T_data, n_perm=2000, title="CD8_v3_clones_bySample", by_cond_or_samp="Patient_ID", cluster_col ="sub_lin2")

# Plot
cd4_clust_order = ["fT_3_CD4_naive","mT_1_CD4_eff_GZMA","mT_15_CD4_CXCL13","mT_6_Treg","mT_9_CD4_eff_CTLA4"]
cd8_clust_order = ["fT_8_CD8_naive","mT_7_CD8_eff_GZMB","mT_12_MAIT","mT_2_CD8_eff_GZMK","mT_5_cycling","mT_4_CD8_eff_LAG3"]

rep_func.plot_expansion_enrichment(plot_df="perm_test_2000_all_and_by_SampleCD4_v3_clones_bySample.csv",
                                   title="CD4_v3_clones_bySample",
                                   feature="Sample",
                                   cluster_order=cd4_clust_order,
                                   figsize=(23, 5),
                                   fontsize=15)

rep_func.plot_expansion_enrichment("perm_test_2000_all_and_by_SampleCD8_v3_clones_bySample.csv",
                                   title="CD8_v3_clones_bySample",
                                   feature="Sample",
                                   cluster_order=cd8_clust_order,
                                   figsize=(23, 5),
                                   fontsize=15)



#=============================================================================================
# Sup Fig 10G) Additional markers of checkpoint T cells
#=============================================================================================
sc.set_figure_params(vector_friendly=True, dpi_save=300)
T_data = mu.read(data_dir+"/T_cells.h5mu")

# Top proteins in mT_4
func.feature_plots_pdf(T_data["ADT"], gene_list=["P:"+x for x in ["CD38","CD39","CD103","CD49b","CD278"]], title="SuppFig6G_Tfeature_plots_prot", use_raw=False, figsize=(7,4), nrows=2, ncols=3, axes_titlesize=20, gridsize=75, hex=True, ticks = [0, 1, 2, 3], tick_labels = ["0", "1", "2", "3+"], vmax=3, cmap = "Blues")

func.feature_plots_pdf(T_data["GEX"], gene_list=["G:"+x for x in ["PRDM1","SELL","KLF2","PDCD1"]], title="SuppFig6G_Tfeature_plots_genes", use_raw=False, figsize=(7,4), nrows=2, ncols=3, axes_titlesize=20, gridsize=75, hex=True, ticks = [0, 1, 2, 3, 4, 5, 6], tick_labels = ["0", "1", "2", "3", "4", "5", "6+"], vmax=6)




#=============================================================================================
# Sup Fig 10H) Additional markers for mMac5
#=============================================================================================
sc.set_figure_params(vector_friendly=True, dpi_save=300)
mMac_data = mu.read(data_dir+"/mMac_Mon_12subsets.h5mu")

func.feature_plots_pdf(mMac_data["GEX"], gene_list=["G:"+x for x in ["IL32","IL18BP","CCL8","IL15RA","IFIT3","IFIT2","ISG20"]], title="SuppFig6H_mMac5_feature_plots_genes1", use_raw=False, figsize=(7,4), nrows=2, ncols=3, axes_titlesize=20, gridsize=75, hex=True, ticks = [0, 1, 2, 3,4], tick_labels = ["0", "1", "2", "3","4+"], vmax=4)
func.feature_plots_pdf(mMac_data["GEX"], gene_list=["G:"+x for x in ["IL15RA","HLA-DQB2","APOL4","SOCS1","CLIC2"]], title="SuppFig6H_mMac5_feature_plots_genes2", use_raw=False, figsize=(7,4), nrows=2, ncols=3, axes_titlesize=20, gridsize=75, hex=True, ticks = [0, 1, 2, 3], tick_labels = ["0", "1", "2", "3+"], vmax=3)


mMac_data["ADT"].obsm["X_umap"] = mMac_data["GEX"][mMac_data["ADT"].obs_names, :].obsm["X_umap"]
func.feature_plots_pdf(mMac_data["ADT"], gene_list=["P:"+x for x in ["CD274","CD38","CD86"]], title="SuppFig6H_mMac5_feature_plots_prot", use_raw=False, figsize=(7,4), nrows=2, ncols=3, axes_titlesize=20, gridsize=75, hex=True, ticks = [0, 1, 2], tick_labels = ["0", "1", "2+"], vmax=2, cmap = "Blues")
func.feature_plots_pdf(mMac_data["ADT"], gene_list=["P:"+x for x in ["CD86"]], title="SuppFig6H_mMac5_feature_plots_prot2", use_raw=False, figsize=(7,4), nrows=2, ncols=3, axes_titlesize=20, gridsize=75, hex=True, ticks = [0, 1], tick_labels = ["0", "1+"], vmax=1, cmap = "Blues")


#=============================================================================================
# Sup Fig 10I, J) Showing mMac5 doesn't change in DA or DE across all comparisons
#=============================================================================================
os.chdir("/projects/Pregnancy_placenta/Results/22_Final_figures/SuppFigure_6/Myeloid_summary_heatmap/")
df = pd.read_csv("/projects/Pregnancy_placenta/Results/17_Summaries/DEGs/DEGs_count_level3_pval0.15_abslog2FC0.58.csv")
df_myeloid = df.loc[df["Lineage"].isin(["DC","HBC","Mon","PAMM"]),:].copy()
df_pivot = df_myeloid.pivot(columns = "Comparison", index=["Lineage","Subset"])
df_pivot.to_excel("DEGs_myeloid.xlsx")

## ==> Go to R code

