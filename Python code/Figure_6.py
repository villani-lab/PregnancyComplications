import os
from unittest.mock import inplace

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import spearmanr
import muon as mu
import scanpy as sc
from statsmodels.stats.multitest import multipletests
import sys

pd.set_option('display.max_columns', 10)
pd.set_option("display.max_rows", 50)
pd.set_option('display.width', 400)
matplotlib.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
sc.set_figure_params(vector_friendly=True, dpi_save=400)

code_path="/projects/Pregnancy_placenta/Code/Publication_code"
sys.path.append(code_path)
import Repertoire_functions as rep_func
import FigFunctions as fig_func
os.chdir("/projects/Pregnancy_placenta/Results/22_Final_figures/Figure_6/")
data_dir = "/projects/Pregnancy_placenta/Results/21_Final_objects/"


#====================================================================================
# Define T cell clones + TCR QC plots
#====================================================================================
T_data = mu.read(data_dir+"/T_wTCR.h5mu")
T_data["GEX"].obs.pop('expanded')
T_data["GEX"].obs.pop('clone_percent')
rep_func.calc_clone_percent(T_data["GEX"], sample_col='Sample', min_cell_expanded = 2, min_prop_expanded = 0.005)



#====================================================================================
# Figure 6B) feature plots of expanded clones
#====================================================================================
samp_cond_conv = pd.read_csv(code_path+"/AuxFiles/sample_condition_list.csv")

T_data["GEX"].obs["expanded"] = pd.Categorical(T_data["GEX"].obs["expanded"].map({False: "No", True: "Yes"}))
cond_conv = {'Cont_unlabored':'Control unlabored', 'Cont_labored':'Control labored', 'Cont_PTB':'Control PTB','PTB':"sPTB", 'TermPE':'Term PE', 'PrePE': 'Preterm PE', 'DM1':'DM1', 'FGR':'FGR'}

plt.close("all")
fig, axes = plt.subplots(2,4, figsize = (7,4))
for i, cond in enumerate(['Cont_unlabored', 'Cont_labored', 'Cont_PTB','PTB', 'TermPE', 'PrePE', 'DM1', 'FGR']):
    curr_samp = list(samp_cond_conv.loc[samp_cond_conv["Condition_group"]==cond, "Sample"])
    curr_data = T_data["GEX"][T_data["GEX"].obs["Sample"].isin(curr_samp),].copy()
    sc.pl.scatter(curr_data, color="expanded", palette=["lightgrey","#65000d"], ax=axes[divmod(i, 4)],  basis="umap", show=False, size=6)
    cont_title = cond_conv[cond]
    axes[divmod(i, 4)].set_title(cont_title, size=12, fontweight="bold")
    axes[divmod(i, 4)].get_legend().remove()
    axes[divmod(i, 4)].set_xlabel("UMAP1", fontsize=7)
    axes[divmod(i, 4)].set_ylabel("UMAP2",fontsize=7)
plt.tight_layout()
fig.savefig("Fig6b_expanded_clones_per_cond.pdf")





###########################################################################
# Fig 6C) Testing enrichment of expanded clones in specific clusters
###########################################################################
import importlib
rep_func = importlib.reload(rep_func)

os.chdir("/projects/Pregnancy_placenta/Results/22_Final_figures/Figure_6/ExpClone_enrichemnet_V3_clone_level/")

T_data = mu.read(data_dir+"/T_cells.h5mu")
# removing gd T cells for this analysis and splitting CD8 and CD4 enrichment
cd4T_data = T_data[T_data["GEX"].obs["sub_lin2"].isin(["mT_1_CD4_eff","fT_3_CD4_naive","mT_6_Treg","mT_15_CD4_CXCL13","mT_9_CD4_eff_CTLA4"])].copy()
cd8T_data = T_data[T_data["GEX"].obs["sub_lin2"].isin([ "fT_8_CD8_naive", 'mT_5_cycling',"mT_7_CD8_eff_GZMB","mT_12_MAIT" ,"mT_2_CD8_eff_GZMK","mT_4_CD8_eff_LAG3"])].copy()

# Compute enrichment
rep_func.test_enrichment_expanded_per_cluster(T_data =cd4T_data, n_perm=2000, title="CD4_cloneLevel", cluster_col ="sub_lin2")
rep_func.test_enrichment_expanded_per_cluster(cd8T_data, n_perm=2000, title="CD8_cloneLevel", cluster_col ="sub_lin2")


# Plot enrichment
cd4_clust_order = ["fT_3_CD4_naive","mT_1_CD4_eff","mT_15_CD4_CXCL13","mT_6_Treg","mT_9_CD4_eff_CTLA4"]
cd8_clust_order = ["fT_8_CD8_naive","mT_7_CD8_eff_GZMB","mT_12_MAIT","mT_2_CD8_eff_GZMK","mT_5_cycling","mT_4_CD8_eff_LAG3"]
cond_order = ["All","Cont_unlabored","Cont_labored","Cont_PTB","PrePE","TermPE","PTB","FGR","DM1"]

rep_func.plot_expansion_enrichment(plot_df="perm_test_2000_all_and_by_condition_CD4_cloneLevel.csv",
                                   title="CD4_cloneLevel_byCond",
                                   feature="Condition",
                                   cluster_order=cd4_clust_order,
                                   feature_order = cond_order,
                                   figsize = (7,3.5),
                                   scale_variable=0.5)


rep_func.plot_expansion_enrichment(plot_df="perm_test_2000_all_and_by_condition_CD8_cloneLevel.csv",
                                   title="CD8_cloneLevel_byCond",
                                   feature="Condition",
                                   cluster_order=cd8_clust_order,
                                   feature_order = cond_order,
                                   figsize = (7,3.5),
                                   scale_variable=0.5)




###########################################################################
# Fig 6D) Cytotoxic and checkpoint genes
###########################################################################
import importlib
fig_func = importlib.reload(fig_func)

T_data = mu.read(data_dir+"/T_cells.h5mu")

fig_func.feature_plots_pdf(T_data["GEX"], gene_list=["G:"+x for x in ["GZMH","GZMB","NKG7","PRF1","TOX","TIGIT","CTLA4","LAG3"]], title="Fig6D_feature_plots", use_raw=False, figsize=(10,4), gridsize=75, hex=True, nrows=2, axes_titlesize=20)




###########################################################################
# Fig 6F) ISG and checkpoint and exhaustion genes
###########################################################################
import importlib
fig_func = importlib.reload(fig_func)

mMac_data = mu.read(data_dir+"/mMac_Mon_12subsets.h5mu")

fig_func.feature_plots_pdf(mMac_data["GEX"], gene_list=["G:"+x for x in ["CXCL10","CXCL9","IDO1","CXCL11"]], title="Fig6F_mMac_feature_plots_ISG", use_raw=False, figsize=(10,4), gridsize=75, hex=True, nrows=2, axes_titlesize=20)
fig_func.feature_plots_pdf(mMac_data["GEX"], gene_list=["G:"+x for x in ["CD274","ETV7","PDCD1LG2","CD38"]], title="Fig6F_mMac_feature_plots_checkpoint", use_raw=False, figsize=(10,4), gridsize=75, hex=True, nrows=2, axes_titlesize=20, ticks = [0, 1, 2, 3], tick_labels = ["0", "1", "2", "3+"], vmax=3)
fig_func.feature_plots_pdf(mMac_data["GEX"], gene_list=["G:"+x for x in ["CXCL10","CXCL9","IDO1","CXCL11","CD274","ETV7","PDCD1LG2","CD38"]], title="Fig6F_mMac_feature_plots_SC", use_raw=False, figsize=(6,4), nrows=2, axes_titlesize=20)

sc.set_figure_params(vector_friendly=True, dpi_save=300)
fig, ax = plt.subplots(figsize=(5, 5))
plt.figure(figsize=(5, 5))
sc.pl.scatter(mMac_data["GEX"], color="sub_lin2", groups=["mMac_5_CXCL10_IDO1"], palette=["red"], basis="umap", size=1, show=False, ax=ax)
for coll in ax.collections:
    coll.set_rasterized(True)
fig.savefig("Fig6E_mMac5_umap.pdf")


###########################################################################
# Fig 6L,M) Co-expression of co-inhibitory genes on myeloid and T cells
###########################################################################
import importlib
fig_func = importlib.reload(fig_func)

t_data = mu.read(data_dir + "T_cells.h5mu")["GEX"]
mMac_data = mu.read(data_dir + "mMac_Mon_12subsets.h5mu")["GEX"]
gene_list = pd.read_excel(code_path + "/AuxFiles/exhaustion_inhibition_gene_list.xlsx", sheet_name="short_list")
subset_pairs = [(t_sub, my_sub) for my_sub in ["mMac_5"] for t_sub in ["mT_4_CD8_eff", "mT_9_CD4_eff"]]

# Plot correlation scatter plots between 2 myeloid subsets and 3 T subsets
fig_func.scatter_co_expression(t_data, mMac_data, subset_pairs, gene_list=coinhib_list, pdf_title = "Fig7L_co_inhibit_myeloid_T_corr.pdf")

# Dotplot of stats
stat_df = pd.read_csv("Fig7L_stats_coinhibitory.csv")
data_for_plot = stat_df.assign(pair=stat_df["T_gene"] + "-" + stat_df["Myeloid_gene"],logp=-np.log10(stat_df["adj_pval"]))

sns.set_theme(rc={"figure.dpi": 300})
sns.set_style("ticks")
from matplotlib.colors import TwoSlopeNorm
from matplotlib.cm import ScalarMappable

plt.close("all")
plt.figure(figsize=(3.5, 5))
norm = TwoSlopeNorm(vmin=-0.7, vcenter=0, vmax=0.7)
ax = sns.scatterplot(data=data_for_plot, x="T_sub", y="pair", hue="corr", size="logp",
    palette="RdBu_r", hue_norm=norm, sizes=(40, 250), legend=False, edgecolor=(data_for_plot["adj_pval"] < 0.05).map({True: "black", False: "none"}), linewidth=2)

ax.set_xlim(-0.3, 1.3)
ax.tick_params(axis='x', labelrotation=90)

sm = ScalarMappable(cmap="RdBu_r", norm=norm)
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax, label="Correlation")
cbar.set_ticks([-0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7])

# values to SHOW in legend
pvals = [0.05, 0.01, 0.001]
# convert to the variable actually used for point size
logp_vals = -np.log10(pvals)
size_min, size_max = 40, 250
vmin, vmax = data_for_plot["logp"].min(), data_for_plot["logp"].max()

def scale_size(v):
    return size_min + (v - vmin) / (vmax - vmin) * (size_max - size_min)

handles = [Line2D([0], [0],marker="o", linestyle="",markerfacecolor="gray", markeredgecolor="black",        markersize=scale_size(lp)**0.5,label=f"{p:g}") for p, lp in zip(pvals, logp_vals)]
ax.legend(handles=handles, title="adj_pval", frameon=False, loc="upper left", bbox_to_anchor=(1.04, 0.3))
plt.tight_layout()
plt.savefig("Fig7L_Dotplot_coInhib_corr.pdf")
