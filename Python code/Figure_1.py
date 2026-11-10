
import muon as mu
import os
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
sc.set_figure_params(vector_friendly=True, dpi_save=300)
import sys
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']

data_dir="/projects/Pregnancy_placenta/Results/21_Final_objects/"
code_dir = "/projects/Pregnancy_placenta/Code/Publication_code/"
sys.path.append(code_dir)
import FigFunctions as func


_DEF_LINEAGES = ['Trophoblasts','EVT','HBC','mMac',"Monocytes","DC",'T','NK','B_Plasma','dStromal','fPerivascular','fFibroblasts', 'Endothelial','HSC','Other']
_DEF_COLORS = ["#60B364","#B2DF8A",  # trophoblasts
               "#800026","#f03b20","#fa9fb5","#E2C087","#fd8d3c","#fed976","#E80C7A",  # immune
               "#1f2b66", "#5384d4", "#5ab5ed", "#8F6FE8",  # stromal/endo
               "grey","lightgrey"]
_DEF_PALETTE = dict(zip(_DEF_LINEAGES, _DEF_COLORS))

os.chdir("/projects/Pregnancy_placenta/Results/22_Final_figures/Figure_1/")



#======================================================
# Fig1B) General UMAP
#======================================================
sc.set_figure_params(vector_friendly=True, dpi_save=300)
data = mu.read(data_dir + "AllCells.h5mu")

fig, ax = plt.subplots(1, 1, figsize=(4, 4))
sc.pl.umap(data["GEX"], color="Lineage", palette=_DEF_PALETTE, show=False, ax=ax)
handles, labels = ax.get_legend_handles_labels()
hmap = dict(zip(labels, handles))
ordered_labels  = [l for l in _DEF_LINEAGES if l in hmap]
ordered_handles = [hmap[l] for l in ordered_labels]
ax.legend(ordered_handles, ordered_labels, frameon=False)
plt.tight_layout()
fig.savefig("Fig1B_global_UMAP_legend.pdf")

fig, ax = plt.subplots(1, 1, figsize=(4, 4))
sc.pl.umap(data["GEX"], color="Lineage", palette=_DEF_PALETTE, show=False, legend_loc="upper right", ax=ax)
fig.savefig("Fig1B_global_UMAP.pdf")


#========================================================
# Fig1C-F) UMAPs of different lineages with demultiplexing
#========================================================
os.chdir("/projects/Pregnancy_placenta/Results/22_Final_figures/Figure_1/UMAP_lineages")
sc.set_figure_params(dpi=500, vector_friendly=True, fontsize=6)
plt.rcParams['axes.labelsize'] = 5

data_objects = {"AllTroph": [data_dir+"/All_Troph_EVT.h5mu", None],
                "VillousTroph": [data_dir+"/Trophoblasts.h5mu", None],
                "EVT": [data_dir + "/EVT.h5mu", None],
                "Fetal_PV":	[data_dir+"/fPV.h5mu", None],
                "Fetal_fib": [data_dir+"fFB.h5mu", None],
                "Decidual_stromal": [data_dir+"dS.h5mu", None],
                "Endothelial": [data_dir+"/Endo.h5mu", ['fEndo_1','fEndo_2','fEndo_3','fEndo_4','mEndo_1','mEndo_2','mEndo_3','mEndo_lymph','fEndo_cycling']],
                "HBC": [data_dir+"/HBC.h5mu", None],
                "HBC_13subsets": [data_dir+"/HBC_13subsets.h5mu", ['HBC_1_EGR1_PLA2G7', 'HBC_4_IL1B_MMP19', 'HBC_6_VIM_VMO1', 'HBC_13_S100A8','HBC_12_CXCL10','HBC_9_CXCL10_IFIT44','HBC_7_lncRNA',  'HBC_10_cycling', 'HBC_11', 'HBC_2_CXCL2_CXCL3', 'HBC_8_neuro', 'HBC_3','HBC_5']],
                "mMac_Mon_12subsets": [data_dir + "mMac_Mon_12subsets.h5mu", ["Mon_1_CD14_CLEC4D","Mon_2_CD14_CXCL1","mMon_4_MS1","fDC_CLECL10A_S100A8","fMon_3_CD16","mMac_1_C3_LTC4S","mMac_2_MMP7_CHI3L1","mMac_7_APOE","mMac_4_NUPR1","mMac_6_CX3CR1","mMac_5_CXCL10_IDO1","mMac_3_MT-CO1"]],
                "Mixed_HBC": [data_dir+"mixed_HBC.h5mu", None],
                "DCs": [data_dir+"/DCs_6subsets.h5mu", ['fDC2', 'fDC3', 'mDC2/4', 'mDC1', 'mDC_mature',"pDC"]],
                "T": [data_dir+"/T_cells.h5mu",  ["fT_8_CD8_naive","fT_3_CD4_naive","mT_1_CD4_eff","mT_7_CD8_eff_GZMB","mT_2_CD8_eff_GZMK","mT_4_CD8_eff_LAG3","mT_9_CD4_eff_CTLA4","mT_5_cycling","mT_15_CD4_CXCL13","mT_6_Treg","mT_12_MAIT","mT_10_gd_TRDV1_eff","T_11_gd_TRDV1_naive","mT_13_gd_DN_TRDV2-9","mT_14_gd_CD8_TRDV2-9"]],
                "NK": [data_dir+"/NK.h5mu", ['mNK_3_GNLY_KIR2DL4', 'mNK_2_CXCR4', 'fNK_4_TNF_KLRG1', 'NK_6_SELL_IL7R', 'NK_5_cycling','mNKT_7','fNKT_1_GZMB_GNLY', 'fNKT_8_PRDM16_KLRC4']],
                "B_plasma": [data_dir+"/B_plasma.h5mu", None]}

for lin in data_objects.keys():
    pdf_file = lin + "_umap_trio.pdf"
    if not(Path(pdf_file).exists()):
        print("**** "+ lin +" ****")
        func.plot_UMAP_MatFet_trio(lin_obj_path=data_objects[lin][0], title=lin)

os.chdir("/projects/Pregnancy_placenta/Results/22_Final_figures/Figure_1/UMAP_lineages/full_names/")
for lin in data_objects.keys():
    pdf_file = lin + "_umap_trio.pdf"
    if not(Path(pdf_file).exists()):
        print("**** "+ lin +" ****")
        if (lin=="DCs"):
            func.plot_UMAP_MatFet_trio(lin_obj_path=data_objects[lin][0], title=lin, cluster_order=data_objects[lin][1], dotsize=15)
        else:
            if (lin == "HBC_13subsets"):
                data = mu.read(data_objects[lin][0])
                clust =['HBC_1_EGR1_PLA2G7', 'HBC_2_CXCL2_CXCL3', 'HBC_3', 'HBC_4_IL1B_MMP19', 'HBC_5', 'HBC_6_VIM_VMO1', 'HBC_7_lncRNA', 'HBC_8_neuro', 'HBC_9_CXCL10_IFIT44','HBC_10_cycling', 'HBC_11', 'HBC_12_CXCL10', 'HBC_13_S100A8']
                func.plot_UMAP_MatFet_trio(lin_obj_path=data_objects[lin][0], title=lin, cluster_order=data_objects[lin][1], cluster_palette=dict(zip(clust ,data["GEX"].uns['sub_lin2_colors'])))
            else:
                func.plot_UMAP_MatFet_trio(lin_obj_path=data_objects[lin][0], title=lin, cluster_order=data_objects[lin][1])




#========================================================
# Fig1G) Mean frequencies per conditions
#========================================================
os.chdir("/projects/Pregnancy_placenta/Results/22_Final_figures/Figure_1/UMAP_lineages/")
data = mu.read(data_dir + "AllCells.h5mu")

# Frequencies per sample
total_per_samp = data["GEX"].obs["Sample"].value_counts().reset_index()
total_per_samp.rename(columns={"count":"total"},inplace=True)
counts_lin_samp = data["GEX"].obs[["Sample","Lineage"]].value_counts().reset_index()
counts_lin_samp = counts_lin_samp.merge(total_per_samp, on="Sample", how="left")
counts_lin_samp["freq"] = (counts_lin_samp["count"] * 100) / counts_lin_samp["total"]

# Mean frequencies per conditions
samp_cond = pd.read_csv(code_dir + "AuxFiles/sample_condition_list.csv")
condition_mean = pd.DataFrame()
for cond in samp_cond["Condition_group"]:
    curr_samp = list(samp_cond.loc[samp_cond["Condition_group"] == cond,"Sample"])
    curr_df = counts_lin_samp.loc[counts_lin_samp["Sample"].isin(curr_samp),:].copy()
    curr_df["Sample"] = curr_df["Sample"].cat.remove_unused_categories()
    curr_table = pd.pivot_table(curr_df[["Lineage","Sample","freq"]], index="Lineage", columns="Sample", values="freq", aggfunc=sum)
    curr_table["Mean"] = curr_table.mean(axis=1)
    condition_mean[cond] = curr_table["Mean"]

# Stacked bar plot of lineages per condition
cmap_lin = {"Myeloid":"#590902","HBC":"#BF3604","T":"#F27405","NK":"#F2AE2E","B_Plasma":"#F2668B","Trophoblasts": "#60B364","EVT":"#638c3e","dStromal":"#022873","fPerivascular":"#3f6cb4","fFibroblasts":"#92cdf0","Endothelial":"#6B0AC9","Other":"#000000","HSC":"grey"}
bar_order = ["HBC","Trophoblasts","Myeloid","fPerivascular","fFibroblasts","T","NK","EVT","B_Plasma","Endothelial","dStromal","HSC","Other"]
cond_order = []
condition_mean = condition_mean.T[bar_order]
condition_mean = condition_mean.loc[['Cont_labored', 'Cont_unlabored', 'Cont_PTB', 'PTB', 'DM1', 'FGR', 'PrePE', 'TermPE'],:]

fig, ax = plt.subplots()
ax = condition_mean.plot.barh(stacked=True, figsize=(6, 4), color=[cmap_lin[c] for c in condition_mean.columns], ax=ax)
handles, labels = ax.get_legend_handles_labels()
idx = [labels.index(l) for l in bar_order if l in labels]
ax.legend([handles[i] for i in idx], [labels[i] for i in idx], loc="upper left", bbox_to_anchor=(1.02, 1), borderaxespad=0, frameon=False, fontsize=10)
ax.grid(False)
#plt.show()
fig.savefig("Fig1D_stacked_lineage_barplot.pdf", bbox_inches="tight")

