import os
import sys

import anndata
import pandas as pd
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
import muon as mu

pd.set_option('display.max_columns', 10)
pd.set_option("display.max_rows", 50)
pd.set_option('display.width', 400)
matplotlib.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']

code_dir="/projects/Pregnancy_placenta/Code_Publication/"
sys.path.append(code_dir)
import FigFunctions as fig_func

# default lineage palette you provided
_DEF_LINEAGES = ['Trophoblasts','EVT','fPV','fFB','dS','Endothelial','HBC','Myeloid','T','NK','B_Plasma','HSC','Other']
_DEF_COLORS = ["#60B364","#B2DF8A",  # trophoblasts
               "#08306b","#2b8cbe","#c6dbef","#8F6FE8",  # stromal/endo
               "#800026","#f03b20","#fd8d3c","#fed976","#fa9fb5",  # immune
               "grey","lightgrey"]
_DEF_PALETTE = dict(zip(_DEF_LINEAGES, _DEF_COLORS))

work_dir = "/projects/Pregnancy_placenta/Results/22_Final_figures/Figure_3/"
os.chdir(work_dir)

#==============================================================
# Running DA
#==============================================================
os.chdir("/projects/Pregnancy_placenta/Results/23_Supp_tables/DA/")
lineages=['Trophoblasts','EVT',
          'fPerivascular',"fFibroblasts", 'dStromal', 'Endothelial',
          'HBC', 'Myeloid','T','NK','B_Plasma',"HSC","Other"]
lin_colors = dict(zip(lineages, ["#60B364","#B2DF8A", # greens for trophoblats
                                 '#08306b','#2b8cbe','#c6dbef','#8F6FE8', # blues for stromal, purple for endo
                                 "#800026","#f03b20", "#fd8d3c","#fed976","#fa9fb5", # reds for immune
                                 "grey","lightgrey"]))
condition_colors = dict(zip(["Cont_labored","Cont_unlabored",'Cont_PTB',"DM1","FGR","PrePE","TermPE","PTB"],
                            ["#CCEA8D", "#A6BC09", "green", "#5DB3D8", "#FFC817", "#EF7E29", "#E53938", "#8F2F82"]))

# Making a myeloid object for all myeloid cells
HBC_data = mu.read("/projects/Pregnancy_placenta/Results/21_Final_objects/GEO_objects/HBC_13subsets.h5mu")
mon_mMac_data = mu.read("/projects/Pregnancy_placenta/Results/21_Final_objects/GEO_objects/mMac_Mon_12subsets.h5mu")
DCs = mu.read("/projects/Pregnancy_placenta/Results/21_Final_objects/GEO_objects/DCs_6subsets.h5mu")
for x in ['X_pca', 'X_pca_harmony', 'X_umap', 'pca_harmony_knn_distances', 'pca_harmony_knn_indices']:
    HBC_data["GEX"].obsm.pop(x)
    mon_mMac_data["GEX"].obsm.pop(x)
    DCs["GEX"].obsm.pop(x)
HBC_data["GEX"].varm.pop("de_final")
DCs["GEX"].varm.pop("de_final")
mon_mMac_data["GEX"].varm.pop("de_final")

# Cleaning obs to only leave these: ['Patient_ID', 'Batch', 'Condition', 'Lineage', 'Cell_subset']:
for x in [ 'Parent_lineage', 'Cell_subset_short', 'Cell_subset_num', 'Fetal/ Maternal', 'Maternal_age', 'GA_delivery', 'Fetal_sex', 'Labor', 'Birthweight_percentile', 'Condition_UnlaboredCont', 'Condition_LaboredCont', 'Condition_PTBCont', 'Condition_sPTB', 'Condition_PrePE', 'Condition_PE', 'Condition_FGR', 'Condition_DM1', 'n_counts', 'n_genes', 'percent_mito']:
    mon_mMac_data["GEX"].obs.pop(x)
    HBC_data["GEX"].obs.pop(x)
    DCs["GEX"].obs.pop(x)

import anndata as an
myeloid_data = an.concat([HBC_data["GEX"], mon_mMac_data["GEX"], DCs["GEX"]])
myeloid_mu = mu.MuData({"GEX": myeloid_data, "ADT": myeloid_data}) # Just a placeholder for ADT object
myeloid_mu.write("AllMyeloid.h5mu")

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
lin_dict = {"B_Plasma":"B_plasma.h5mu",
            "dS":"dS.h5mu",
            "Endo":"Endo.h5mu",
            "EVT":"EVT.h5mu",
            "fFB":"fFB.h5mu",
            "fPV":"fPV.h5mu",
            "AllMyeloid":"AllMyeloid.h5mu",
            "NK":"NK.h5mu",
            "T":"T_cells.h5mu",
            "Trophoblasts":"Trophoblasts.h5mu"}


for lin in lin_dict.keys():
    print(lin)
    if lin == "AllMyeloid":
        data = mu.read(lin_dict[lin])
    else:
        data = mu.read("/projects/Pregnancy_placenta/Results/21_Final_objects/GEO_objects/"+lin_dict[lin])
    fig_func.DA_condition_groups(data,
                                 cluster_col='Cell_subset',
                                 comparisons = comparisons,
                                 samp_condition_file= code_dir+"AuxFiles/sample_condition_list.csv",
                                 lineage= lin)
lin="AllCells"
data = mu.read("/projects/Pregnancy_placenta/Results/21_Final_objects/AllCells.h5mu")
fig_func.DA_condition_groups(data,
                             cluster_col='Lineage',
                             comparisons = comparisons,
                             samp_condition_file= code_dir+"AuxFiles/sample_condition_list.csv",
                             lineage= lin)

# Summarizing all DA - stats and frequencies
da_stats = pd.DataFrame()
da_freq = pd.DataFrame()
for lin in ["AllCells", 'B_Plasma', 'dS', 'Endo', 'EVT', 'fFB', 'fPV', 'AllMyeloid', 'NK', 'T', 'Trophoblasts']:
    print(lin)
    df = pd.read_csv(lin+"_comparisons_stats.csv")
    df["Lineage"]=lin
    da_stats = pd.concat([da_stats, df])

    freq = pd.read_csv(lin+"_frequencies_per_pt.csv")
    freq = freq.T
    freq.columns = freq.loc["Patient_ID",:]
    freq.drop("Patient_ID", inplace=True)
    freq["Lineage"]=lin
    da_freq = pd.concat([da_freq, freq])

# merging these to one table
stats_cols = ["mean_case", "mean_cont", "log(FC)", "P-val", "qval"]
pivot_df = da_stats.pivot(index="Cluster", columns="Comparison", values=stats_cols)
pivot_df = pivot_df.swaplevel(0, 1, axis=1) # switch column levels: stat/comparison -> comparison/stat
pivot_df = pivot_df.sort_index(axis=1, level=0) # sort by comparison
pivot_df.columns = [f"{comparison} : {stat}" for comparison, stat in pivot_df.columns]  # Flattening the multi-layer column names
final = da_freq.merge(pivot_df, left_index=True, right_index=True, how="outer")

final.to_excel("Table_S7.xlsx") # Table S7
da_freq.to_csv("all_frequencies.csv") # All DAs
da_stats.to_csv("all_frequencies.csv") # All stats


#======================================================================================================================
# Running DE
# In this section we will create pseudo-bulk files that will be used by DEseq which will be run in R
#======================================================================================================================
fig_func = importlib.reload(fig_func)
os.chdir("/projects/Pregnancy_placenta/Results/23_Supp_tables/DEG_comparisons/")

for lin in lin_dict.keys():
    print(lin)
    data = mu.read("/projects/Pregnancy_placenta/Results/21_Final_objects/"+lin_dict[lin])
    fig_func.create_PB_files(data = data["GEX"],
                             sample_label="Sample",
                             cluster_label='sub_lin2',
                             condition_label="Condition",
                             other_meta=["Batch"],
                             title=lin,
                             doublet_clust_to_rm=[],
                             doublet_column="")




#==============================================================
# Figure 3A - DA of all condition against unlabored controls
#==============================================================
comparisons = ["PTB_vs_Cont_unlabored",
               "PrePE_vs_Cont_unlabored",
               "TermPE_vs_Cont_unlabored",
               "DM1_vs_Cont_unlabored",
               "FGR_vs_Cont_unlabored"]

import importlib
fig_func = importlib.reload(fig_func)

# main figure
fig_func.multiDA_BarGrid(
    comparisons=comparisons,
    pdf_name="3A_DA_Unlabored_onlySig.pdf",
    path_to_stats="/projects/Pregnancy_placenta/Results/22_Final_figures/DE_DA_results/",
    da_pval_thresh=0.05,  # This threshold is used to have a transparent fill for populations with pval< <thresh>
    onlyMinDAPop=True,   # include only populations with a (p_val<0.05) AND (FC>minDA_absFC_thresh) in ANY comparison
    minDA_absFC_thresh=0, # threshold applied to abs(log(FC))
    share_y=False,
    xticks_fontsize=10,
    width=8,
    global_order=['SCT', 'VCT_3', 'VCTfusion_2', 'EVT_1', 'EVT_2',
        'HBC_1', 'HBC_4', 'HBC_7', 'HBC_11', 'HBC_13',
        'fDC2', 'mDC2/4','Mon_1_CD14', 'mMon_4_MS1',
        'mMac_1', 'mMac_2', 'mMac_3', 'mMac_4', 'mMac_6',
        'mT_1_CD4_eff', 'fNK_4', 'mNK_2', 'mB_2',
        'fPV_2_1', 'fPV_2_2', 'fPV_3_2', 'fPV_4', 'fFB_1',
        'dS_1_2', 'dS_3', 'dS_4_1', 'dS_4_2',
        'dS_lncRNA', 'fEndo_1', 'fEndo_2', 'fEndo_cycling', 'mEndo_1', 'mEndo_2', 'mEndo_lymph'])

#==============================================================
# Figure 3B - DA in myeloid cells
#==============================================================
# ===> see R code



#==============================================================
# Figure 3C - DE of all condition against unlabored controls
#==============================================================
comparisons = ["PTB_vs_Cont_unlabored",
               "PrePE_vs_Cont_unlabored",
               "TermPE_vs_Cont_unlabored",
               "DM1_vs_Cont_unlabored",
               "FGR_vs_Cont_unlabored"]
fig_func.multiDE_BarGrid(comparisons=comparisons,
                         pval_adj_thresh=0.15,
                         logFC_thresh=0.58,
                         pdf_name="3C_DE_unlabored_min50DEGs.pdf",
                         path_to_stats=code_dir+"AuxFiles/",
                         xticks_fontsize=10,
                         mimDEGs_thresh=50, # include only populations with a minimum <mimDEGs_thresh> DEGs
                         onlyMinDEGsPop=True,
                         DEGs_transparency_thresh=0, # A threshshold below which the bar will be partially transparent
                         width=8,
                         global_order=['SCT', 'VCT_1', 'VCT_2', 'VCT_3', 'VCT_cycling_1', 'VCT_cycling_2', 'VCTccc', 'VCTfusion_1', 'VCTfusion_2', 'HBC_1', 'HBC_2', 'HBC_3', 'HBC_4', 'HBC_5', 'HBC_6', 'HBC_7', 'HBC_9','HBC_11', 'HBC_13','Mon_1_CD14', 'Mon_2_CD14', 'fMon_3_CD16', 'fPV_2_1', 'fFB_1', 'dS_1_1', 'dS_3', 'fEndo_1', 'fEndo_2', 'fEndo_3'])


#==============================================================
# Figure 3D - Pathway enrichment across conditions
#==============================================================
# ===> see R code



