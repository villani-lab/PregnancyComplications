import os
import sys
import pandas as pd
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt

pd.set_option('display.max_columns', 10)
pd.set_option("display.max_rows", 50)
pd.set_option('display.width', 400)
matplotlib.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']

code_dir="/projects/Pregnancy_placenta/Code_Publication/Python/"
sys.path.append(code_dir + "Functions/")
import FigFunctions as fig_func
import FC_Dotplot_Functions as fc_func

work_dir = "/projects/Pregnancy_placenta/Results/22_Final_figures/Figure_4/"
os.chdir(work_dir)

####################################################
# Fig4A - DE/ DA bar plots
####################################################
fig_func.doubleBarPlot_DA_DE(comp="PTB_vs_Cont_PTB",
                       pval_adj_thresh=0.15,
                       logFC_thresh=0.58,
                       pdf=True,
                       path = code_dir + "AuxFiles/",
                       onlySig=True,
                       figsize=(7,4),
                       DEGs_transparency_thresh=50,
                       mimDEGs_thresh=50)

####################################################
# Fig4B - PTB Pathway enrichment heatmap
####################################################
# ==> See R code


####################################################
# Fig4C - PTB gene dotplots
####################################################
os.chdir("/projects/Pregnancy_placenta/Results/22_Final_figures/Figure_4/")
base_path = "/projects/Pregnancy_placenta/Results/13_DE_PAL/"

# Immune cells
#====================
gene_groups = {"Inflammation": ["CXCL8","IL1B","NFKB1","OSM","RIPK2","ICAM1"],
               "Stress": ["HSP90AA1","HSPB1","HSPH1","DNAJB1", "HSPE1","GRPEL1", "LDHB"],
               "Angiogenesis": ["CYCS","ATF3","CCND1","APLNR","CPNE2","PROM1"],
               "T_function": ["IFNG","ICOS","CD83","NAMPT"],
               "Cell_migration":["BAG3","GPR183","PLAUR"],
               "ECM":["COL3A1","FN1"],
               "ISGs": ["HLA-DQA2", "HLA-DMA", "HLA-DMB", "ISG15", "TRIM22"]}

# Map celltypes -> Excel file path
file_map = {
    # HBC
    "HBC_4_IL1B_MMP19": os.path.join(base_path, "Myeloid/HBC/Level_2/DEGs_HBC_V5_Level2_PTB_vs_Cont_PTB.xlsx"),
    "HBC_6":             os.path.join(base_path, "Myeloid/HBC/Level_2/DEGs_HBC_V5_Level2_PTB_vs_Cont_PTB.xlsx"),
    "HBC_1":             os.path.join(base_path, "Myeloid/HBC/Level_2/DEGs_HBC_V5_Level2_PTB_vs_Cont_PTB.xlsx"),
    # PAMM
    "mMac_upper":        os.path.join(base_path,"Myeloid/PAMM/Level_1.5//DEGs_PAMM_V5_Level1.5_PTB_vs_Cont_PTB.xlsx"),
    "CD14_monocyte":      os.path.join(base_path,"Myeloid/Mon/Level_1.5//DEGs_Mon_V5_Level1.5_PTB_vs_Cont_PTB.xlsx"),
    "CD16_monocyte":      os.path.join(base_path,"Myeloid/Mon/Level_1.5//DEGs_Mon_V5_Level1.5_PTB_vs_Cont_PTB.xlsx"),
    # DCs
    "cDC": os.path.join(base_path, "Myeloid/DC/Level_1.5/DEGs_DC_V5_Level1.5_PTB_vs_Cont_PTB.xlsx"),
    # T
    "mCD4_eff": "/projects/Pregnancy_placenta/Results/13_DE_PAL/T/Level_2/DEGs_T_V6_Level2_PTB_vs_Cont_PTB.xlsx",
    "mCD8_eff_GZMK": "/projects/Pregnancy_placenta/Results/13_DE_PAL/T/Level_2/DEGs_T_V6_Level2_PTB_vs_Cont_PTB.xlsx",
    # Trophoblasts
    "VCT_1":             os.path.join(base_path, "Trophoblasts/Level_3/DEGs_Troph_V6_Level3_PTB_vs_Cont_PTB.xlsx"),
    "VCTfusion_1":       os.path.join(base_path, "Trophoblasts/Level_3/DEGs_Troph_V6_Level3_PTB_vs_Cont_PTB.xlsx"),
    # Endo
    "Endo_m":         os.path.join(base_path, "Endothelial/Level_2/DEGs_Endo_V4_Level2_PTB_vs_Cont_PTB.xlsx"),
    "Endo_f":         os.path.join(base_path, "Endothelial/Level_2/DEGs_Endo_V4_Level2_PTB_vs_Cont_PTB.xlsx")
}

import importlib
fc_func = importlib.reload(fc_func)
fc_func.plot_FC_dotplot(file_map, gene_groups, title="PTB vs PTB_cont", pdf_name="Fig4C_gene_dotplot.pdf")


## change name of genes "G:" and names of clusters