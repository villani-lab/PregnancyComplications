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

code_path="/projects/Pregnancy_placenta/Code/Publication_code"
sys.path.append(code_path)
import FigFunctions as fig_func
import FC_Dotplot_Functions as fc_func

work_dir = "/projects/Pregnancy_placenta/Results/22_Final_figures/SuppFigure_6/"
os.chdir(work_dir)

#==================================================================
# SupFig 6 A-C: Additional comparisons for DA and DE
#==================================================================
fig_func.doubleBarPlot_DA_DE(comp="DM1_vs_Cont_unlabored",
                       pval_adj_thresh=0.15,
                       logFC_thresh=0.58,
                       pdf=True,
                       path = code_path + "/AuxFiles/",
                       onlySig=True,
                       figsize=(7,4),
                       DEGs_transparency_thresh=50,
                       mimDEGs_thresh=50)

fig_func.doubleBarPlot_DA_DE(comp="DM1_vs_Cont_labored",
                       pval_adj_thresh=0.15,
                       logFC_thresh=0.58,
                       pdf=True,
                       path = code_path + "/AuxFiles/",
                       onlySig=True,
                       figsize=(7,4),
                       DEGs_transparency_thresh=50,
                       mimDEGs_thresh=50)

fig_func.doubleBarPlot_DA_DE(comp="FGR_vs_Cont_labored",
                       pval_adj_thresh=0.15,
                       logFC_thresh=0.58,
                       pdf=True,
                       path = code_path + "/AuxFiles/",
                       onlySig=True,
                       figsize=(7,4),
                       DEGs_transparency_thresh=50,
                       mimDEGs_thresh=50)

fig_func.doubleBarPlot_DA_DE(comp="FGR_vs_Cont_unlabored",
                       pval_adj_thresh=0.15,
                       logFC_thresh=0.58,
                       pdf=True,
                       path = code_path + "/AuxFiles/",
                       onlySig=True,
                       figsize=(7,4),
                       DEGs_transparency_thresh=50,
                       mimDEGs_thresh=50)
