import seaborn as sns
import os
import anndata as an
import sys
import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
import muon as mu
from muon import MuData
import matplotlib
code_dir="/projects/Pregnancy_placenta/Code_Publication/"
sys.path.append(code_dir)
import FigFunctions as func

sc.set_figure_params(vector_friendly=True, dpi_save=300)
matplotlib.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']

os.chdir("/projects/Pregnancy_placenta/Results/22_Final_figures/Figure_2/")


############################################################
# Fig 2B) EVT dotplot
###########################################################
EVT_subcluster = mu.read("/projects/Pregnancy_placenta/Results/11_Subclustering_B1-15/EVT/RD1_EVT_webbrowser/RD1_EVT_V6.h5mu")

# gene expression
EVT_adata = EVT_subcluster['GEX'].copy()

# proteiin
EVT_pdata = EVT_subcluster['ADT'].copy()
EVT_pdata.obs["sub_lin"] = ("EVT" + EVT_pdata.obs["final_clusters"].astype(str)).astype("category")

strip_categories = pd.DataFrame({
    "markers": [
        "G:FSTL3", "G:FN1", "G:HLA-G",
        "P:CD73","P:CD155",
        "G:TAC3", "G:CRABP2",
        #"G:TAC3", "G:HTRA4",
        "G:ENG", "G:INHBA",
         "G:KCNQ1OT1", "G:HERC1","G:CORO6",
        "G:RPS12","G:RPL39",
        "G:COL5A1", "G:KRT17",
        "P:CD49b",'P:CD49f',"P:CD146",
        "G:PAPPA","G:NOTUM",
        "G:IGFBP3", "G:TFPI2", "G:GPX3"
    ],
    "types": [
        "RNA", "RNA", "RNA",
        "Protein","Protein",
        "RNA", "RNA",
        #"RNA", "RNA",
         "RNA", "RNA",
        "RNA", "RNA",
        "RNA", "RNA","RNA",
        "RNA", "RNA",
        "Protein", "Protein","Protein",
        "RNA", "RNA",
        "RNA", "RNA", "RNA"
    ],
    "blocks": [
        "EVT markers", "EVT markers", "EVT markers",
        "Protein markers","Protein markers",
        "EVT1 marker", "EVT1 marker",
        #"EVT1 and 3", "EVT1 and 3",
        "EVT3", "EVT3",
         "EVT4", "EVT4",
        "EVT4", "EVT4","EVT4",
        "EVT5", "EVT5",
        "EVT5_Protein", "EVT5_Protein", "EVT5_Protein",
        "EVT5", "EVT5",
        "EVT6", "EVT6", "EVT6"
    ]
})


import importlib
func = importlib.reload(func)

func.dot_plot_with_umap_stack(
    EVT_adata,
    prot_anndata= EVT_pdata,
    groupby="sub_lin",
    strip_categories=strip_categories,
    shown_categories=["EVT1", "EVT2", "EVT3", "EVT4", "EVT5", "EVT6"],
    legend_fontsize=16,
    legend_fontoutline=0,
    xtick_fontsize=16,
    category_size=16,
    legend_length=23,
    umap_aspect=1.3,
    legend_dotsize=80,
   #standard_scale="var",
    figsize=(6, 3),
    cbar_orientation="vertical",
    cbar_loc=[1.05, 0.3, 0.015, 0.61],
    legend_loc=[.95, 0.05, 0.2, 0.1],
    save_fig="EVT_UMAP_dotplot.pdf",
    dot_size_exponent=1)

