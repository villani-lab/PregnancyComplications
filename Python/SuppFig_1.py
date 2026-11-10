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
code_dir="/projects/Pregnancy_placenta/Code_Publication/Python/"
sys.path.append(code_dir + "Functions/")
import FigFunctions as func
#import Marker_dotplot_fig_RN as func_dotplot
import Func_Marker_dotplot_fig as func_dotplot
os.chdir("/projects/Pregnancy_placenta/Results/22_Final_figures/SuppFigure_1/")


#=================================================================
# Supp Fig 1B) UMAP showing fetal-maternal demultiplexing results
#=================================================================
import importlib
func = importlib.reload(func)

_DEF_LINEAGES = ['Trophoblasts','EVT','HBC','mMac',"Monocytes","DC",'T','NK','B_Plasma','dStromal','fPerivascular','fFibroblasts', 'Endothelial','HSC','Other']
_DEF_COLORS = ["#60B364","#B2DF8A",  # trophoblasts
               "#800026","#f03b20","#fa9fb5","#E2C087","#fd8d3c","#fed976","#E80C7A",  # immune
               "#1f2b66", "#5384d4", "#5ab5ed", "#8F6FE8",  # stromal/endo
               "grey","lightgrey"]
_DEF_PALETTE = dict(zip(_DEF_LINEAGES, _DEF_COLORS))

data = mu.read(data_dir+"AllCells.h5mu")


func = importlib.reload(func)
func.plot_UMAP_MatFet_trio(
    lin_obj_path=data,
    title="AllCells",
    cluster_col="Lineage",
    cluster_order=_DEF_LINEAGES,
    xlabel="% cells",
    xlabel_fontsize=12,
    legend_fontsize=10,
    legend_title_fontsize=10,
    figsize=(10, 3),
    cluster_palette=_DEF_PALETTE,
    add_bar_count=True,
    count_ticks=[0,100000,200000,300000])


# Supp Fig 1B) Dotplot All cells
#================================
os.chdir("/projects/Pregnancy_placenta/Results/22_Final_figures/SuppFigure_1/")
all_data = mu.read(data_dir+"/AllCells.h5mu")
all_gex = all_data['GEX']  # gene expression
all_ADT = all_data['ADT'] # protein
strip_categories = func_dotplot.read_markers_dotplot(file_path= code_dir + "AuxFiles/marker genes_v3.xlsx",
                          lineage="All")

clust_order = ['Trophoblasts',"EVT","Endothelial","dStromal","fPerivascular","fFibroblasts","HBC","Myeloid","T","NK","B_Plasma",'HSC','Other']
all_gex.obs["Lineage"] = pd.Categorical(all_gex.obs["Lineage"], ordered=True, categories=clust_order)
all_ADT.obs["Lineage"] = all_gex.obs.loc[all_ADT.obs_names, "Lineage"]

genes = list(strip_categories.loc[strip_categories["types"]=="RNA","markers"])
proteins = list(strip_categories.loc[strip_categories["types"]=="Protein","markers"])
all_gex = all_gex[:,genes]
all_ADT = all_ADT[:,proteins]
cats = list(all_gex.obs["Lineage"].cat.categories)
all_gex.uns["Lineage_colors"] = [_DEF_PALETTE[c] for c in cats]

func_dotplot.dot_plot_with_umap_stack(
    rna_anndata=all_gex,
    prot_anndata= all_ADT,
    groupby="Lineage",
    strip_categories=strip_categories,
    standard_scale="var",
    legend_fontsize=16,
    legend_fontoutline=0,
    xtick_fontsize=12,
    category_size=12,
    umap_title_fontsize=18,
    legend_length=20,
    umap_title = "All cells",
    umap_aspect=1,
    legend_dotsize=80,
    figsize=(20, 5),
    cbar_orientation="vertical",
    cbar_loc=[1.05, 0.3, 0.015, 0.61],
    legend_loc=[.95, 0.05, 0.2, 0.1],
    save_fig="supp1_All_cells.pdf",
    dpi = 600,
umap_legend_loc=None)


import importlib
func_dotplot = importlib.reload(func_dotplot)

pb_all = pd.read_csv("All_cells_prot_genes_exp.csv")

func_dotplot.dot_plot_with_umap_stack_plot_only(
    pb_all,
    rna_anndata=all_gex,
    prot_anndata= all_ADT,
    groupby="Lineage",
    strip_categories=strip_categories,
    legend_fontsize=16,
    legend_fontoutline=0,
    xtick_fontsize=12,
    category_size=12,
    umap_title_fontsize=18,
    legend_length=20,
    umap_title = "All cells",
    umap_aspect=1,
    legend_dotsize=80,
    figsize=(20, 5),
    cbar_orientation="vertical",
    cbar_loc=[1.05, 0.3, 0.015, 0.61],
    legend_loc=[.95, 0.05, 0.2, 0.1],
    save_fig="supp1_All_cells.pdf",
    dpi = 600,
    umap_legend_loc="None")


#========================================================================
# Supp Fig 1G, I, L, O, G - demultiplexing results on UMAP and bar pltos
#========================================================================
# see code in Figure_1.py under "# Fig1C-F) UMAPs of different lineages with demultiplexing"


#=======================================================================
# Supp Fig 1E) Dotplot T cells
#=======================================================================
os.chdir("/projects/Pregnancy_placenta/Results/22_Final_figures/SuppFigure_1/")

T_cell_data = mu.read(data_dir+"/T_cells.h5mu")
T_gex = T_cell_data['GEX'].copy()  # gene expression
T_ADT = T_cell_data['ADT'].copy() # protein
strip_categories = func_dotplot.read_markers_dotplot(file_path= code_dir + "AuxFiles/marker genes_v3.xlsx", lineage="T")

# ordering clusters
clust_order = ["fT_8_CD8_naive","fT_3_CD4_naive","mT_1_CD4_eff","mT_7_CD8_eff_GZMB","mT_2_CD8_eff_GZMK","mT_4_CD8_eff_LAG3","mT_9_CD4_eff_CTLA4","mT_5_cycling","mT_15_CD4_CXCL13","mT_6_Treg","mT_12_MAIT","mT_10_gd_TRDV1_eff","T_11_gd_TRDV1_naive","mT_13_gd_DN_TRDV2-9","mT_14_gd_CD8_TRDV2-9"]
T_gex.obs["sub_lin2"] = pd.Categorical(T_gex.obs["sub_lin2"], ordered=True, categories=clust_order)
T_ADT.obs["sub_lin2"].replace({"mT_9_CD4_act_CTLA4":"mT_9_CD4_eff_CTLA4", "mT_4_CD8_act":"mT_4_CD8_eff_LAG3"}, inplace=True)
T_ADT.obs["sub_lin2"] = pd.Categorical(T_ADT.obs["sub_lin2"], ordered=True, categories=clust_order)

func_dotplot.dot_plot_with_umap_stack(
    T_gex,
    prot_anndata= T_ADT,
    groupby="sub_lin2",
    strip_categories=strip_categories,
    standard_scale="var",
    legend_fontsize=16,
    legend_fontoutline=0,
    xtick_fontsize=12,
    category_size=12,
    umap_title_fontsize=18,
    legend_length=20,
    umap_title = "Tcells",
    umap_aspect=1,
    legend_dotsize=80,
    figsize=(20, 5),
    cbar_orientation="vertical",
    cbar_loc=[1.05, 0.3, 0.015, 0.61],
    legend_loc=[.95, 0.05, 0.2, 0.1],
    save_fig="supp1_Tcells.pdf",
    dpi = 600)




#=======================================================================
# Supp Fig 1F) Dotplot NK cells
#=======================================================================
NK_cell_data = mu.read(data_dir+"/NK.h5mu")
NK_gex = NK_cell_data['GEX'].copy() # gene expression
NK_ADT = NK_cell_data['ADT'].copy() # protein
strip_categories = func_dotplot.read_markers_dotplot(file_path= code_dir + "AuxFiles/marker genes_v3.xlsx",
                          lineage="NK")

clust_order = ['mNK_3_GNLY_KIR2DL4', 'mNK_2_CXCR4', 'fNK_4_TNF_KLRG1', 'NK_6_SELL_IL7R', 'NK_5_cycling','mNKT_7','fNKT_1_GZMB_GNLY', 'fNKT_8_PRDM16_KLRC4']
NK_gex.obs["sub_lin2"] = pd.Categorical(NK_gex.obs["sub_lin2"], ordered=True, categories=clust_order)
NK_ADT.obs["sub_lin2"] = pd.Categorical(NK_ADT.obs["sub_lin2"], ordered=True, categories=clust_order)

func_dotplot.dot_plot_with_umap_stack(
    NK_gex,
    prot_anndata= NK_ADT,
    groupby="sub_lin2",
    strip_categories=strip_categories,
    standard_scale="var",
    legend_fontsize=16,
    legend_fontoutline=0,
    xtick_fontsize=12,
    category_size=12,
    umap_title_fontsize=18,
    legend_length=20,
    umap_title = "NK cells",
    umap_aspect=1,
    legend_dotsize=100,
    figsize=(14, 5),
    cbar_orientation="vertical",
    cbar_loc=[1.05, 0.3, 0.015, 0.61],
    legend_loc=[.95, 0.05, 0.2, 0.1],
    save_fig="supp1_NK_cells.pdf",
    dpi = 600)



#=======================================================================
# Supp Fig 1H) Dotplot T cells
#=======================================================================
DC_cell_data = mu.read(data_dir+"/DCs_6subsets.h5mu")
DC_gex = DC_cell_data['GEX'].copy() # gene expression
DC_ADT = DC_cell_data['ADT'].copy() # protein
strip_categories = func_dotplot.read_markers_dotplot(file_path= code_dir + "AuxFiles/marker genes_v3.xlsx",
                          lineage="DCs")

clust_order = ['fDC2', 'fDC3', 'mDC2/4', 'mDC1', 'mDC_mature',"pDC"]
DC_gex.obs["short_sublin2"] = pd.Categorical(DC_gex.obs["short_sublin2"], ordered=True, categories=clust_order)
DC_ADT.obs["short_sublin2"] = pd.Categorical(DC_ADT.obs["short_sublin2"], ordered=True, categories=clust_order)

# reset GEX colors (now it is red-grey)
DC_gex.uns.pop('sub_lin2_colors')
sc.pl.umap(DC_gex, color="short_sublin2")

fig, ax, pb_all = func_dotplot.dot_plot_with_umap_stack_efficient(
    DC_gex,
    prot_anndata= DC_ADT,
    groupby="short_sublin2",
    strip_categories=strip_categories,
    standard_scale="var",
    legend_fontsize=16,
    legend_fontoutline=0,
    xtick_fontsize=12,
    category_size=12,
    umap_title_fontsize=18,
    legend_length=20,
    umap_title = "DC cells",
    umap_aspect=1,
    umap_dotsize=15,
    legend_dotsize=100,
    figsize=(14, 5),
    cbar_orientation="vertical",
    cbar_loc=[1.05, 0.3, 0.015, 0.61],
    legend_loc=[.95, 0.05, 0.2, 0.1],
    save_fig="supp1_DC_cells.pdf",
    dpi = 600)



#=======================================================================
# Supp Fig 1J) Dotplot mMacs
#=======================================================================
mMac_data = mu.read(data_dir+"/mMac_Mon_12subsets.h5mu")
mMac_gex = mMac_data['GEX'].copy() # gene expression
mMac_ADT = mMac_data['ADT'].copy() # protein
strip_categories = func_dotplot.read_markers_dotplot(file_path= code_dir + "AuxFiles/marker genes_v3.xlsx",
                          lineage="Mon_mMac")

clust_order = ["Mon_1_CD14_CLEC4D","Mon_2_CD14_CXCL1","mMon_4_MS1","fDC_CLECL10A_S100A8","fMon_3_CD16","mMac_1_C3_LTC4S","mMac_2_MMP7_CHI3L1","mMac_7_APOE","mMac_4_NUPR1","mMac_6_CX3CR1","mMac_5_CXCL10_IDO1","mMac_3_MT-CO1"]
mMac_gex.obs["sub_lin2"] = pd.Categorical(mMac_gex.obs["sub_lin2"], ordered=True, categories=clust_order)
mMac_ADT.obs["sub_lin2"] = mMac_gex.obs.loc[mMac_ADT.obs_names, "sub_lin2"]

# reset GEX colors (now it is red-grey)
mMac_gex.uns.pop('sub_lin2_colors')
sc.pl.umap(mMac_gex, color="sub_lin2")
import importlib
func_dotplot = importlib.reload(func_dotplot)

func_dotplot.dot_plot_with_umap_stack_efficient(
    mMac_gex,
    prot_anndata= mMac_ADT,
    groupby="sub_lin2",
    strip_categories=strip_categories,
    standard_scale="var",
    legend_fontsize=16,
    legend_fontoutline=0,
    xtick_fontsize=12,
    category_size=12,
    umap_title_fontsize=18,
    legend_length=30,
    umap_title = "mMac mon cells",
    umap_aspect=1.5,
    umap_dotsize=15,
    legend_dotsize=100,
    figsize=(25, 5),
    cbar_orientation="vertical",
    cbar_loc=[1.05, 0.3, 0.015, 0.61],
    legend_loc=[.95, 0.05, 0.2, 0.1],
    save_fig="supp1_mMac_Mon_cells.pdf",
    dpi = 600)



# Supp Fig 1K) MS1 100 top gene score
#=======================================================================
ms1_sig = list(pd.read_csv("MS1_signature_top100.txt", header=None)[0])
ms1_sig = ["G:"+x for x in ms1_sig]
func.signature_score_per_cell(data = mdata["GEX"], gene_set=ms1_sig, score_title="MS1_100")

clust_order = ['mMon_4_MS1','Mon_1_CD14_CLEC4D', 'Mon_2_CD14_CXCL1','fMon_3_CD16', 'fDC_CLECL10A_S100A8', 'mMac_1_C3_LTC4S', 'mMac_2_MMP7_CHI3L1', 'mMac_3_MT-CO1', 'mMac_4_NUPR1', 'mMac_5_CXCL10_IDO1', 'mMac_6_CX3CR1', 'mMac_7_APOE']
# Add boxplot of mean MS score per patient (y axis) for every cluster (x axis)
score= "MS1_100"
with PdfPages("MS1_top100_marker_genes_score.pdf") as pdf:
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(6, 3))
    obs_gex = mdata["GEX"].obs
    sc.pl.umap(adata=mdata["GEX"], color=score, cmap="bwr", vcenter=0, size=5, ax=axes[0], show=False)

    mean_score = (mdata["GEX"].obs.groupby(["Cell_subset"], observed=True)[score].mean().reset_index())
    sns.barplot(data=mean_score, x="Cell_subset", y=score, order=clust_order, ax=axes[1], color="lightgrey")
    axes[1].set_ylabel(f"Mean {score} per cell subset")
    axes[1].grid(False)
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    pdf.savefig(bbox_inches="tight", dpi=300)


# Supp Fig 1L) Unified UMAP of HBC, DCs monocytes and mMacs
#=======================================================================

data = mu.read("/projects/Pregnancy_placenta/Results/11_Subclustering_B1-15/Myeloid/Myeloid_V5.h5mu")

data_hbc = mu.read("/projects/Pregnancy_placenta/Results/21_Final_objects/GEO_objects/HBC_13subsets.h5mu")
data_mMac = mu.read("/projects/Pregnancy_placenta/Results/21_Final_objects/GEO_objects/mMac_Mon_12subsets.h5mu")
data_dc = mu.read("/projects/Pregnancy_placenta/Results/21_Final_objects/GEO_objects/DCs_6subsets.h5mu")

data.shape[0]
data_dc.shape[0] + data_mMac.shape[0] + data_hbc.shape[0]

data["GEX"].obs[['sub_lin', 'sub_lin2']]
len(set(data["GEX"].obs['sub_lin2']))
sc.pl.umap(data["GEX"], color = "sub_lin")

# ---- UMAP colored by subset, with colors grouped by lineage and legends split left/right ----
import numpy as np
from matplotlib.lines import Line2D
adata = data["GEX"]
col = "Cell_subset"

# Transfer annotations from "Cell_subset" in the final objects to the unified adata, matching on cell barcode (obs_names)
lineage_objs = {"mMac": data_mMac, "DC": data_dc, "HBC": data_hbc}
subset_map = {}
lineage_subsets = {}
for lin, obj in lineage_objs.items():
    s = obj["GEX"].obs["Cell_subset"].astype(str)
    assert s.index.is_unique, f"Duplicate barcodes in {lin} object"
    subset_map.update(s.to_dict())
    lineage_subsets[lin] = list(dict.fromkeys(s))  # unique subsets, in order of appearance

adata.obs[col] = adata.obs_names.map(subset_map)
print("Cells without a Cell_subset annotation:", adata.obs[col].isna().sum(), "of", adata.n_obs)
adata.obs[col] = adata.obs[col].fillna("Unassigned")

monocytes = ['Mon_1_CD14_CLEC4D','Mon_2_CD14_CXCL1',  'fMon_3_CD16', 'mMon_4_MS1']
# Legend order: by the serial number after the first underscore (HBC_2 before HBC_10); names without a serial number go last
import re
def serial_key(name):
    m = re.search(r"_(\d+)(?:_|$)", name)
    return (0, int(m.group(1)), name) if m else (1, 0, name)

hbc_list = sorted(lineage_subsets["HBC"], key=serial_key)
dc_list = sorted(lineage_subsets["DC"], key=str.lower)
mon_list = list(monocytes)  # in the order listed above; monocytes are split out of the mMac/Mon group
mmac_list = sorted([s for s in lineage_subsets["mMac"] if s not in monocytes], key=serial_key)
print(hbc_list, mmac_list, sep="\n")
print(len(hbc_list), len(dc_list), len(mmac_list), len(mon_list))  # expect 13, 6, 8, 4
assert set(monocytes) <= set(lineage_subsets["mMac"]), "Monocyte names not found in the mMac object"

# One hue family per lineage, shades sampled from a colormap; HBC (right blob) cool, mMac/Mon (left blob) warm, DC blues/purples
def shades(cmap, n, lo=0.15, hi=0.95):
    return [matplotlib.colors.to_hex(plt.get_cmap(cmap)(x)) for x in np.linspace(lo, hi, n)]

palette = {}
palette.update(dict(zip(mmac_list, shades("Spectral", len(mmac_list), 0.0, 0.45))))   # reds/oranges/yellows
palette.update(dict(zip(mon_list, shades("copper", len(mon_list), 0.5, 1.0))))        # lighter browns (avoid the near-black start of viridis in HBC)
palette.update(dict(zip(dc_list, shades("cool", len(dc_list), 0.0, 0.8))))            # cyan -> magenta
palette.update(dict(zip(hbc_list, shades("viridis", len(hbc_list), 0.0, 0.9))))       # purple -> teal -> green -> yellow
palette["Unassigned"] = "#BBBBBB"

adata.obs[col] = adata.obs[col].astype(str).astype("category")
adata.uns[f"{col}_colors"] = [palette[c] for c in adata.obs[col].cat.categories]

fig, ax = plt.subplots(figsize=(6, 6))
sc.pl.umap(adata, color=col, ax=ax, legend_loc=None, frameon=False, size=1, title="", show=False)

def make_legend(names, title, loc, anchor):
    handles = [Line2D([], [], marker="o", ls="", ms=7, color=palette[n], label=n) for n in names]
    leg = ax.legend(handles=handles, title=title, loc=loc, bbox_to_anchor=anchor, frameon=False,
                    fontsize=8, title_fontsize=9, ncol=1, handletextpad=0.2, labelspacing=0.4)
    ax.add_artist(leg)
    return leg

# Left blob (HBC) legend to the left of the plot, right blob (mMac/Mon + DCs) legends to the right
legends = [
    make_legend(hbc_list, "HBC", "upper right", (-0.02, 1.0)),
    make_legend(mmac_list, "mMac", "upper left", (1.02, 1.0)),
    make_legend(mon_list, "Monocytes", "center left", (1.02, 0.35)),
    make_legend(dc_list, "DC", "lower left", (1.02, 0.0)),
]
# Explicitly include the legends in the saved bounding box so they are not clipped
fig.savefig("UMAP_sub_lin_split_legend.pdf", bbox_inches="tight", bbox_extra_artists=legends, dpi=300)




#=======================================================================
# Supp Fig 1K) Dotplot HBCs
#=======================================================================
HBC_data = mu.read(data_dir+"/HBC_13subsets.h5mu")
HBC_gex = HBC_data['GEX'].copy() # gene expression
HBC_ADT = HBC_data['ADT'].copy() # protein
strip_categories = func_dotplot.read_markers_dotplot(file_path= code_dir + "AuxFiles/marker genes_v3.xlsx",  lineage="HBC")

clust_order = ['HBC_1_EGR1_PLA2G7', 'HBC_4_IL1B_MMP19', 'HBC_6_VIM_VMO1', 'HBC_13_S100A8','HBC_12_CXCL10','HBC_9_CXCL10_IFIT44','HBC_7_lncRNA',  'HBC_10_cycling', 'HBC_11', 'HBC_2_CXCL2_CXCL3', 'HBC_8_neuro', 'HBC_3','HBC_5']
HBC_gex.obs["sub_lin2"] = pd.Categorical(HBC_gex.obs["sub_lin2"], ordered=True, categories=clust_order)
HBC_ADT.obs["sub_lin2"] = pd.Categorical(HBC_ADT.obs["sub_lin2"], ordered=True, categories=clust_order)

# reset GEX colors (now it is red-grey)
#HBC_gex.uns.pop('sub_lin2_colors')
sc.pl.umap(HBC_gex, color="sub_lin2")
import importlib
func_dotplot = importlib.reload(func_dotplot)

func_dotplot.dot_plot_with_umap_stack_efficient(
    HBC_gex,
    prot_anndata= HBC_ADT,
    groupby="sub_lin2",
    strip_categories=strip_categories,
    standard_scale="var",
    legend_fontsize=16,
    legend_fontoutline=0,
    xtick_fontsize=12,
    category_size=12,
    umap_title_fontsize=18,
    legend_length=50,
    umap_title = "HBCs",
    umap_aspect=0.8,
    umap_dotsize=15,
    legend_dotsize=100,
    figsize=(20, 5),
    cbar_orientation="vertical",
    cbar_loc=[1.05, 0.3, 0.015, 0.61],
    legend_loc=[.95, 0.05, 0.2, 0.1],
    save_fig="supp1_HBC_cells.pdf",
    dpi = 600)



sc.pl.umap(HBC_gex, color="short_sublin2", show=False)
plt.tight_layout()
plt.show()



#==================================================================
# Supp fig 1 M) DEGs HBC and mMacs
#==================================================================
os.chdir("/projects/Pregnancy_placenta/Results/22_Final_figures/SuppFigure_1/")
all_data = mu.read("/projects/Pregnancy_placenta/Results/21_Final_objects/AllCells.h5mu")
pg.de_analysis(all_data["GEX"], cluster='Lineage', de_key='de_final')
adata = all_data["GEX"].copy()
adata.raw = adata.copy()
pb.pseudobulk.calc_pseudobulk_stats(adata=adata, save_path="./allCells_PB", save_name="AllCells_V5", cluster_label='Lineage', sample_label="Sample", de_res_label='de_final', oneAVAtab=True, save_format='h5ad')

df = pd.read_excel("AllCells_V5_ava_pseudobulk_filt.xlsx")  # /projects/Pregnancy_placenta/Results/23_Supp_tables/DEG_markers/allCells_PB/
df = df.loc[df["Comparison"]=="HBCvsmMac"].copy()
df = df.sort_values("logFC")
mMac_genes = list(df.iloc[0:6,]["gene"])
HBC_genes = list(df.iloc[-6:,]["gene"])

func.feature_plots_pdf(all_data["GEX"], mMac_genes + HBC_genes, title="mMac_HBC_genes_hex", nrows=2, ncols=6, figsize=(10, 3), hex=True, legend_fontsize=8)
func.feature_plots_pdf(all_data["GEX"], mMac_genes + HBC_genes, title="mMac_HBC_genes3", nrows=2, ncols=6, figsize=(10, 3),  legend_fontsize=8, dot_size=0.03)



#=======================================================================
# Supp Fig 1O) Dotplot B/ Plasma cells
#=======================================================================
B_plasma_data = mu.read(data_dir+"/B_plasma.h5mu")
B_gex = B_plasma_data['GEX'].copy() # gene expression
B_ADT = B_plasma_data['ADT'].copy() # protein
strip_categories = func_dotplot.read_markers_dotplot(file_path= code_dir + "AuxFiles/marker genes_v3.xlsx",  lineage="B_Plasma")

clust_order = ['fB_1_naive', 'fB_3_cycling', 'mB_2_memory', 'mPlasma']
B_gex.obs["sub_lin2"] = pd.Categorical(B_gex.obs["sub_lin2"], ordered=True, categories=clust_order)
B_ADT.obs["sub_lin2"] = B_gex.obs.loc[B_ADT.obs_names, "sub_lin2"]

# reset GEX colors (now it is red-grey)
B_gex.uns['sub_lin2_colors']=['#1f77b4', '#ff7f0e', '#279e68', '#b5bd61']
sc.pl.umap(B_gex, color="sub_lin2")

func_dotplot.dot_plot_with_umap_stack_efficient(
    B_gex,
    prot_anndata= B_ADT,
    groupby="sub_lin2",
    strip_categories=strip_categories,
    standard_scale="var",
    legend_fontsize=16,
    legend_fontoutline=0,
    xtick_fontsize=12,
    category_size=12,
    umap_title_fontsize=18,
    legend_length=30,
    umap_title = "B_plasma",
    umap_aspect=2,
    umap_dotsize=15,
    legend_dotsize=100,
    figsize=(10, 4),
    cbar_orientation="vertical",
    cbar_loc=[1.05, 0.3, 0.015, 0.61],
    legend_loc=[.95, 0.05, 0.2, 0.1],
    save_fig="supp1_B_plasma_dotplot.pdf",
    dpi = 600)




# Trophoblasts
#====================
troph_subcluster = mu.read("Results/21_Final_objects/Trophoblasts.h5mu")
troph_subcluster["ADT"].obs['sub_lin2'] = troph_subcluster["GEX"][troph_subcluster["ADT"].obs_names,:].obs['sub_lin2'].copy()
troph_subcluster["ADT"].obs['sub_lin'] = troph_subcluster["GEX"][troph_subcluster["ADT"].obs_names,:].obs['sub_lin'].copy()
troph_adata = troph_subcluster['GEX'].copy()

troph_pdata = troph_subcluster['ADT'].copy()

# read markers
strip_cat_villous_troph = func_marker_dotplt.read_markers_dotplot(file_path = marker_file_path,
                                        lineage = "Villous_Troph")

# Plot
os.chdir("/projects/Pregnancy_placenta/Results/22_Final_figures_PAL/SupFig_Marker_Dotplots")
func_marker_dotplt.dot_plot_with_umap_stack(
    troph_adata,
    prot_anndata=troph_pdata,
    groupby="sub_lin2",
    strip_categories=strip_cat_villous_troph,
    shown_categories=['VCT_1', 'VCT_2', 'VCT_3', 'VCT_4', 'VCT_cycling_1', 'VCT_cycling_2', 'VCT_cycling_3', 'VCT_cycling_4', 'VCTfusion_1',
                      'VCTfusion_2', 'VCTccc', 'SCT'],
    legend_fontsize=14,
    legend_fontoutline=0,
    xtick_fontsize=16,
    category_size=16,
    legend_length=35,
    umap_title="Villous Troph Sub Lin.",
    umap_aspect=1.,
    legend_dotsize=80,
    cbar_orientation="vertical",
    cbar_loc=[1.05, 0.3, 0.015, 0.61],
    umap_title_fontsize=17,
    legend_loc=[0.955, 0.05, 0.10, 0.1],
    figsize=(18, 7),
    dpi=600
)




marker_dotplot_code_dir = code_dir + "Functions/"
if marker_dotplot_code_dir not in sys.path:
    sys.path.append(marker_dotplot_code_dir)

import Func_Marker_dotplot_fig as func_marker_dotplt
func_marker_dotplt = importlib.reload(func_marker_dotplt)

marker_file_path = code_dir + "AuxFiles/marker genes_v4_PAL.xlsx"
marker_dotplot_output_dir = "/projects/Pregnancy_placenta/Results/22_Final_figures_PAL/SupFig_Marker_Dotplots"


# Sup Fig 1P - Endothelial cells
#================================
endo_subcluster = mu.read(data_dir+"Endo.h5mu")
strip_cat_endo = func_marker_dotplt.read_markers_dotplot(
    file_path=marker_file_path,
    lineage="Endo")

os.chdir(marker_dotplot_output_dir)
func_marker_dotplt.dot_plot_with_umap_stack(
    endo_subcluster['GEX'],
    prot_anndata=endo_subcluster['ADT'],
    groupby="sub_lin2",
    strip_categories=strip_cat_endo,
    shown_categories=["fEndo_1", "fEndo_2", "fEndo_3", "fEndo_4", "mEndo_1", "mEndo_2", "mEndo_3", "mEndo_lymph", "fEndo_cycling"],
    legend_fontsize=16,
    legend_fontoutline=0,
    standard_scale="var",
    xtick_fontsize=16,
    category_size=16,
    legend_length=50,
    umap_title="All Endo",
    umap_aspect=1.2,
    legend_dotsize=80,
    cbar_orientation="vertical",
    cbar_loc=[1.05, 0.3, 0.015, 0.61],
    umap_title_fontsize=17,
    legend_loc=[0.955, 0.05, 0.10, 0.1],
    figsize=(18.5, 7),
    save_fig="SupFig_Endo_ALL_sublin2_V3.pdf",
    dpi=600)


# Sup Fig 1R - fFB
#==================
fFB_subcluster = mu.read(data_dir+"fFB.h5mu")
fFB_subcluster["ADT"].obs['short_sublin2'] = fFB_subcluster["GEX"][
    fFB_subcluster["ADT"].obs_names, :].obs['short_sublin2'].copy()
strip_cat_fFB = func_marker_dotplt.read_markers_dotplot(
    file_path=marker_file_path,
    lineage="fFB")

os.chdir(marker_dotplot_output_dir)
func_marker_dotplt.dot_plot_with_umap_stack(
    fFB_subcluster['GEX'],
    prot_anndata=fFB_subcluster['ADT'],
    groupby="short_sublin2",
    strip_categories=strip_cat_fFB,
    shown_categories=["fFB_1", "fFB_2", "fFB_3_1", "fFB_3_2", "fFB_4_1", "fFB_4_2", "fFB_cycling"],
    legend_fontsize=16,
    legend_fontoutline=0,
    xtick_fontsize=16,
    category_size=16,
    legend_length=20,
    umap_aspect=1.3,
    legend_dotsize=80,
    standard_scale="var",
    umap_title="fFB",
    cbar_orientation="vertical",
    cbar_loc=[1.05, 0.3, 0.015, 0.61],
    umap_title_fontsize=17,
    legend_loc=[0.955, 0.05, 0.10, 0.1],
    figsize=(18, 7),
    save_fig="SupFig_fFB_sublin2_dotplot_V2.pdf",
    dpi=600)


# Sup Fig 1S - Trophoblasts
#==========================
troph_subcluster = mu.read(data_dir+"Trophoblasts.h5mu")
troph_subcluster["ADT"].obs['sub_lin2'] = troph_subcluster["GEX"][
    troph_subcluster["ADT"].obs_names, :].obs['sub_lin2'].copy()
troph_subcluster["ADT"].obs['sub_lin'] = troph_subcluster["GEX"][
    troph_subcluster["ADT"].obs_names, :].obs['sub_lin'].copy()
troph_adata = troph_subcluster['GEX'].copy()
troph_pdata = troph_subcluster['ADT'].copy()
strip_cat_villous_troph = func_marker_dotplt.read_markers_dotplot(
    file_path=marker_file_path,
    lineage="Villous_Troph")

os.chdir(marker_dotplot_output_dir)
func_marker_dotplt.dot_plot_with_umap_stack(
    troph_adata,
    prot_anndata=troph_pdata,
    groupby="sub_lin2",
    strip_categories=strip_cat_villous_troph,
    shown_categories=['VCT_1', 'VCT_2', 'VCT_3', 'VCT_4', 'VCT_cycling_1', 'VCT_cycling_2', 'VCT_cycling_3', 'VCT_cycling_4', 'VCTfusion_1',
                      'VCTfusion_2', 'VCTccc', 'SCT'],
    legend_fontsize=14,
    legend_fontoutline=0,
    xtick_fontsize=16,
    category_size=16,
    legend_length=40,
    umap_title="Villous Troph Sub Lin.",
    umap_aspect=1,
    legend_dotsize=80,
    standard_scale="var",
    cbar_orientation="vertical",
    cbar_loc=[1.05, 0.3, 0.015, 0.61],
    umap_title_fontsize=17,
    legend_loc=[0.955, 0.05, 0.10, 0.1],
    figsize=(18, 7),
    save_fig="SupFig_Troph_ALL_sublin2_V3.pdf",
    dpi=600)


# Sup Fig 1T - fPV
#==================
fPV_subcluster = mu.read(data_dir+"fPV.h5mu")
strip_cat_fPV = func_marker_dotplt.read_markers_dotplot(
    file_path=marker_file_path,
    lineage="fPV")

os.chdir(marker_dotplot_output_dir)
func_marker_dotplt.dot_plot_with_umap_stack(
    fPV_subcluster['GEX'],
    prot_anndata=fPV_subcluster['ADT'],
    groupby="sub_lin2",
    strip_categories=strip_cat_fPV,
    shown_categories=["fPV_1_1", "fPV_1_2", "fPV_2_1", "fPV_2_2", "fPV_3_1", "fPV_3_2", "fPV_4_ISG15_IFI6", "fPV_cycling"],
    legend_fontsize=16,
    legend_fontoutline=0,
    xtick_fontsize=16,
    category_size=16,
    legend_length=35,
    umap_aspect=1.2,
    umap_title="fPV Sub Lin.",
    legend_dotsize=80,
    standard_scale="var",
    cbar_orientation="vertical",
    cbar_loc=[1.05, 0.3, 0.015, 0.61],
    umap_title_fontsize=17,
    legend_loc=[0.955, 0.05, 0.10, 0.1],
    figsize=(18, 7),
    save_fig="SupFig2D_fPVsublin_dotplot_V3.pdf",
    dpi=600)


# Sup Fig 1U - Decidual Stromal
#==============================
dS_subcluster = mu.read(data_dir+"dS.h5mu")
strip_cat_dS = func_marker_dotplt.read_markers_dotplot(
    file_path=marker_file_path,
    lineage="dS")

os.chdir(marker_dotplot_output_dir)
func_marker_dotplt.dot_plot_with_umap_stack(
    dS_subcluster['GEX'],
    prot_anndata=dS_subcluster['ADT'],
    groupby="sub_lin2",
    umap_dotsize=8,
    strip_categories=strip_cat_dS,
    shown_categories=['dS_1_1', 'dS_1_2', 'dS_2_1', 'dS_2_2', 'dS_2_3', 'dS_3_PRL', 'dS_4_1', 'dS_4_2', 'dS_5_MMP11', 'dS_lncRNA'],
    legend_fontsize=16,
    legend_fontoutline=0,
    xtick_fontsize=18,
    category_size=16,
    standard_scale="var",
    legend_length=35,
    umap_title="dS Sub Lin.",
    umap_aspect=1.1,
    legend_dotsize=80,
    cbar_orientation="vertical",
    cbar_loc=[1.05, 0.3, 0.015, 0.61],
    umap_title_fontsize=17,
    legend_loc=[0.955, 0.05, 0.10, 0.1],
    figsize=(18, 7),
    save_fig="SupFig_dS_ALL_sublin2_V3.pdf",
    dpi=600)


# Supp Figure 1X-Y
#=====================

os.chdir("/projects/Pregnancy_placenta/Results/25_Review/R2_C3/")
code_dir="/projects/Pregnancy_placenta/Code_Publication/Python/"

data = mu.read("/projects/Pregnancy_placenta/Results/21_Final_objects/GEO_objects/AllCells.h5mu")
obs = data["GEX"].obs[["Patient_ID", "Lineage", "Cell_subset", "Labor"]]  # only metadata is needed - no copying of the expression matrices
obs = obs[~obs["Lineage"].isin(["HSC", "Other"])].copy()

# (Supp Figure 1X) Barplot by condition
#====================================================================================================================
# Condition colors (same as in the paper figures)
conditions = ["Cont_labored", "Cont_unlabored", "Cont_PTB", "PTB", "DM1", "FGR", "PrePE", "TermPE"]
condition_colors = dict(zip(["Cont_labored", "Cont_unlabored", "Cont_PTB", "DM1", "FGR", "PrePE", "TermPE", "PTB"],
                            ["#CCEA8D", "#A6BC09", "green", "#5DB3D8", "#FFC817", "#EF7E29", "#E53938", "#8F2F82"]))
lineages = ["Trophoblasts", "EVT", "fPerivascular", "fFibroblasts", "dStromal", "Endothelial", "HBC", "Myeloid", "T", "NK", "B_Plasma"]

# Patients with several conditions are counted once per condition: the merge duplicates their cells (one row per condition)
pt_cond = pd.read_csv(code_dir + "AuxFiles/patient_condition_list.csv")
for pt in pt_cond["Patient_ID"][pt_cond["Patient_ID"].duplicated()].unique():
    print(f"{pt}: counted in conditions {list(pt_cond.loc[pt_cond['Patient_ID'] == pt, 'Condition_group'])}")
obs = obs.merge(pt_cond.rename(columns={"Condition_group": "Condition"}), on="Patient_ID", how="left")

# % of cells per subset from each condition (each bar sums to 100)
counts = pd.crosstab(obs["Cell_subset"], obs["Condition"])[conditions]
percent = counts.div(counts.sum(axis=1), axis=0) * 100

# Sort subsets by lineage order
subset_lineage = obs.drop_duplicates("Cell_subset").set_index("Cell_subset")["Lineage"].astype(str)
subset_lineage = subset_lineage.replace({"mMac": "Myeloid", "Monocytes": "Myeloid", "DC": "Myeloid"})  # not in the `lineages` list
order = sorted(percent.index, key=lambda s: (lineages.index(subset_lineage[s]), s))
percent = percent.loc[order]

fig, ax = plt.subplots(figsize=(max(8, 0.3 * len(percent)), 5))
percent.plot.bar(stacked=True, color=[condition_colors[c] for c in conditions], width=0.85, ax=ax)
ax.set_ylabel("% of cells in subset")
ax.set_xlabel("")
ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1), frameon=False)
ax.grid(False)
fig.savefig("Stacked_barplt_Cells_condition.pdf", bbox_inches="tight")


# ((Supp Figure 1Y) Barplot by labor status
#====================================================================================================================
# % of cells per subset from each labor status (each bar sums to 100)
counts = pd.crosstab(obs["Cell_subset"], obs["Labor"])
percent = counts.div(counts.sum(axis=1), axis=0) * 100

# Sort subsets by lineage order
subset_lineage = obs.drop_duplicates("Cell_subset").set_index("Cell_subset")["Lineage"].astype(str)
subset_lineage = subset_lineage.replace({"mMac": "Myeloid", "Monocytes": "Myeloid", "DC": "Myeloid"})  # not in the `lineages` list
order = sorted(percent.index, key=lambda s: (lineages.index(subset_lineage[s]), s))
percent = percent.loc[order]

fig, ax = plt.subplots(1,1,figsize=(max(8, 0.3 * len(percent)), 5))
percent.plot.bar(stacked=True, width=0.85, ax=ax)
ax.set_ylabel("% of cells in subset")
ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1), frameon=False)
ax.grid(False)
fig.savefig("Stacked_barplot_LaborStatus.pdf", bbox_inches="tight")


