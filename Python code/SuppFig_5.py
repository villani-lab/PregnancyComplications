# supp fig 5 has all the comparisons against all other controls

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

sys.path.append("/projects/Pregnancy_placenta/Code/Publication_code")
import FigFunctions as fig_func

# default lineage palette you provided
_DEF_LINEAGES = ['Trophoblasts','EVT','fPV','fFB','dS','Endothelial','HBC','Myeloid','T','NK','B_Plasma','HSC','Other']
_DEF_COLORS = ["#60B364","#B2DF8A",  # trophoblasts
               "#08306b","#2b8cbe","#c6dbef","#8F6FE8",  # stromal/endo
               "#800026","#f03b20","#fd8d3c","#fed976","#fa9fb5",  # immune
               "grey","lightgrey"]
_DEF_PALETTE = dict(zip(_DEF_LINEAGES, _DEF_COLORS))

work_dir = "/projects/Pregnancy_placenta/Results/22_Final_figures/SuppFigure_5/"
os.chdir(work_dir)


#========================================================================================
# Sup Fig 5A) DA across all comparisons
#========================================================================================
comparisons = ["PTB_vs_Cont_labored",
               "PrePE_vs_Cont_labored",
               "TermPE_vs_Cont_labored",
               "DM1_vs_Cont_labored",
               "FGR_vs_Cont_labored",
               "PrePE_vs_Cont_PTB",
               "Cont_labored_vs_Cont_unlabored"]


fig_func.multiDA_BarGrid(
    comparisons=comparisons,
    pdf_name="Supp_DA_labored_and_PTBcont_onlySig.pdf",
    path_to_stats="/projects/Pregnancy_placenta/Code/Publication_code/AuxFiles/",
    da_pval_thresh=0.05,  # This threshold is used to have a transparent fill for populations with pval< <thresh>
    onlyMinDAPop=True,   # include only populations with a (p_val<0.05) AND (FC>minDA_absFC_thresh) in ANY comparison
    minDA_absFC_thresh=0, # threshold applied to abs(log(FC))
    share_y=False,
    xticks_fontsize=10,
    width=8,
    global_order=['VCT_1', 'VCT_2', 'VCT_cycling_3', 'VCT_cycling_4', 'VCTfusion_2', 'EVT_1', 'EVT_6',
                  'HBC_1', 'HBC_2', 'HBC_5', 'HBC_7', 'HBC_8', 'HBC_11', 'HBC_13',
                  'mDC2/4', 'pDC', 'mMac_1', 'mMac_2', 'mMac_3', 'mMac_4', 'mMac_6',
                   'mT_1_CD4_eff', 'fT_3_CD4_naive', 'mT_6_Treg','mT_12_MAIT',
                  'fNKT_1', 'mNK_2', 'fNK_4',
                  'fPV_2_1', 'fPV_2_2', 'fPV_4', 'fFB_1', 'fFB_3_1', 'fFB_3_2', 'fFB_4_1', 'fFB_cycling',
                  'dS_4_2', 'dS_5', 'dS_lncRNA', 'fEndo_1', 'mEndo_3', 'fEndo_4', 'fEndo_cycling', 'mEndo_lymph'])


#========================================================================================
# Sup Fig 5B) DE across all comparisons
#========================================================================================
fig_func.multiDE_BarGrid(comparisons=comparisons,
                         pval_adj_thresh=0.15,
                         logFC_thresh=0.58,
                         pdf_name="Supp_DE_labored_and_PTBcont_onlySig_min50DEGs.pdf",
                         path_to_stats="/projects/Pregnancy_placenta/Code/Publication_code/AuxFiles/",
                         xticks_fontsize=10,
                         mimDEGs_thresh=50, # include only populations with a minimum <mimDEGs_thresh> DEGs
                         onlyMinDEGsPop=True,
                         DEGs_transparency_thresh=0, # A threshshold below which the bar will be partially transparent
                         width=8)
                         #global_order=['SCT', 'VCT_1', 'VCT_2', 'VCT_3', 'VCT_cycling_1', 'VCT_cycling_2', 'VCTfusion_1', 'EVT_2', 'HBC_1',  'HBC_2', 'HBC_3', 'HBC_4', 'HBC_5', 'HBC_6', 'HBC_8', 'HBC_9','HBC_10', 'HBC_11', 'Mon_1_CD14', 'mMac_1', 'mMac_2', 'mMac_3', 'mMac_4', 'mMac_6', 'mNK_2', 'fPV_1_1', 'fPV_2_1', 'fPV_3_2', 'fFB_1', 'dS_1_1', 'dS_3', 'fEndo_1', 'mEndo_1', 'fEndo_3'])


