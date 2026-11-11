library(openxlsx)
library(tidyr)
library(stringr)
library(ComplexHeatmap)
library(colorRamp2)

code_dir="/projects/Pregnancy_placenta/Code_Publication/R/"
source(paste0(code_dir, "Functions/FigFunctions_DA.R"))
setwd("/projects/Pregnancy_placenta/Results/22_Final_figures/Figure_3/")

all_da_stat = read.csv("/projects/Pregnancy_placenta/Results/22_Final_figures/DE_DA_results/all_DA_stats.csv")
myeloid_da_stat = all_da_stat[all_da_stat$Lineage=="Myeloid",]
myeloid_list= prep_df_stat_for_plotting(myeloid_da_stat)

log2fc_int_list = list(Myeloid = myeloid_list$log2_df)
qvalue_int_list = list(Myeloid = myeloid_list$qvalue_df)
pvalue_int_list = list(Myeloid = myeloid_list$pvalue_df)
myeloid_order=c("HBC_1","HBC_2","HBC_3","HBC_4","HBC_5","HBC_6","HBC_7","HBC_8","HBC_9","HBC_10","HBC_11","HBC_12","HBC_13","Mon_1_CD14","Mon_2_CD14","mMon_4_MS1","fMon_3_CD16","mMac_1","mMac_2","mMac_3","mMac_4","mMac_5","mMac_6","mMac_7","fDC2","mDC1","fDC3","mDC2/4","mDC_mature","fDC","pDC")

row_split_by=data.frame(c("HBC","HBC","HBC","HBC","HBC","HBC","HBC","HBC","HBC","HBC","HBC","HBC","HBC","Mon","Mon","Mon","Mon","mMac","mMac","mMac","mMac","mMac","mMac","mMac","DCs","DCs","DCs","DCs","DCs","DCs","DCs"), row.names = c("HBC_1","HBC_2","HBC_3","HBC_4","HBC_5","HBC_6","HBC_7","HBC_8","HBC_9","HBC_10","HBC_11","HBC_12","HBC_13","Mon_1_CD14","Mon_2_CD14","mMon_4_MS1","fMon_3_CD16","mMac_1","mMac_2","mMac_3","mMac_4","mMac_5","mMac_6","mMac_7","fDC2","mDC1","fDC3","mDC2/4","mDC_mature","fDC","pDC"))
colnames(row_split_by)="split"
row_split_by[,"split"] = factor(row_split_by[,"split"], levels=c("HBC","Mon","mMac","DCs"), ordered = TRUE)

DA_heatmaps(lineages_of_interest="Myeloid",
                 log2FC_list=log2fc_int_list,
                 qval_list=qvalue_int_list,
                 pval_list=pvalue_int_list,
                 plot_row_order = TRUE,
                 desired_row_order=myeloid_order,  # New argument for custom order
                 #order_comp2 = c("Cont_labored", "Cont_unlabored", "Cont_PTB"),
                 plot_type = "both",
                 save_pdf = TRUE,
                 output_prefix = "3B_myeloid_DA_heatmap",
                 output_dir = ".",
                 pdf_width = 7,
                 pdf_height = 8,
                 row_split_by=row_split_by)




