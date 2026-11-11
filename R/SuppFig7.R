source("/projects/Pregnancy_placenta/Code_R/Figures/figure_functions.R")
setwd("/projects/Pregnancy_placenta/Results/22_Final_figures/SuppFigure_7/DEG_count/")

# 1) DE summary
DEGs_count_df = read.csv(file = "/projects/Pregnancy_placenta/Code/Publication_code/AuxFiles/DEGs_count_pval0.15_abslog2FC0.58.csv")
DEGs_myeloid = DEGs_count_df[DEGs_count_df[["Lineage"]] %in% c("DC","HBC","Mon", "PAMM"),]
lin_subset_degs = unique(DEGs_myeloid[,c("Subset", "Lineage")])
colnames(lin_subset_degs) = c("subset","Lineage")
rownames(lin_subset_degs) = lin_subset_degs$subset
clust_order=c("HBC_1_PLA2G7","HBC_2","HBC_3","HBC_4","HBC_5","HBC_6","HBC_7_lncRNA","HBC_8_neuro","HBC_9_CXCL10_IFIT1","HBC_10_cycling","HBC_11","HBC_12_CXCL10","HBC_13_S100A8",
                                 "mMac_1_C3_LTC4S","mMac_2_MMP7_CHI3L1","mMac_3_MT-CO1","mMac_4_NUPR1","mMac_5_CXCL10_IDO1","mMac_6_CX3CR1","mMac_7_APOE","Mon_1_CD14_CLEC4D","Mon_2_CD14_CXCL1","fMon_3_CD16","mMon_4_MS1","fDC_CLECL10A_S100A8","pDC")
DEGS_heatmp_wrapper(DEGs_count_df=DEGs_myeloid, lin_subset_degs=lin_subset_degs,
                    output_str="myeloid_level3_blackWhite",
                    clust_order=clust_order)


# 2) DA summary
DA_all_res = read.csv("/projects/Pregnancy_placenta/Code/Publication_code/AuxFiles/DA_all_pval.csv", check.names = F, row.names = "index")       # based on p-value
#DA_all_res = read.csv("/projects/Pregnancy_placenta/Code/Publication_code/AuxFiles/DA_all_qval.csv", check.names = F, row.names = "index")       # based on p-value

# Making q-values binary
DA_all_res[DA_all_res < 0.05] = 0 # significant
DA_all_res[DA_all_res >= 0.05] = 1 # not significant

# Filtering to myeloid only
#myeloid_clust = c("fDC_CLECL10A_S100A8","fDC2","fDC3","HBC_1_PLA2G7","HBC_10_cycling", "HBC_11","HBC_2","HBC_3","HBC_4","HBC_5", "HBC_6","HBC_7_lncRNA","HBC_8_neuro","HBC_9_CXCL10_IFIT1", "mDC1",  "mDC2/4","mMon_4_MS1","Mon_1_CD14_CLEC4D","Mon_2_CD14_CXCL1","fMon_3_CD16",  "HBC_12_CXCL10" ,"HBC_13_S100A8","pDC","mMac_1_C3_LTC4S","mMac_2_MMP7_CHI3L1", "mMac_3_MT-CO1","mMac_4_NUPR1","mMac_5_CXCL10_IDO1","mMac_6_CX3CR1","mMac_7_APOE", "mDC_mature")
DA_all_res = DA_all_res[clust_order,]

source("/projects/Pregnancy_placenta/Code_R/Figures/figure_functions.R")
mat_values = c(0, max(DA_all_res))
plot_heatmap(df=DA_all_res, lin_subset_df=lin_subset_degs, output_prefix="myeloid_DA_summary", mat_name = "pval",mat_values=mat_values, mat_colors= c("black", "white"), cluster_rows=F, pdf_width=7, pdf_height=15)


