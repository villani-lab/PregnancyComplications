library(ggplot2)
library(reshape2)
library(xlsx)
library(ComplexHeatmap)
library(tidyr)
library(stringr)
library(dplyr)

code_path="/projects/Pregnancy_placenta/Code_R/Figures/"
code_dir="/projects/Pregnancy_placenta/Code_Publication/R/"
source(paste0(code_dir, "Functions/figure_functions.R"))
setwd("/projects/Pregnancy_placenta/Results/22_Final_figures/Figure_5/")


#========================================================
# Figure 5B - pathway enrichment in term- and preterm-PE
#========================================================

df_canonical = read.table(paste0(code_path, "AUX_files/IPA_res_sublin_canonical_pathways.txt"), sep="\t", header =1, quote = "", comment.char = "")
df_disease = read.table(paste0(code_path, "AUX_files/IPA_res_sublin_diseases.txt"), sep="\t", header =1, quote = "", comment.char = "")
df = rbind(cbind(df_canonical, Type="Canonical_Pathways"), cbind(df_disease, Type="Diseases_and_biofunctions"))
PE_df = df[df$Comparison %in% c("PrePE_vs_Cont_unlabored", "TermPE_vs_Cont_unlabored"),]
PE_df[PE_df["Lineage"]=="HBC","Lineage"] = "Myeloid"

plot_heatmap_PW(df_comp_merged=PE_df,
                        zscore_min_cp = 2,
                        padj_thresh_cp = 0.05,
                        zscore_min_db = 2,
                        padj_thresh_db = 0.001,
                        save_pdf = TRUE,
                        output_prefix = paste0("Fig5B_PE_vs_unlabored_PWannot"),
                        orderColBy = "Cells",
                        order_vector=c("HBC_2_TermPE","HBC_3_TermPE","CD14_monocyte_TermPE","HBC_4_TermPE","HBC_13_TermPE","HBC_2_PrePE","HBC_3_PrePE","CD14_monocyte_PrePE","HBC_1_PrePE","SCT_TermPE","VCT_1_TermPE","VCT_3_TermPE","VCTfusion_1_TermPE","VCTfusion_2_TermPE","VCT_cycling_1_TermPE","VCT_cycling_2_TermPE","SCT_PrePE","VCT_1_PrePE","VCT_3_PrePE","VCTfusion_1_PrePE","VCTfusion_2_PrePE","VCTccc_PrePE","fEndo_1_PrePE","fEndo_2_PrePE","fEndo_3_PrePE"),
                        pdf_width = 15,
                        pdf_height = 15,
                        cluster_names_col="Cluster",
                        filter_more="N",
                        plot_title="PE vs unlabored_cont",
                        pw_annot_file = paste0(code_path,"/AUX_files/PE_merged_PW_annotation_v6.xlsx"))



