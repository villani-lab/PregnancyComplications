# 1)  Summarizing pathways By lineage (Rachelly's approach)
# 2)  Summarizing pathways regardless of cell subsets (for Andrea)
library(ggplot2)
library(reshape2)
library(xlsx)
library(ComplexHeatmap)
library(tidyr)
library(stringr)
library(dplyr)
source("/projects/Pregnancy_placenta/Code_R/Functions_IPA.R")

setwd("/projects/Pregnancy_placenta/Results/20_Figures/Pathways/")

#========================================================================
# 1) Summarize parent sublineages - Canonical Pathways / Diseases
#========================================================================
pw_dir = "/projects/Pregnancy_placenta/Code_PAL/DE_IPA_Code_Rmarkdown_PAL/IPA_results/Parent_lin/"
for (ipa_type in c("canonical_pathways", "diseases"))
{
  print(ipa_type)
  df_parent = data.frame()
  for(comp_no in c("1","2","3"))
  {
    df_curr = prep_IPA_df(dir_to_pathway = paste0(pw_dir,"/Parent_lin_test3_CompNum_",comp_no,"/"), pattern_int = ipa_type)
    df_curr$CompNum = comp_no
    df_parent = rbind(df_parent, df_curr)
  }
  write.table(df_parent, paste0("IPA_res_parent_lin_",ipa_type,".txt"))
}

# subset populations that don't have a "sub_lin" title
parent_sub = subset(df_parent,Cluster %in% c("SCT","VCTccc","Endo_mL","Endo_fProlif"))


# Summarize sublineages - Canonical Pathways / Diseases
#---------------------------------
pw_dir = "/projects/Pregnancy_placenta/Code_PAL/DE_IPA_Code_Rmarkdown_PAL/IPA_results/sub_lin/"

df_sublin = data.frame()
for(comp_no in c("1","2","3","4","5"))
{
  df_sublin_curr = prep_IPA_df(dir_to_pathway = paste0(pw_dir,"/Sublin_test3_CompNum_",comp_no,"/"), pattern_int = ipa_type)
  df_sublin_curr$CompNum = comp_no
  df_sublin = rbind(df_sublin, df_sublin_curr)
}

# Add cell subsets from parent  that are missing a "sub_lin"
df_sublin = rbind(df_sublin, parent_sub)


# Add lineage info to the shared_dist df
df_sublin["Lineage"] = "-"
df_sublin[grep("VCT|SCT",df_sublin$Cluster), "Lineage"] = "Trophoblasts"
df_sublin[grep("HBC",df_sublin$Cluster), "Lineage"] = "HBC"
df_sublin[grep("monocyte|cDC|mMac",df_sublin$Cluster), "Lineage"] = "Myeloid"
df_sublin[grep("Endo",df_sublin$Cluster), "Lineage"] = "Endothelial"
df_sublin[grep("dS|fPV|fFB",df_sublin$Cluster), "Lineage"] = "Stromal"

write.table(df_sublin, paste0("IPA_res_sublin_",ipa_type,".txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)




##############
# Summarize mid-level for endothelial cells and T cells
##############
ipa_type="diseases"
pw_dir = "/projects/Pregnancy_placenta/Code_PAL/DE_IPA_Code_Rmarkdown_PAL/IPA_results/Parent_lin/"

df_parent = data.frame()
for(comp_no in c("1","2","3"))
{
  df_curr = prep_IPA_df(dir_to_pathway = paste0(pw_dir,"/Parent_lin_test3_CompNum_",comp_no,"/"), pattern_int = ipa_type)
  df_curr$CompNum = comp_no
  df_parent = rbind(df_parent, df_curr)
}

# subset populations that don't have a "sub_lin" title
parent_sub = subset(df_parent,Cluster %in% c("SCT","VCTccc","Endo_mL","Endo_fProlif"))




###########################################################################################################
# 2) PTB vs Cont_PTB
###########################################################################################################
setwd('/projects/Pregnancy_placenta/Results/20_Figures/Pathways/PTB_vs_PTBcont/')
source("/projects/Pregnancy_placenta/Code_R/Functions_IPA.R")

# Merged - canonical pathways and diseases and biofunctions
df_canonical = read.table("../IPA_res_sublin_canonical_pathways.txt", sep="\t", header =1)
df_disease = read.table("../IPA_res_sublin_diseases.txt", sep="\t", header =1)
df = rbind(cbind(df_canonical, Type="Canonical_Pathways"), cbind(df_disease, Type="Diseases_and_biofunctions"))
PTB_df = df[(df$Comparison == "PTB_vs_Cont_PTB"),]

# rename cluysters by dictionary
PTB_df["Cluster"] = recode(PTB_df[,"Cluster"], !!!c("HBC_10_cycling"="HBC_10", "HBC_4_IL1B_MMP19"="HBC_4", "HBC_2_CXCL2_CXCL3"="HBC_2", "HBC_9_CXCL10_IFIT44"="HBC_9","mixed_HBC_13_S100A8"="HBC_13"))

plot_heatmap_PW_vs_comp(df_comp_merged=PTB_df,
                        zscore_min_cp = 2,
                        padj_thresh_cp = 0.01,
                        zscore_min_db = 3,
                        padj_thresh_db = 0.0001,
                        save_pdf = TRUE,
                        output_prefix = paste0("PTB_vs_PTBcont_PWannot_v3.pdf"),
                        orderColBy = "Cells",
                        pdf_width = 15,
                        pdf_height = 24,
                        cluster_names_col="Cluster",
                        filter_more="N",
                        plot_title="PTB vs PTB_cont",
                        pw_annot_file = "PTB_merged_Pathway_annotation_v3.xlsx")



plot_heatmap_PW_vs_comp(df_comp_merged=PTB_df,
                        zscore_min_cp = 2,
                        padj_thresh_cp = 0.01,
                        zscore_min_db = 3,
                        padj_thresh_db = 0.0001,
                        save_pdf = TRUE,
                        output_prefix = paste0("PTB_vs_PTBcont_PWannot_v4.pdf"),
                        orderColBy = "Cells",
                        pdf_width = 15,
                        pdf_height = 15,
                        cluster_names_col="Cluster",
                        filter_more="N",
                        plot_title="PTB vs PTBcont",
                        pw_annot_file = "PTB_merged_Pathway_annotation_v4_minimal.xlsx")

plot_heatmap_PW_vs_comp(df_comp_merged=PTB_df,
                        zscore_min_cp = 2,
                        padj_thresh_cp = 0.01,
                        zscore_min_db = 3,
                        padj_thresh_db = 0.0001,
                        save_pdf = TRUE,
                        output_prefix = paste0("PTB_vs_PTBcont_PWannot_v5.pdf"),
                        orderColBy = "Cells",
                        pdf_width = 15,
                        pdf_height = 15,
                        cluster_names_col="Cluster",
                        filter_more="N",
                        plot_title="PTB vs PTBcont",
                        pw_annot_file = "PTB_merged_Pathway_annotation_v5.xlsx")


###########################################################################################################
# 3) PE vs. unlabored
###########################################################################################################
setwd("/projects/Pregnancy_placenta/Results/20_Figures/Pathways/PE_vs_unlabored/")
source("/projects/Pregnancy_placenta/Code_R/Functions_IPA.R")

# Merged - canonical pathways and diseases and biofunctions
df_canonical = read.table("../IPA_res_sublin_canonical_pathways.txt", sep="\t", header =1)
df_disease = read.table("../IPA_res_sublin_diseases.txt", sep="\t", header =1)
df = rbind(cbind(df_canonical, Type="Canonical_Pathways"), cbind(df_disease, Type="Diseases_and_biofunctions"))
PE_df = df[df$Comparison %in% c("PrePE_vs_Cont_unlabored", "TermPE_vs_Cont_unlabored"),]
PE_df["Cluster"] = recode(PE_df[,"Cluster"], !!!c("HBC_10_cycling"="HBC_10", "HBC_4_IL1B_MMP19"="HBC_4", "HBC_2_CXCL2_CXCL3"="HBC_2", "HBC_9_CXCL10_IFIT44"="HBC_9","mixed_HBC_13_S100A8"="HBC_13"))

plot_heatmap_PW_vs_comp(df_comp_merged=PE_df,
                        zscore_min_cp = 2,
                        padj_thresh_cp = 0.05,
                        zscore_min_db = 2,
                        padj_thresh_db = 0.001,
                        save_pdf = TRUE,
                        output_prefix = paste0("PE_vs_unlabored_PWannot_v2b.pdf"),
                        orderColBy = "Cells",
                        pdf_width = 15,
                        pdf_height = 24,
                        cluster_names_col="Cluster",
                        filter_more="N",
                        plot_title="PE vs unlabored_cont",
                        pw_annot_file = "PE_merged_PW_annotation_v2.xlsx")

plot_heatmap_PW_vs_comp(df_comp_merged=PE_df,
                        zscore_min_cp = 2,
                        padj_thresh_cp = 0.05,
                        zscore_min_db = 2,
                        padj_thresh_db = 0.001,
                        save_pdf = TRUE,
                        output_prefix = paste0("PE_vs_unlabored_PWannot_v3.pdf"),
                        orderColBy = "Cells",
                        pdf_width = 15,
                        pdf_height = 24,
                        cluster_names_col="Cluster",
                        filter_more="N",
                        plot_title="PE vs unlabored_cont",
                        pw_annot_file = "PE_merged_PW_annotation_v3.xlsx")

plot_heatmap_PW_vs_comp(df_comp_merged=PE_df,
                        zscore_min_cp = 2,
                        padj_thresh_cp = 0.05,
                        zscore_min_db = 2,
                        padj_thresh_db = 0.001,
                        save_pdf = TRUE,
                        output_prefix = paste0("PE_vs_unlabored_PWannot_v4.pdf"),
                        orderColBy = "Cells",
                        pdf_width = 15,
                        pdf_height = 15,
                        cluster_names_col="Cluster",
                        filter_more="N",
                        plot_title="PE vs unlabored_cont",
                        pw_annot_file = "/projects/Pregnancy_placenta/Results/20_Figures/Pathways/PE_vs_unlabored/myeloid/PE_merged_PW_annotation_v4.xlsx")


# Lets have all myeloid cells under "myeloid"
source("/projects/Pregnancy_placenta/Code_R/Functions_IPA.R")
PE_df[PE_df["Lineage"]=="HBC","Lineage"] = "Myeloid"
plot_heatmap_PW_vs_comp(df_comp_merged=PE_df,
                        zscore_min_cp = 2,
                        padj_thresh_cp = 0.05,
                        zscore_min_db = 2,
                        padj_thresh_db = 0.001,
                        save_pdf = TRUE,
                        output_prefix = paste0("PE_vs_unlabored_PWannot_v4_ColSort"),
                        orderColBy = "Cells",
                        order_vector=c("HBC_2_TermPE","HBC_3_TermPE","CD14_monocyte_TermPE","HBC_4_TermPE","HBC_13_TermPE","HBC_2_PrePE","HBC_3_PrePE","CD14_monocyte_PrePE","HBC_1_PrePE","SCT_TermPE","VCT_1_TermPE","VCT_3_TermPE","VCTfusion_1_TermPE","VCTfusion_2_TermPE","VCTp_1_TermPE","VCTp_2_TermPE","SCT_PrePE","VCT_1_PrePE","VCT_3_PrePE","VCTfusion_1_PrePE","VCTfusion_2_PrePE","VCTccc_PrePE","Endo_f1_PrePE","Endo_f2_PrePE","Endo_f3_PrePE"),
                        pdf_width = 15,
                        pdf_height = 15,
                        cluster_names_col="Cluster",
                        filter_more="N",
                        plot_title="PE vs unlabored_cont",
                        pw_annot_file = "/projects/Pregnancy_placenta/Results/20_Figures/Pathways/PE_vs_unlabored/PE_merged_PW_annotation_v4.xlsx")

plot_heatmap_PW_vs_comp(df_comp_merged=PE_df,
                        zscore_min_cp = 2,
                        padj_thresh_cp = 0.05,
                        zscore_min_db = 2,
                        padj_thresh_db = 0.001,
                        save_pdf = TRUE,
                        output_prefix = paste0("PE_vs_unlabored_PWannot_v5"),
                        orderColBy = "Cells",
                        order_vector=c("HBC_2_TermPE","HBC_3_TermPE","CD14_monocyte_TermPE","HBC_4_TermPE","HBC_13_TermPE","HBC_2_PrePE","HBC_3_PrePE","CD14_monocyte_PrePE","HBC_1_PrePE","SCT_TermPE","VCT_1_TermPE","VCT_3_TermPE","VCTfusion_1_TermPE","VCTfusion_2_TermPE","VCTp_1_TermPE","VCTp_2_TermPE","SCT_PrePE","VCT_1_PrePE","VCT_3_PrePE","VCTfusion_1_PrePE","VCTfusion_2_PrePE","VCTccc_PrePE","Endo_f1_PrePE","Endo_f2_PrePE","Endo_f3_PrePE"),
                        pdf_width = 15,
                        pdf_height = 15,
                        cluster_names_col="Cluster",
                        filter_more="N",
                        plot_title="PE vs unlabored_cont",
                        pw_annot_file = "/projects/Pregnancy_placenta/Results/20_Figures/Pathways/PE_vs_unlabored/PE_merged_PW_annotation_v5.xlsx")

PE_df[PE_df["Lineage"]=="HBC","Lineage"] = "Myeloid"
plot_heatmap_PW_vs_comp(df_comp_merged=PE_df,
                        zscore_min_cp = 2,
                        padj_thresh_cp = 0.05,
                        zscore_min_db = 2,
                        padj_thresh_db = 0.001,
                        save_pdf = TRUE,
                        output_prefix = paste0("PE_vs_unlabored_PWannot_v6"),
                        orderColBy = "Cells",
                        order_vector=c("HBC_2_TermPE","HBC_3_TermPE","CD14_monocyte_TermPE","HBC_4_TermPE","HBC_13_TermPE","HBC_2_PrePE","HBC_3_PrePE","CD14_monocyte_PrePE","HBC_1_PrePE","SCT_TermPE","VCT_1_TermPE","VCT_3_TermPE","VCTfusion_1_TermPE","VCTfusion_2_TermPE","VCTp_1_TermPE","VCTp_2_TermPE","SCT_PrePE","VCT_1_PrePE","VCT_3_PrePE","VCTfusion_1_PrePE","VCTfusion_2_PrePE","VCTccc_PrePE","Endo_f1_PrePE","Endo_f2_PrePE","Endo_f3_PrePE"),
                        pdf_width = 15,
                        pdf_height = 15,
                        cluster_names_col="Cluster",
                        filter_more="N",
                        plot_title="PE vs unlabored_cont",
                        pw_annot_file = "/projects/Pregnancy_placenta/Results/20_Figures/Pathways/PE_vs_unlabored/PE_merged_PW_annotation_v6.xlsx")





# PW heatmap for HBC/ Myeloid cells only
#=============================================
setwd("/projects/Pregnancy_placenta/Results/20_Figures/Pathways/PE_vs_unlabored/myeloid/")
# Merged - canonical pathways and diseases and biofunctions
df_canonical = read.table("/projects/Pregnancy_placenta/Results/20_Figures/Pathways/IPA_res_sublin_canonical_pathways.txt", sep="\t", header =1)
df_disease = read.table("/projects/Pregnancy_placenta/Results/20_Figures/Pathways/IPA_res_sublin_diseases.txt", sep="\t", header =1)
df = rbind(cbind(df_canonical, Type="Canonical_Pathways"), cbind(df_disease, Type="Diseases_and_biofunctions"))
PE_df = df[df$Comparison %in% c("PrePE_vs_Cont_unlabored", "TermPE_vs_Cont_unlabored"),]
PE_df["Cluster"] = recode(PE_df[,"Cluster"], !!!c("HBC_10_cycling"="HBC_10", "HBC_4_IL1B_MMP19"="HBC_4", "HBC_2_CXCL2_CXCL3"="HBC_2", "HBC_9_CXCL10_IFIT44"="HBC_9","mixed_HBC_13_S100A8"="HBC_13"))

# Filter for myeloid cells only
myeloid_df = PE_df[(PE_df["Lineage"]=="Myeloid") | PE_df["Lineage"]=="HBC",]
myeloid_df["Lineage"]="Myeloid" # I prefer to show all the cells together and split by condition
plot_heatmap_PW_vs_comp(df_comp_merged=myeloid_df,
                        zscore_min_cp = 1,
                        padj_thresh_cp = 0.05,
                        zscore_min_db = 1,
                        padj_thresh_db = 0.05,
                        save_pdf = TRUE,
                        output_prefix = paste0("PE_vs_unlabored_v4_myeloid_lowerThresh.pdf"),
                        orderColBy = "Cells",
                        pdf_width = 15,
                        pdf_height = 15,
                        cluster_names_col="Cluster",
                        filter_more="N",
                        plot_title="PE vs unlabored_cont",
                        pw_annot_file = "PE_merged_PW_annotation_v4_myeloid.xlsx")
