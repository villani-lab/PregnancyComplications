# Supplemental Figure 9: DM1 and FGR

# Section 1: Load helper functions

code_dir="/projects/Pregnancy_placenta/Code_Publication/R/"
source(paste0(code_dir, "Functions/IPA_pathway_heatmap_functions.R"))


# Section 2: Prepare Pathway Dataframes
## Section 2A: Work with Diseases

df_parent_test3_comp1 = prep_IPA_df(dir_to_pathway = "IPA_results/Parent_lin/Parent_lin_test3_CompNum_1",
           pattern_int = "diseases")
df_parent_test3_comp1$CompNum = 1
df_parent_test3_comp2 = prep_IPA_df(dir_to_pathway = "IPA_results/Parent_lin/Parent_lin_test3_CompNum_2",
           pattern_int = "diseases")
df_parent_test3_comp2$CompNum = 2
df_parent_test3_comp3 = prep_IPA_df(dir_to_pathway = "IPA_results/Parent_lin/Parent_lin_test3_CompNum_3",
           pattern_int = "diseases")
df_parent_test3_comp3$CompNum = 3

df_parent_test3_comp4 = prep_IPA_df(dir_to_pathway = "IPA_results/Parent_lin/Parent_lin_test3b_CompNum_1/",
           pattern_int = "diseases")
df_parent_test3_comp4$CompNum = 4

df_parent_test3_merged = rbind(df_parent_test3_comp1,df_parent_test3_comp2,df_parent_test3_comp3,df_parent_test3_comp4)

# Add in main lin
df_parent_test3_merged$Main_lin = "-"
df_parent_test3_merged[grepl("VCT|SCT|VCTccc|SCT",df_parent_test3_merged$Cluster),]$Main_lin = "Troph"
df_parent_test3_merged[grepl("HBC",df_parent_test3_merged$Cluster),]$Main_lin = "HBC"
df_parent_test3_merged[grepl("Monocytes",df_parent_test3_merged$Cluster),]$Main_lin = "Monocytes"
df_parent_test3_merged[grepl("Endo_f",df_parent_test3_merged$Cluster),]$Main_lin = "EndoF"

df_parent_test3_merged[grepl("Endo_m",df_parent_test3_merged$Cluster),]$Main_lin = "EndoM"
df_parent_test3_merged[grepl("fFB",df_parent_test3_merged$Cluster),]$Main_lin = "fFB"
df_parent_test3_merged[grepl("fPV",df_parent_test3_merged$Cluster),]$Main_lin = "fPV"
df_parent_test3_merged[grepl("NK",df_parent_test3_merged$Cluster),]$Main_lin = "NK"

df_parent_test3_merged[grepl("DC",df_parent_test3_merged$Cluster),]$Main_lin = "DC"


df_parent_test3_merged[grepl("PAMM",df_parent_test3_merged$Cluster),]$Main_lin = "PAMM"
df_parent_test3_merged[grepl("Tcells",df_parent_test3_merged$Cluster),]$Main_lin = "Tcells"
df_parent_test3_merged[grepl("B_cells",df_parent_test3_merged$Cluster),]$Main_lin = "Bcells"



df_parent_test3_merged$Cluster = paste0(df_parent_test3_merged$Cluster,"_L1")

parent_in_sublin = subset(df_parent_test3_merged,Cluster %in% c("SCT_L1","VCTccc_L1"))

df_sublin_test3_comp1 = prep_IPA_df(dir_to_pathway = "IPA_results/sub_lin/Sublin_test3_CompNum_1/",pattern_int = "diseases")
df_sublin_test3_comp1$CompNum = 1

df_sublin_test3_comp2 = prep_IPA_df(dir_to_pathway = "IPA_results/sub_lin/Sublin_test3_CompNum_2/",pattern_int = "diseases")
df_sublin_test3_comp2$CompNum = 2

df_sublin_test3_comp3 = prep_IPA_df(dir_to_pathway = "IPA_results/sub_lin/Sublin_test3_CompNum_3/",pattern_int = "diseases")
df_sublin_test3_comp3$CompNum = 3

df_sublin_test3_comp4 = prep_IPA_df(dir_to_pathway = "IPA_results/sub_lin/Sublin_test3_CompNum_4/",pattern_int = "diseases")
df_sublin_test3_comp4$CompNum = 4

df_sublin_test3_comp5 = prep_IPA_df(dir_to_pathway = "IPA_results/sub_lin/Sublin_test3_CompNum_5/",pattern_int = "diseases")
df_sublin_test3_comp5$CompNum = 5
df_sublin_test3_merged = rbind(df_sublin_test3_comp1,df_sublin_test3_comp2,df_sublin_test3_comp3,df_sublin_test3_comp4,df_sublin_test3_comp5)

df_sublin_test3_merged$Main_lin = "-"
df_sublin_test3_merged[grepl("VCT|SCT|VCTccc",df_sublin_test3_merged$Cluster),]$Main_lin = "Troph"
df_sublin_test3_merged[grepl("HBC",df_sublin_test3_merged$Cluster),]$Main_lin = "HBC"
df_sublin_test3_merged[grepl("monocyte",df_sublin_test3_merged$Cluster),]$Main_lin = "Monocytes"
df_sublin_test3_merged[grepl("Endo_f",df_sublin_test3_merged$Cluster),]$Main_lin = "EndoF"
df_sublin_test3_merged[grepl("Endo_m",df_sublin_test3_merged$Cluster),]$Main_lin = "EndoM"
df_sublin_test3_merged[grepl("fFB",df_sublin_test3_merged$Cluster),]$Main_lin = "fFB"
df_sublin_test3_merged[grepl("fPV",df_sublin_test3_merged$Cluster),]$Main_lin = "fPV"
df_sublin_test3_merged[grepl("cDC",df_sublin_test3_merged$Cluster),]$Main_lin = "DC"
df_sublin_test3_merged[grepl("mMac",df_sublin_test3_merged$Cluster),]$Main_lin = "PAMM"
df_sublin_test3_merged[grepl("dS",df_sublin_test3_merged$Cluster),]$Main_lin = "dS"

df_sublin_test3_merged$Cluster = paste0(df_sublin_test3_merged$Cluster,"_L2")



df_sublin_test3_merged = rbind(df_sublin_test3_merged,parent_in_sublin)

saveRDS(df_sublin_test3_merged,file = "rds_objects/IPA_pathways/final_objs/sublin_diseases_V1.rds")


## Section 2B: Work with Canonical Pathways
df_parent_test3_comp1 = prep_IPA_df(dir_to_pathway = "IPA_results/Parent_lin/Parent_lin_test3_CompNum_1",
           pattern_int = "canonical_pathways")
df_parent_test3_comp1$CompNum = 1
df_parent_test3_comp2 = prep_IPA_df(dir_to_pathway = "IPA_results/Parent_lin/Parent_lin_test3_CompNum_2",
           pattern_int = "canonical_pathways")
df_parent_test3_comp2$CompNum = 2
df_parent_test3_comp3 = prep_IPA_df(dir_to_pathway = "IPA_results/Parent_lin/Parent_lin_test3_CompNum_3",
           pattern_int = "canonical_pathways")
df_parent_test3_comp3$CompNum = 3

df_parent_test3_comp4 = prep_IPA_df(dir_to_pathway = "IPA_results/Parent_lin/Parent_lin_test3b_CompNum_1/",
           pattern_int = "canonical_pathways")
df_parent_test3_comp4$CompNum = 4

df_parent_test3_merged = rbind(df_parent_test3_comp1,df_parent_test3_comp2,df_parent_test3_comp3,df_parent_test3_comp4)

# Add in main lin
df_parent_test3_merged$Main_lin = "-"
df_parent_test3_merged[grepl("VCT|SCT|VCTccc|SCT",df_parent_test3_merged$Cluster),]$Main_lin = "Troph"
df_parent_test3_merged[grepl("HBC",df_parent_test3_merged$Cluster),]$Main_lin = "HBC"
df_parent_test3_merged[grepl("Monocytes",df_parent_test3_merged$Cluster),]$Main_lin = "Monocytes"
df_parent_test3_merged[grepl("Endo_f",df_parent_test3_merged$Cluster),]$Main_lin = "EndoF"

df_parent_test3_merged[grepl("Endo_m",df_parent_test3_merged$Cluster),]$Main_lin = "EndoM"
df_parent_test3_merged[grepl("fFB",df_parent_test3_merged$Cluster),]$Main_lin = "fFB"
df_parent_test3_merged[grepl("fPV",df_parent_test3_merged$Cluster),]$Main_lin = "fPV"
df_parent_test3_merged[grepl("NK",df_parent_test3_merged$Cluster),]$Main_lin = "NK"

df_parent_test3_merged[grepl("DC",df_parent_test3_merged$Cluster),]$Main_lin = "DC"


df_parent_test3_merged[grepl("PAMM",df_parent_test3_merged$Cluster),]$Main_lin = "PAMM"
df_parent_test3_merged[grepl("Tcells",df_parent_test3_merged$Cluster),]$Main_lin = "Tcells"
df_parent_test3_merged[grepl("B_cells",df_parent_test3_merged$Cluster),]$Main_lin = "Bcells"


df_parent_test3_merged$Cluster = paste0(df_parent_test3_merged$Cluster,"_L1")


df_sublin_test3_comp1 = prep_IPA_df(dir_to_pathway = "IPA_results/sub_lin/Sublin_test3_CompNum_1/",pattern_int = "canonical_pathways")
df_sublin_test3_comp1$CompNum = 1

df_sublin_test3_comp2 = prep_IPA_df(dir_to_pathway = "IPA_results/sub_lin/Sublin_test3_CompNum_2/",pattern_int = "canonical_pathways")
df_sublin_test3_comp2$CompNum = 2

df_sublin_test3_comp3 = prep_IPA_df(dir_to_pathway = "IPA_results/sub_lin/Sublin_test3_CompNum_3/",pattern_int = "canonical_pathways")
df_sublin_test3_comp3$CompNum = 3

df_sublin_test3_comp4 = prep_IPA_df(dir_to_pathway = "IPA_results/sub_lin/Sublin_test3_CompNum_4/",pattern_int = "canonical_pathways")
df_sublin_test3_comp4$CompNum = 4

df_sublin_test3_comp5 = prep_IPA_df(dir_to_pathway = "IPA_results/sub_lin/Sublin_test3_CompNum_5/",pattern_int = "canonical_pathways")
df_sublin_test3_comp5$CompNum = 5
df_sublin_test3_merged = rbind(df_sublin_test3_comp1,df_sublin_test3_comp2,df_sublin_test3_comp3,df_sublin_test3_comp4,df_sublin_test3_comp5)

parent_in_sublin = subset(df_parent_test3_merged,Cluster %in% c("SCT_L1","VCTccc_L1"))


df_sublin_test3_merged$Main_lin = "-"
df_sublin_test3_merged[grepl("VCT|SCT|VCTccc",df_sublin_test3_merged$Cluster),]$Main_lin = "Troph"
df_sublin_test3_merged[grepl("HBC",df_sublin_test3_merged$Cluster),]$Main_lin = "HBC"
df_sublin_test3_merged[grepl("monocyte",df_sublin_test3_merged$Cluster),]$Main_lin = "Monocytes"
df_sublin_test3_merged[grepl("Endo_f",df_sublin_test3_merged$Cluster),]$Main_lin = "EndoF"
df_sublin_test3_merged[grepl("Endo_m",df_sublin_test3_merged$Cluster),]$Main_lin = "EndoM"
df_sublin_test3_merged[grepl("fFB",df_sublin_test3_merged$Cluster),]$Main_lin = "fFB"
df_sublin_test3_merged[grepl("fPV",df_sublin_test3_merged$Cluster),]$Main_lin = "fPV"
df_sublin_test3_merged[grepl("cDC",df_sublin_test3_merged$Cluster),]$Main_lin = "DC"
df_sublin_test3_merged[grepl("mMac",df_sublin_test3_merged$Cluster),]$Main_lin = "PAMM"
df_sublin_test3_merged[grepl("dS",df_sublin_test3_merged$Cluster),]$Main_lin = "dS"

df_sublin_test3_merged$Cluster = paste0(df_sublin_test3_merged$Cluster,"_L2")




parent_in_sublin = subset(df_parent_test3_merged,Cluster %in% c("SCT_L1","VCTccc_L1"))
df_sublin_test3_merged = rbind(df_sublin_test3_merged,parent_in_sublin)

saveRDS(df_sublin_test3_merged,file = "rds_objects/IPA_pathways/final_objs/sublin_canonical_pathways_V1.rds")

# Section 3: Prepare DF for plotting (start here if var prepped)

## load each object
df_diseases = readRDS("rds_objects/IPA_pathways/final_objs/sublin_diseases_V1.rds")
df_canonical = readRDS("rds_objects/IPA_pathways/final_objs/sublin_canonical_pathways_V1.rds")


df_diseases$Database = "disease"
df_canonical$Database = "canonical"

df_merged = rbind(df_diseases,df_canonical)

## Fix the main lin nomenclature

df_merged$Main_lin_fixed = df_merged$Main_lin

df_merged$Main_lin_fixed[df_merged$Main_lin_fixed=="Troph"] = "Trophoblasts"
df_merged$Main_lin_fixed[df_merged$Main_lin_fixed=="HBC"] = "HBC"
df_merged$Main_lin_fixed[df_merged$Main_lin_fixed=="Monocytes"] = "Myeloid"
df_merged$Main_lin_fixed[df_merged$Main_lin_fixed=="DC"] = "Myeloid"
df_merged$Main_lin_fixed[df_merged$Main_lin_fixed=="EndoF"] = "Endothelial"
df_merged$Main_lin_fixed[df_merged$Main_lin_fixed=="EndoM"] = "Endothelial"
df_merged$Main_lin_fixed[df_merged$Main_lin_fixed=="PAMM"] = "Myeloid"
df_merged$Main_lin_fixed[df_merged$Main_lin_fixed=="dS"] = "Stromal"
df_merged$Main_lin_fixed[df_merged$Main_lin_fixed=="fFB"] = "Stromal"
df_merged$Main_lin_fixed[df_merged$Main_lin_fixed=="fPV"] = "Stromal"

## fix the sub_lin names

df_merged$Cluster_fixed <- substr(df_merged$Cluster, 1, nchar(df_merged$Cluster) - 3)

df_merged[grepl("HBC_2_CXCL2_CXCL3",df_merged$Cluster_fixed),]$Cluster_fixed ="HBC_2"
df_merged[grepl("HBC_10_cycling",df_merged$Cluster_fixed),]$Cluster_fixed ="HBC_10"
df_merged[grepl("VCTp_1",df_merged$Cluster_fixed),]$Cluster_fixed ="VCT_cycling_1"
df_merged[grepl("VCTp_2",df_merged$Cluster_fixed),]$Cluster_fixed ="VCT_cycling_2"
df_merged[grepl("HBC_4_IL1B_MMP19",df_merged$Cluster_fixed),]$Cluster_fixed ="HBC_4"
df_merged[grepl("HBC_7_lnRNA",df_merged$Cluster_fixed),]$Cluster_fixed ="HBC_7"
df_merged[grepl("HBC_8",df_merged$Cluster_fixed),]$Cluster_fixed ="HBC_8"
df_merged[grepl("HBC_9_CXCL10_IFIT44",df_merged$Cluster_fixed),]$Cluster_fixed ="HBC_9"
df_merged[grepl("mixed_HBC_13_S100A8",df_merged$Cluster_fixed),]$Cluster_fixed ="HBC_13"

df_merged[grepl("Endo_f1",df_merged$Cluster_fixed),]$Cluster_fixed ="fEndo_1"
df_merged[grepl("Endo_f2",df_merged$Cluster_fixed),]$Cluster_fixed ="fEndo_2"
df_merged[grepl("Endo_f3",df_merged$Cluster_fixed),]$Cluster_fixed ="fEndo_3"
df_merged[grepl("Endo_m1",df_merged$Cluster_fixed),]$Cluster_fixed ="mEndo1"

df_merged[grepl("dS1",df_merged$Cluster_fixed),]$Cluster_fixed ="dS_1"
df_merged[grepl("dS3",df_merged$Cluster_fixed),]$Cluster_fixed ="dS_3"

df_merged[grepl("fFB1",df_merged$Cluster_fixed),]$Cluster_fixed ="fFB_1"
df_merged[grepl("fPV2",df_merged$Cluster_fixed),]$Cluster_fixed ="fPV_2"
df_merged[grepl("fPV3",df_merged$Cluster_fixed),]$Cluster_fixed ="fPV_3"

# SupFig 9D: DM1 and FGR

plot_heatmap_single_condition_mainlin(
  df_comp_merged = df_merged,
  cond_int = "DM1",
  control_int = c("Cont_unlabored"),
  padj_or_pval_ast = "padj",
  filter_all = "N",
  zscore_disease = 2,
  padj_disease = 0.05,
  zscore_canonical = 2,
  padj_canonical = 0.25,
  pathway_annotation_file = "input_excel_sheets/misc_excel/PW_annotation_groups_V1.xlsx",
  filter_annotation_group = c("None"),
  omit_unannotated_pathways = TRUE,
  filter_pathway = c("Chronic Myeloid Leukemia Signaling","Tuberculosis Active Signaling Pathway"),
  save_pdf = TRUE,
  pdf_width = 8,
  pdf_height=15,
  output_prefix = "SupFig6D_DM1_vs_Cont_unlabored_V1.pdf",
  output_dir = "final_plots/SupFig6_DM1_FGR"
)




## Section 6E- FGR vs Cont Labored and Unlabored
plot_heatmap_single_condition_mainlin(
  df_comp_merged = df_merged,
  cond_int = "FGR",
  control_int = c("Cont_unlabored",'Cont_labored'),
  padj_or_pval_ast = "padj",
  filter_all = "N",
  zscore_disease = 2,
  padj_disease = 0.001,
  zscore_canonical = 2,
  padj_canonical = 0.25,
  pathway_annotation_file = "input_excel_sheets/misc_excel/PW_annotation_groups_V1.xlsx",
  filter_annotation_group = c("None"),
  omit_unannotated_pathways = TRUE,
  filter_pathway = c("Chronic Myeloid Leukemia Signaling","Tuberculosis Active Signaling Pathway", "Growth of tumor"),
  save_pdf = TRUE,
  pdf_width = 8,
  pdf_height=15,
  output_prefix = "SupFig6E_FGR_vs_Cont_unlabored_labored_V2.pdf",
  output_dir = "final_plots/SupFig6_DM1_FGR"
)


# ======================================================================
# Part 2: pathway enrichment heatmaps for DM1 and FGR
# ======================================================================

library(ggplot2)
library(reshape2)
library(xlsx)
library(ComplexHeatmap)
library(tidyr)
library(writexl)
library(stringr)
library(dplyr)

code_path="/projects/Pregnancy_placenta/Code_R/Figures/"
code_path_PAL="/projects/Pregnancy_placenta/Code_PAL/Preg_placenta_R_codes/Figures"
source(paste0(code_dir, "Functions/DE_DA_functions.R"))
setwd("/projects/Pregnancy_placenta/Results/22_Final_figures_PAL/SupFig_DM1_FGR/")



#========================================================
# modify PW Files
#========================================================
pw_annot_file_PE <- read_excel(paste0(code_path,"/AUX_files/PE_merged_PW_annotation_v6.xlsx")) %>%
  select(Pathway, Annotation) %>%
  mutate(Annotation = na_if(Annotation, "NA"))
pw_annot_file_PE = pw_annot_file_PE[,c("Pathway","Annotation")]
pw_annot_file_PE$Annotation=="NA"

pw_annot_file_PTB = read_excel(paste0(code_path,"/AUX_files/PTB_merged_Pathway_annotation_v5.xlsx")) %>%
  select(Pathway, Annotation) %>%
  mutate(Annotation = na_if(Annotation, "NA"))
pw_annot_file_PTB = pw_annot_file_PTB[,c("Pathway","Annotation")]

# Rename annotation columns to track origin
pw_PE  <- pw_annot_file_PE  %>% rename(Annotation_PE  = Annotation)
pw_PTB <- pw_annot_file_PTB %>% rename(Annotation_PTB = Annotation)

# Full join on Pathway
merged_pw <- full_join(pw_PE, pw_PTB, by = "Pathway") %>%
  arrange(Pathway)

## LOADING MY OWN MERGED DM1 results
merged_pw = read_excel(paste0(code_path_PAL,"/AUX_files/DM1_merged_PW_annotation_v3.xlsx"))%>%
  select(Pathway, Annotation) %>%
  mutate(Annotation = na_if(Annotation, "NA"))
merged_pw= merged_pw[,c("Pathway","Annotation")]

#========================================================
# SupFig 9 (panel TBD) - pathway enrichment in DM1
#========================================================

df_canonical = read.table(paste0(code_path, "AUX_files/IPA_res_sublin_canonical_pathways.txt"), sep="\t", header =1)
df_disease = read.table(paste0(code_path, "AUX_files/IPA_res_sublin_diseases.txt"), sep="\t", header =1)

df = rbind(cbind(df_canonical, Type="Canonical_Pathways"), cbind(df_disease, Type="Diseases_and_biofunctions"))
DM1_df <- df[df$Comparison %in% "DM1_vs_Cont_unlabored",]
DM1_df[DM1_df["Lineage"]=="HBC","Lineage"] <- "Myeloid"


# get unique DM1 pathways
dm1_paths <- unique(DM1_df$Canonical_Pathways)

# add in_DM1 column to merged_pw
merged_pw <- merged_pw %>%
  mutate(in_DM1 = ifelse(Pathway %in% dm1_paths, "Y", "N"))

# find DM1 pathways not already present in merged_pw
missing_dm1 <- setdiff(dm1_paths, merged_pw$Pathway)

# make rows to append
missing_rows <- data.frame(
  Pathway = missing_dm1,
  #Annotation_PE = NA,
  #Annotation_PTB = NA,
  Annotation = NA,
  in_DM1 = "Y",
  stringsAsFactors = FALSE
)

# append missing rows
merged_pw <- bind_rows(merged_pw, missing_rows) %>%
  arrange(Pathway)

# merge columns, prioritie PE vs PTB

## ignoring this for now as we want to see which pathways appear but not annoted
#merged_pw_sub <- merged_pw %>%
#  filter(
#    in_DM1 == "Y",
#    !(is.na(Annotation_PE) & is.na(Annotation_PTB))
#  )
merged_pw_sub = merged_pw

#merged_pw_sub <- merged_pw_sub %>%
#  mutate(
    #Annotation = coalesce(Annotation_PE, Annotation_PTB)
#  )

#adding this to see whats not annotated
merged_pw_sub[is.na(merged_pw_sub$Annotation),]$Annotation = "Not_Annotated"

merged_pw_sub = merged_pw_sub[,c("Pathway","Annotation")]


write_xlsx(merged_pw_sub,paste0(code_path_PAL,"/AUX_files/DM1_merged_PW_annotation_v2_PLUS_NOT_ASSIGNED.xlsx") )

# ----- RUN pathway heatmap -----
plot_heatmap_PW(df_comp_merged=DM1_df,
                        zscore_min_cp =1,
                        padj_thresh_cp = 0.25,
                        zscore_min_db = 2,
                        padj_thresh_db = 0.05,
                        save_pdf = TRUE,
                        #output_prefix = paste0("SupFig_DM1_both_adjPval_0.25_v1"),
                output_prefix = paste0("SupFig_DM1_both_adjPval_0.05_ALL_PATHS_V5_AE"),
                        orderColBy = "Cells",
                       # order_vector=c("HBC_2_TermPE","HBC_3_TermPE","CD14_monocyte_TermPE","HBC_4_TermPE","HBC_13_TermPE","HBC_2_PrePE","HBC_3_PrePE","CD14_monocyte_PrePE","HBC_1_PrePE","SCT_TermPE","VCT_1_TermPE","VCT_3_TermPE","VCTfusion_1_TermPE","VCTfusion_2_TermPE","VCT_cycling_1_TermPE","VCT_cycling_2_TermPE","SCT_PrePE","VCT_1_PrePE","VCT_3_PrePE","VCTfusion_1_PrePE","VCTfusion_2_PrePE","VCTccc_PrePE","fEndo_1_PrePE","fEndo_2_PrePE","fEndo_3_PrePE"),
                        pdf_width = 15,
                      split_cols = c(DM1="black"),
                        pdf_height = 28,
                        cluster_names_col="Cluster",
                        filter_more="N",
                row_display_mode = "pathway_annotation_group",
                        #plot_title="PE vs unlabored_cont",
                       pw_annot_file = paste0(code_path_PAL,"/AUX_files/DM1_merged_PW_annotation_v4.xlsx"))


#========================================================
# SupFig 9 (panel TBD) - pathway enrichment in FGR
#========================================================
merged_pw = read_excel(paste0(code_path_PAL,"/AUX_files/DM1_merged_PW_annotation_v3.xlsx"))%>%
  select(Pathway, Annotation) %>%
  mutate(Annotation = na_if(Annotation, "NA"))
merged_pw= merged_pw[,c("Pathway","Annotation")]

df_canonical = read.table(paste0(code_path, "AUX_files/IPA_res_sublin_canonical_pathways.txt"), sep="\t", header =1)
df_disease = read.table(paste0(code_path, "AUX_files/IPA_res_sublin_diseases.txt"), sep="\t", header =1)

df = rbind(cbind(df_canonical, Type="Canonical_Pathways"), cbind(df_disease, Type="Diseases_and_biofunctions"))
FGR_df <- df[df$Comparison %in% c("FGR_vs_Cont_unlabored","FGR_vs_Cont_labored"),]
FGR_df [FGR_df ["Lineage"]=="HBC","Lineage"] <- "Myeloid"


# get unique DM1 pathways
fgr_paths <- unique(FGR_df$Canonical_Pathways)

# add in_DM1 column to merged_pw
merged_pw_FGR <- merged_pw %>%
  mutate(in_FGR = ifelse(Pathway %in% fgr_paths, "Y", "N"))

# find FGR pathways not already present in merged_pw
missing_fgr <- setdiff(fgr_paths, merged_pw_FGR$Pathway)

# make rows to append
missing_rows <- data.frame(
  Pathway = missing_fgr,
  #Annotation_PE = NA,
  #Annotation_PTB = NA,
  in_FGR= "Y",
  Annotation = NA,
  stringsAsFactors = FALSE
)

# append missing rows
merged_pw_FGR <- bind_rows(merged_pw_FGR, missing_rows) %>%
  arrange(Pathway)

# merge columns, prioritie PE vs PTB
#merged_pw_FGR_sub <- merged_pw_FGR %>%
#  filter(
#    in_FGR == "Y",
#    !(is.na(Annotation_PE) & is.na(Annotation_PTB))
#  )
merged_pw_FGR_sub = merged_pw_FGR
#merged_pw_FGR_sub <- merged_pw_FGR_sub %>%
#  mutate(
#    Annotation = coalesce(Annotation_PE, Annotation_PTB)
#  )

#merged_pw_FGR_sub = merged_pw_FGR_sub[,c("Pathway","Annotation")]
merged_pw_FGR_sub[is.na(merged_pw_FGR_sub$Annotation),]$Annotation = "Not_Annotated"

merged_pw_FGR_sub = merged_pw_FGR_sub[,c("Pathway","Annotation")]

write_xlsx(merged_pw_sub,paste0(code_path_PAL,"/AUX_files/FGR_merged_PW_annotation_v3.xlsx") )


# ----- RUN pathway heatmap -----
plot_heatmap_PW(df_comp_merged=FGR_df,
                        zscore_min_cp = 2,
                        padj_thresh_cp = 0.25,
                        zscore_min_db = 3,
                        padj_thresh_db = 0.001,
                        save_pdf = TRUE,
                        #output_prefix = paste0("SupFig_FGR_CP_adjPval_0.05_DB_adjPval_0.001_test4"),
                      output_prefix = paste0("SupFig_FGR_SMFM_thresh_V1"),
                        orderColBy = "Cells",
                       # order_vector=c("HBC_2_TermPE","HBC_3_TermPE","CD14_monocyte_TermPE","HBC_4_TermPE","HBC_13_TermPE","HBC_2_PrePE","HBC_3_PrePE","CD14_monocyte_PrePE","HBC_1_PrePE","SCT_TermPE","VCT_1_TermPE","VCT_3_TermPE","VCTfusion_1_TermPE","VCTfusion_2_TermPE","VCT_cycling_1_TermPE","VCT_cycling_2_TermPE","SCT_PrePE","VCT_1_PrePE","VCT_3_PrePE","VCTfusion_1_PrePE","VCTfusion_2_PrePE","VCTccc_PrePE","fEndo_1_PrePE","fEndo_2_PrePE","fEndo_3_PrePE"),
                        pdf_width = 15,
                      split_cols = c(Cont_labored="Grey", Cont_unlabored="black"),
                        pdf_height = 28,
                        cluster_names_col="Cluster",
                        filter_more="N",
                row_display_mode = "pathway_annotation_group",
                split_by_cond_OR_ctrl = "Ctrl",
                        #plot_title="PE vs unlabored_cont",
                       pw_annot_file = paste0(code_path_PAL,"/AUX_files/FGR_merged_PW_annotation_v3.xlsx"))

