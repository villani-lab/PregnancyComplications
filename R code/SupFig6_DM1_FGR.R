# Supplemental Figure 6: DM1 and FGR

# Section 1: Load helper functions

source("helpful_functions/funcs_supfig5_supfig6.R")


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

# SupFig 6D: DM1 and FGR

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
