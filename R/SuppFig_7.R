# Supplemental Figure 7: other control analyses

# Section 1: Load helper functions

code_dir="/projects/Pregnancy_placenta/Code_Publication/R/"
source(paste0(code_dir, "Functions/funcs_supfig5_supfig6.R"))


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

lin_colors = c(Myeloid = "#590902", HBC = "#BF3604", Trophoblasts = "#60B364", Stromal = "#022873", Endothelial = "#6B0AC9")

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

# Sup Fig7C- sPTB vs labored and unlabored

plot_heatmap_single_condition_mainlin(
  df_comp_merged = df_merged,
  cond_int = "PTB",
  control_int = c("Cont_unlabored",'Cont_labored'),
  padj_or_pval_ast = "padj",
  filter_all = "N",
  zscore_disease = 2,
  padj_disease = 0.0001,
  zscore_canonical = 2,
  padj_canonical = 0.01,
  pathway_annotation_file = "input_excel_sheets/misc_excel/PW_annotation_groups_V1.xlsx",
  filter_annotation_group = c("None"),
  omit_unannotated_pathways = TRUE,
  filter_pathway = c("Chronic Myeloid Leukemia Signaling","Tuberculosis Active Signaling Pathway", "Growth of tumor","Activation of lymphocytes"),
  save_pdf = TRUE,
  pdf_width = 10,
  pdf_height=15,
  output_prefix = "SupFig5C_sPTB_vs_labored_unlabored_V1",
  output_dir = "final_plots/SupFig5_other_controls/"
)


# SupFig7D- PrePE vs Cont PTB

plot_heatmap_single_control_mainlin(
  df_comp_merged = df_merged,
  control_int = "Cont_PTB",
  cond_int = c("PrePE"),
  lineage_column = "Main_lin_fixed",
  lineage_colors = lin_colors,
  padj_or_pval_ast = "padj",
  filter_all = "N",
   zscore_disease = 2,
  padj_disease = 0.05,
  zscore_canonical = 2,
  padj_canonical = 0.05,
  filter_annotation_group = c("None"),
  omit_unannotated_pathways = TRUE,
  filter_pathway = c("Growth of tumor"),
  pathway_annotation_file = "input_excel_sheets/misc_excel/PW_annotation_groups_V1.xlsx",
  save_pdf = TRUE,
  pdf_width = 10,
  pdf_height=12,
  output_prefix = "SupFig5D_PrePE_vs_Cont_PTB_V1",
  output_dir = "final_plots/SupFig5_other_controls/"

)



# SupFig7E - TermPE and PrePE vs Control LAbored
plot_heatmap_single_control_mainlin(
  df_comp_merged = df_merged,
  control_int = "Cont_labored",
  cond_int = c("PrePE","TermPE"),
  lineage_column = "Main_lin_fixed",
  lineage_colors = lin_colors,
  padj_or_pval_ast = "padj",
  filter_all = "N",
   zscore_disease = 2,
  padj_disease = 0.25,
  zscore_canonical = 2,
  padj_canonical = 0.25,
  filter_annotation_group = c("None"),
  omit_unannotated_pathways = TRUE,
  filter_pathway = c("Growth of tumor"),
  pathway_annotation_file = "input_excel_sheets/misc_excel/PW_annotation_groups_V1.xlsx",
 save_pdf = TRUE,
  pdf_width = 10,
  pdf_height=10,
  output_prefix = "SupFig5E_TermPE_PrePE_vs_Cont_Labored_V1",
  output_dir = "final_plots/SupFig5_other_controls/"
)


# ======================================================================
# Part 2: Venn diagram and shared-pathway Z-score scatter (labored vs unlabored controls)
# ======================================================================

# Purpose of Script ----
# Supp Figure 7: compare pathway results using unlabored and labored controls (Venn diagram and shared-pathway Z-score scatter).
# Run from the FigS6_code folder, either all at once or section by section.
#
# Documentation:
# This script loads the study inputs, matches the control comparisons (unlabored vs labored),
# applies |Z-score| > 2 and adjusted P < 0.25, and saves the Venn diagram, Z-score scatter and
# source tables in `output/`. Results are counted separately for each pathway, cell population
# and complication, across 103 manuscript-selected IPA terms, 35 populations (including SCT and
# VCTccc), and five complications. The scatter requires adjusted P < 0.25 and numeric Z-scores
# against both controls; colored points also pass both Z-score cutoffs. Reusable input checks
# and formatting are in Functions/SuppFig_7_Venn_helper.R.

# Preparing Packages and Functions ----
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readxl)
  library(systemfonts)
})

source(paste0(code_dir, "Functions/SuppFig_7_Venn_helper.R"))

# Load Variables ----
project_dir <- "../../.."   # Folder containing IPA_results and rds_objects
output_dir <- "output"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

padj_cutoff <- 0.25
zscore_cutoff <- 2
font_family <- "Arial"

complications <- c(
  PTB = "sPTB",
  PrePE = "Preterm PE",
  TermPE = "Term PE",
  FGR = "FGR",
  DM1 = "DM1"
)

control_colors <- c(Unlabored = "#CD534C", Labored = "#0073C2")
complication_shapes <- c("sPTB" = 21, "Preterm PE" = 24, "Term PE" = 22,
                        "FGR" = 23, "DM1" = 25)

# Load Pathway Data ----
canonical_file <- file.path(project_dir,
                           "rds_objects/IPA_pathways/final_objs/sublin_canonical_pathways_V1.rds")
disease_file <- file.path(project_dir,
                         "rds_objects/IPA_pathways/final_objs/sublin_diseases_V1.rds")

canonical_pathways <- readRDS(canonical_file)
disease_functions <- readRDS(disease_file)
pathway_inventory <- read.delim("inputs/pathway_inventory.tsv", encoding = "UTF-8")
cell_palette <- read.delim("inputs/cell_colors.tsv", encoding = "UTF-8")
label_positions <- read.delim("inputs/label_positions.tsv", encoding = "UTF-8")

# Keep the 103 terms shown in manuscript figures.
selected_pathways <- pathway_inventory %>%
  filter(grepl("^(Figure|SupFigure)", Source)) %>%
  mutate(Pathway = clean_ipa_names(Pathway),
         Pathway_id = paste(Database, Pathway, sep = "::")) %>%
  group_by(Database, Pathway, Pathway_id) %>%
  summarise(Paper_sources = paste(sort(unique(Source)), collapse = "; "),
            .groups = "drop") %>%
  arrange(Database, Pathway)

# The saved sublineage objects already include parent-source SCT and VCTccc.
ipa_results <- prepare_ipa_results(canonical_pathways, disease_functions)
populations <- unique(ipa_results[c("Cell_type", "Cell_label", "Source_level")])
populations <- populations[order(match(populations$Cell_type,
                                      c("SCT_L1", "VCTccc_L1"))), ]

comparison_headers <- unique(ipa_results[c(
  "Database", "Cell_type", "Cell_label", "Source_level", "Condition",
  "Control", "Comparison", "CompNum", "Clust_Comparison"
)])

ipa_results <- ipa_results %>%
  filter(Pathway_id %in% selected_pathways$Pathway_id,
         Condition %in% names(complications),
         Control %in% c("Cont_unlabored", "Cont_labored"))

# Match each selected statistic to its original TEST3 export before using it.
source_check <- check_ipa_statistics(ipa_results, project_dir)
ipa_results <- source_check$data

# IPA stores -log10(adjusted P); convert it back to adjusted P.
ipa_results$Adjusted_p <- 10^(-ipa_results$log10_padj)
ipa_results$Complication <- unname(complications[ipa_results$Condition])

# Match Unlabored and Labored Results ----
# One row = one pathway, cell population and complication.
unlabored_results <- ipa_results %>%
  filter(Control == "Cont_unlabored") %>%
  select(Pathway_id, Cell_type, Condition,
         Unlabored_Zscore = Zscore, Unlabored_Adjusted_p = Adjusted_p)

labored_results <- ipa_results %>%
  filter(Control == "Cont_labored") %>%
  select(Pathway_id, Cell_type, Condition,
         Labored_Zscore = Zscore, Labored_Adjusted_p = Adjusted_p)

pathway_results <- expand.grid(
  Pathway_id = selected_pathways$Pathway_id,
  Cell_type = populations$Cell_type,
  Condition = names(complications),
  stringsAsFactors = FALSE
) %>%
  left_join(select(selected_pathways, Pathway_id, Pathway, Database), by = "Pathway_id") %>%
  left_join(select(populations, Cell_type, Cell_label), by = "Cell_type") %>%
  mutate(Complication = unname(complications[Condition])) %>%
  left_join(unlabored_results, by = c("Pathway_id", "Cell_type", "Condition")) %>%
  left_join(labored_results, by = c("Pathway_id", "Cell_type", "Condition"))

# Apply the Significance Thresholds ----
# X = |Z| > 2 AND adjusted P < 0.25. Everything else, including NA, is "-".
# The same pathway counts separately in each cell population and complication.
pathway_results <- pathway_results %>%
  mutate(
    Unlabored_flag = if_else(
      is.finite(Unlabored_Zscore) & is.finite(Unlabored_Adjusted_p) &
        abs(Unlabored_Zscore) > zscore_cutoff & Unlabored_Adjusted_p < padj_cutoff,
      "X", "-"
    ),
    Labored_flag = if_else(
      is.finite(Labored_Zscore) & is.finite(Labored_Adjusted_p) &
        abs(Labored_Zscore) > zscore_cutoff & Labored_Adjusted_p < padj_cutoff,
      "X", "-"
    ),
    Category = case_when(
      Unlabored_flag == "X" & Labored_flag == "X" ~ "Both",
      Unlabored_flag == "X" ~ "Unlabored only",
      Labored_flag == "X" ~ "Labored only",
      TRUE ~ "Neither"
    ),
    Direction = case_when(
      Category != "Both" ~ NA_character_,
      sign(Unlabored_Zscore) == sign(Labored_Zscore) ~ "Same sign",
      TRUE ~ "Opposite signs"
    )
  )

# Count Pathway Results ----
venn_counts <- pathway_results %>%
  summarise(Unlabored_only = sum(Category == "Unlabored only"),
            Both = sum(Category == "Both"),
            Labored_only = sum(Category == "Labored only"),
            Total_unlabored = sum(Unlabored_flag == "X"),
            Total_labored = sum(Labored_flag == "X"))

counts_by_complication <- pathway_results %>%
  group_by(Complication) %>%
  summarise(Unlabored_only = sum(Category == "Unlabored only"),
            Both = sum(Category == "Both"),
            Labored_only = sum(Category == "Labored only"),
            Total = sum(Category != "Neither"), .groups = "drop") %>%
  arrange(match(Complication, complications))

counts_by_complication <- bind_rows(
  counts_by_complication,
  data.frame(Complication = "All complications",
             as.list(colSums(counts_by_complication[-1])))
)
print(venn_counts)

# Prepare the Z-score Scatter ----
# All plotted points have numeric Z and adjusted P < 0.25 against BOTH controls.
# Color the shared qualifying results; show the remaining points in gray.
scatter_results <- pathway_results %>%
  filter(is.finite(Unlabored_Zscore), is.finite(Labored_Zscore),
         is.finite(Unlabored_Adjusted_p), is.finite(Labored_Adjusted_p),
         Unlabored_Adjusted_p < padj_cutoff, Labored_Adjusted_p < padj_cutoff)

shared_results <- scatter_results %>%
  filter(Category == "Both") %>%
  arrange(match(Condition, names(complications)), Cell_label, Pathway)

gray_results <- scatter_results %>%
  filter(Category != "Both")

# Label every shared result. Moving a label never changes the point's Z-score.
plot_labels <- shared_results %>%
  left_join(label_positions, by = c("Pathway_id", "Cell_type", "Condition"))

cell_colors <- setNames(cell_palette$Subset_color, cell_palette$Cell_label)
check_plot_inputs(scatter_results, plot_labels, cell_colors, font_family)
label_lines <- make_label_lines(plot_labels, font_family)

# Plot the Venn Diagram ----
# The helper draws equal-size circles and places the calculated counts inside.
venn_plot <- plot_pathway_venn(
  counts = venn_counts,
  colors = control_colors,
  title = "Pathway overlap across control groups",
  family = font_family
)

# Plot the Shared Pathway Z-scores ----
scatter_plot <- ggplot(scatter_results,
                       aes(x = Unlabored_Zscore, y = Labored_Zscore)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#A0A0A0", linewidth = 0.14) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#A0A0A0", linewidth = 0.14) +
  geom_point(data = gray_results, color = "#D8D8D8", alpha = 0.7, size = 0.8) +
  geom_segment(data = label_lines,
               aes(x = Start_x, y = Start_y, xend = End_x, yend = End_y),
               inherit.aes = FALSE, color = "#505050", linewidth = 0.14) +
  geom_point(data = shared_results, aes(fill = Cell_label, shape = Complication),
             color = "#303030", size = 2, stroke = 0.2) +
  geom_text(data = label_lines, aes(x = Label_x, y = Label_y, label = Label),
            inherit.aes = FALSE, family = font_family, size = 8 / .pt,
            lineheight = 1.2) +
  scale_fill_manual(values = cell_colors, breaks = names(cell_colors)) +
  scale_shape_manual(values = complication_shapes, breaks = names(complication_shapes)) +
  scale_x_continuous(breaks = c(-4, -2, 0, 2, 4)) +
  scale_y_continuous(breaks = c(-4, -2, 0, 2, 4)) +
  coord_fixed(xlim = c(-5.5, 5.5), ylim = c(-5.5, 5.5), expand = FALSE) +
  labs(title = "Shared pathway directionality across control groups",
       x = "Complication vs Unlabored Control (Z-score)",
       y = "Complication vs Labored Control (Z-score)",
       fill = "Cell population", shape = "Complication") +
  theme_classic(base_family = font_family, base_size = 8) +
  theme(axis.line = element_blank(),
        panel.border = element_rect(color = "#333333", fill = NA, linewidth = 0.18),
        axis.title = element_text(size = 8.5),
        axis.text = element_text(size = 8, color = "black"),
        axis.ticks = element_line(linewidth = 0.18),
        axis.ticks.length = unit(3, "pt"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 10.5),
        legend.position = "bottom",
        legend.box = "vertical",
        legend.box.just = "center",
        legend.spacing.y = unit(0, "pt"),
        legend.margin = margin(0, 0, 0, 0),
        legend.key.height = unit(12, "pt"),
        legend.key.width = unit(12, "pt"),
        legend.title = element_text(face = "bold", size = 8),
        legend.text = element_text(size = 8)) +
  guides(fill = guide_legend(order = 1, nrow = 1, override.aes = list(shape = 21)),
         shape = guide_legend(order = 2, nrow = 1, override.aes = list(fill = "#AAAAAA")))

# Save Figures ----
ggsave(venn_plot,
       filename = file.path(output_dir, "FigS6_pathway_overlap_Venn.pdf"),
       device = cairo_pdf, family = font_family,
       width = 280 / 72, height = 246 / 72, units = "in")

# Keep the square plot, title and bottom legend rows centered on the PDF page.
save_scatter_pdf(scatter_plot,
                 filename = file.path(output_dir, "FigS6_shared_pathway_scatter.pdf"),
                 family = font_family)

# Save Source Data and Summary Tables ----
write_tsv(pathway_results, file.path(output_dir, "Pathway_cell_complication_results.tsv"))
write_tsv(venn_counts, file.path(output_dir, "Venn_counts.tsv"))
write_tsv(shared_results, file.path(output_dir, "Shared_results.tsv"))
write_tsv(gray_results, file.path(output_dir, "Grey_background_results.tsv"))
write_tsv(counts_by_complication, file.path(output_dir, "Counts_by_complication.tsv"))
write_tsv(selected_pathways, file.path(output_dir, "Pathway_universe.tsv"))
write_tsv(populations, file.path(output_dir, "Included_populations.tsv"))
write_tsv(ipa_results, file.path(output_dir, "Source_results_long.tsv"))
write_tsv(comparison_headers, file.path(output_dir, "Original_comparison_headers.tsv"))
write_tsv(plot_labels, file.path(output_dir, "Named_results.tsv"))
write_tsv(label_lines, file.path(output_dir, "Scatter_label_layout.tsv"))
write_tsv(cell_palette, file.path(output_dir, "Cell_colors.tsv"))
write_tsv(data.frame(Complication = names(complication_shapes), Shape = complication_shapes),
          file.path(output_dir, "Complication_shapes.tsv"))

save_complication_tables(pathway_results, selected_pathways, populations,
                         complications, output_dir)

# Save Source Checks ----
write_tsv(source_check$checks, file.path(output_dir, "Raw_RDS_verification.tsv"))
input_files <- c(canonical_file, disease_file, source_check$files,
                 "inputs/pathway_inventory.tsv", "inputs/cell_colors.tsv",
                 "inputs/label_positions.tsv", file.path(code_dir, "SuppFig_7.R"),
                 file.path(code_dir, "Functions/SuppFig_7_Venn_helper.R"))
write_tsv(data.frame(File = input_files, MD5 = unname(tools::md5sum(input_files))),
          file.path(output_dir, "Input_checksums.tsv"))
write_tsv(data.frame(Adjusted_p_cutoff = padj_cutoff, Absolute_Z_cutoff = zscore_cutoff,
                     Pathways = nrow(selected_pathways), Populations = nrow(populations),
                     Complications = length(complications), Shared_results = nrow(shared_results)),
          file.path(output_dir, "Analysis_settings.tsv"))
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))

message("Figure S6 PDFs and tables saved in: ", normalizePath(output_dir))
