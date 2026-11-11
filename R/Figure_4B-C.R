# Figure 4B) PTB pathway enrichment heatmap
# Heatmap of IPA pathways enrichemnt results in PTB vs cont_PTB across cell types

code_dir="/projects/Pregnancy_placenta/Code_Publication/R/"
source(paste0(code_dir, "Functions/figure_functions.R"))
code_path = "/projects/Pregnancy_placenta/Code_R/Figures/"
setwd("/projects/Pregnancy_placenta/Results/22_Final_figures/Figure_4/")

# Merged - canonical pathways and diseases and biofunctions
df_canonical = read.table(paste0(code_path, "/AUX_files/IPA_res_sublin_canonical_pathways.txt"), sep="\t", header =1)
df_disease = read.table(paste0(code_path, "/AUX_files/IPA_res_sublin_diseases.txt"), sep="\t", header =1)
df = rbind(cbind(df_canonical, Type="Canonical_Pathways"), cbind(df_disease, Type="Diseases_and_biofunctions"))
PTB_df = df[(df$Comparison == "PTB_vs_Cont_PTB"),]

plot_heatmap_PW(df_comp_merged=PTB_df,
                        zscore_min_cp = 2,
                        padj_thresh_cp = 0.01,
                        zscore_min_db = 3,
                        padj_thresh_db = 0.0001,
                        save_pdf = TRUE,
                        output_prefix = paste0("Fig4B_PTB_PW_heatmap.pdf"),
                        order_vector=c("HBC_1","HBC_10","HBC_11","HBC_13","HBC_2","HBC_4","HBC_6","HBC_9","CD14_monocyte","CD16_monocyte","cDC","mMac_upper","fPV2","fPV3","mEndo_1","VCT_cycling_1","VCT_cycling_2","VCTfusion_1","VCT_1"),
                        orderColBy = "Cells",
                        pdf_width = 15,
                        pdf_height = 24,
                        cluster_names_col="Cluster",
                        filter_more="N",
                        plot_title="PTB vs PTB_cont",
                        pw_annot_file = paste0(code_path, "/AUX_files/PTB_merged_Pathway_annotation_v5.xlsx"))


# ======================================================================
# Figure 4C
# ======================================================================

# Figure 4C) PTB gene dotplots

# Load packages ----
suppressPackageStartupMessages({
  library(openxlsx)
  library(readxl)
  library(tidyr)
  library(stringr)
  library(ComplexHeatmap)
  library(colorRamp2)
  library(dplyr)
  library(circlize)
  library(ComplexUpset)
  library(grid) # for adjusting matrix box size
  })




#### load functions ----
#setwd("/projects/Pregnancy_placenta/Results/13_DE_PAL")
source(paste0(code_dir, "Functions/funcs_gene_pathway.R"))


# load directory of files ----
setwd("/projects/Pregnancy_placenta/Results/13_DE_PAL")

# load file going over where things are located ----
parent_lin = read.csv("/projects/Pregnancy_placenta/Results/13_DE_PAL/output_csv_IPA/redo_Parentlin/Parent_Test_3.csv")
sub_lin = read.csv("/projects/Pregnancy_placenta/Results/13_DE_PAL/output_csv_IPA/redo_sublin/sublin_Test_3.csv")

all_lin_loc = rbind(data.frame(parent_lin),data.frame(sub_lin))

# Load directory of variables----
IPA_pathway_genes_directory = "/projects/Pregnancy_placenta/Results/13_DE_PAL/output_csv_IPA/output_IPA_genes"

#### Workspace----
list_files = list.files(path = "/projects/Pregnancy_placenta/Results/20_Figures/Pathways/IPA_genes/")


list_files = list.files("/projects/Pregnancy_placenta/Results/20_Figures/Pathways/IPA_genes/")


# fPV2 vs fPV3 - Apoptosis and Necrosis
mats = get_pathway_gene_heatmap_mats(
  path_of_interest = "Necrosis",
  list_files = list_files,
  df_loc = all_lin_loc,
  file_regex = "PTB",
  loc_of_IPA_gene_files = "/projects/Pregnancy_placenta/Results/20_Figures/Pathways/IPA_genes/",
  cell_types_keep = c("fPV2",'fPV3'),
  comparison_regex = "Cont_PTB"
)

ht = plot_pathway_gene_heatmap(
  log2fc_mat = mats$log2fc_mat,
  ast_mat = mats$ast_mat,
  path_of_interest = "Necrosis",
  column_title = "Necrosis - fPV2 vs fPV3 in PTB vs ContPTB",
  annotation_values = ifelse(grepl("TermPE", colnames(mats$log2fc_mat)), "TermPE", "PrePE"),
  annotation_colors = c("TermPE" = "#9ecae1", "PrePE" = "#fdae6b"),
  column_labels = sub("-.*", "", colnames(mats$log2fc_mat))
)


# add comparison
apoptosis_df <- mats$pathway_gene_df
apoptosis_df$comparison<-paste0("Apoptosis-",apoptosis_df$comparison)

necrosis_df<-mats$pathway_gene_df
necrosis_df$comparison<-paste0("Necrosis-",necrosis_df$comparison)

pathway_info_merged <- rbind(apoptosis_df,necrosis_df)

## draw Upset Plotf
upset_input_df <- make_complex_upset_input_genes(
  pathway_info_merged,
  as_logical = FALSE
)
comparisons = unique(pathway_info_merged$comparison)
upset(upset_input_df,comparisons,
      set_sizes=(
        upset_set_size()
        + geom_text(aes(label=..count..), hjust=1.1, stat='count')))

## save unique comparison
csv_path_int = "/projects/Pregnancy_placenta/Code_PAL/Preg_placenta_R_codes/Figures/pathway_csv_ht_results/Fig4_PTB_vs_Cont_PTB_apoptosis_necrosis_fPV2_fPV3.csv"
write.csv(pathway_info_merged ,file = csv_path_int)

make_gene_overlap_summary_csv(
  csv_file = csv_path_int,
  output_dir ="/projects/Pregnancy_placenta/Code_PAL/Preg_placenta_R_codes/Figures/upset_plot_results/figure4",
  output_prefix = "fig4_necrosis_apoptosis_comp_fPV2_fPV3",
   save_long_csv = FALSE,
)


# fPV2 vs fPV3 - Angiogenesis
mats = get_pathway_gene_heatmap_mats(
  path_of_interest = "Angiogenesis",
  list_files = list_files,
  df_loc = all_lin_loc,
  file_regex = "PTB",
  loc_of_IPA_gene_files = "/projects/Pregnancy_placenta/Results/20_Figures/Pathways/IPA_genes/",
  cell_types_keep = c("fPV2",'fPV3'),
  comparison_regex = "Cont_PTB"
)

ht = plot_pathway_gene_heatmap(
  log2fc_mat = mats$log2fc_mat,
  ast_mat = mats$ast_mat,
  path_of_interest = "Angiogenesis",
  column_title = "Angiogenesis - fPV2 vs fPV3 in PTB vs ContPTB",
  annotation_values = ifelse(grepl("TermPE", colnames(mats$log2fc_mat)), "TermPE", "PrePE"),
  annotation_colors = c("TermPE" = "#9ecae1", "PrePE" = "#fdae6b"),
  column_labels = sub("-.*", "", colnames(mats$log2fc_mat))
)


## draw Upset Plotf
upset_input_df <- make_complex_upset_input_genes(
  mats$pathway_gene_df,
  as_logical = FALSE
)
comparisons = unique(mats$pathway_gene_df$comparison)
upset(upset_input_df,comparisons,
      set_sizes=(
        upset_set_size()
        + geom_text(aes(label=..count..), hjust=1.1, stat='count')))

## save uniue comparison
csv_path_int = "/projects/Pregnancy_placenta/Code_PAL/Preg_placenta_R_codes/Figures/pathway_csv_ht_results/Fig4_PTB_vs_Cont_PTB_angiogenesis_fPV2_fPV3.csv"
write.csv(mats$pathway_gene_df ,file = csv_path_int)

make_gene_overlap_summary_csv(
  csv_file = csv_path_int,
  output_dir ="/projects/Pregnancy_placenta/Code_PAL/Preg_placenta_R_codes/Figures/upset_plot_results/figure4",
  output_prefix = "fig4_PTB_ContPTB_angiogenesis_fPV2_fPV3",
   save_long_csv = FALSE,
)


# HBC - Angiogenesis
mats = get_pathway_gene_heatmap_mats(
  path_of_interest = "Angiogenesis",
  list_files = list_files,
  df_loc = all_lin_loc,
  file_regex = "PTB",
  loc_of_IPA_gene_files = "/projects/Pregnancy_placenta/Results/20_Figures/Pathways/IPA_genes/",
  cell_types_keep = c("HBC"),
  comparison_regex = "Cont_PTB"
)

ht = plot_pathway_gene_heatmap(
  log2fc_mat = mats$log2fc_mat,
  ast_mat = mats$ast_mat,
  path_of_interest = "Angiogenesis",
  column_title = "Angiogenesis - HBC in PTB vs ContPTB",
  annotation_values = ifelse(grepl("TermPE", colnames(mats$log2fc_mat)), "TermPE", "PrePE"),
  annotation_colors = c("TermPE" = "#9ecae1", "PrePE" = "#fdae6b"),
  column_labels = sub("-.*", "", colnames(mats$log2fc_mat))
)


## draw Upset Plotf
upset_input_df <- make_complex_upset_input_genes(
  mats$pathway_gene_df,
  as_logical = FALSE
)
comparisons = unique(mats$pathway_gene_df$comparison)
upset(upset_input_df,comparisons,
      set_sizes=(
        upset_set_size()
        + geom_text(aes(label=..count..), hjust=1.1, stat='count')))

## save uniue comparison
csv_path_int = "/projects/Pregnancy_placenta/Code_PAL/Preg_placenta_R_codes/Figures/pathway_csv_ht_results/Fig4_PTB_vs_Cont_PTB_angiogenesis_HBC.csv"
write.csv(mats$pathway_gene_df ,file = csv_path_int)

make_gene_overlap_summary_csv(
  csv_file = csv_path_int,
  output_dir ="/projects/Pregnancy_placenta/Code_PAL/Preg_placenta_R_codes/Figures/upset_plot_results/figure4",
  output_prefix = "fig4_PTB_ContPTB_angiogenesis_HBC",
   save_long_csv = FALSE,
)


# HBC - migration of cells
mats = get_pathway_gene_heatmap_mats(
  path_of_interest = "Migration of cells",
  list_files = list_files,
  df_loc = all_lin_loc,
  file_regex = "PTB",
  loc_of_IPA_gene_files = "/projects/Pregnancy_placenta/Results/20_Figures/Pathways/IPA_genes/",
  cell_types_keep = c("HBC"),
  comparison_regex = "Cont_PTB"
)

ht = plot_pathway_gene_heatmap(
  log2fc_mat = mats$log2fc_mat,
  ast_mat = mats$ast_mat,
  path_of_interest = "Migration of cells",
  column_title = "Migration of cells - HBC in PTB vs ContPTB",
  annotation_values = ifelse(grepl("TermPE", colnames(mats$log2fc_mat)), "TermPE", "PrePE"),
  annotation_colors = c("TermPE" = "#9ecae1", "PrePE" = "#fdae6b"),
  column_labels = sub("-.*", "", colnames(mats$log2fc_mat))
)


## draw Upset Plotf
upset_input_df <- make_complex_upset_input_genes(
  mats$pathway_gene_df,
  as_logical = FALSE
)
comparisons = unique(mats$pathway_gene_df$comparison)
upset(upset_input_df,comparisons,
      set_sizes=(
        upset_set_size()
        + geom_text(aes(label=..count..), hjust=1.1, stat='count')))

## save uniue comparison
csv_path_int = "/projects/Pregnancy_placenta/Code_PAL/Preg_placenta_R_codes/Figures/pathway_csv_ht_results/Fig4_PTB_vs_Cont_PTB_migration_of_cells_HBC.csv"
write.csv(mats$pathway_gene_df ,file = csv_path_int)

make_gene_overlap_summary_csv(
  csv_file = csv_path_int,
  output_dir ="/projects/Pregnancy_placenta/Code_PAL/Preg_placenta_R_codes/Figures/upset_plot_results/figure4",
  output_prefix = "fig4_PTB_ContPTB_migration_of_cells_HBC",
   save_long_csv = FALSE,
)

# HBC - Cell Movement
mats = get_pathway_gene_heatmap_mats(
  path_of_interest = "Cell movement",
  list_files = list_files,
  df_loc = all_lin_loc,
  file_regex = "PTB",
  loc_of_IPA_gene_files = "/projects/Pregnancy_placenta/Results/20_Figures/Pathways/IPA_genes/",
  cell_types_keep = c("HBC"),
  comparison_regex = "Cont_PTB"
)

ht = plot_pathway_gene_heatmap(
  log2fc_mat = mats$log2fc_mat,
  ast_mat = mats$ast_mat,
  path_of_interest = "Cell movement",
  column_title = "Cell movement - HBC in PTB vs ContPTB",
  annotation_values = ifelse(grepl("TermPE", colnames(mats$log2fc_mat)), "TermPE", "PrePE"),
  annotation_colors = c("TermPE" = "#9ecae1", "PrePE" = "#fdae6b"),
  column_labels = sub("-.*", "", colnames(mats$log2fc_mat))
)


## draw Upset Plotf
upset_input_df <- make_complex_upset_input_genes(
  mats$pathway_gene_df,
  as_logical = FALSE
)
comparisons = unique(mats$pathway_gene_df$comparison)
upset(upset_input_df,comparisons,
      set_sizes=(
        upset_set_size()
        + geom_text(aes(label=..count..), hjust=1.1, stat='count')))

## save uniue comparison
csv_path_int = "/projects/Pregnancy_placenta/Code_PAL/Preg_placenta_R_codes/Figures/pathway_csv_ht_results/Fig4_PTB_vs_Cont_PTB_cell_movement_HBC.csv"
write.csv(mats$pathway_gene_df ,file = csv_path_int)

make_gene_overlap_summary_csv(
  csv_file = csv_path_int,
  output_dir ="/projects/Pregnancy_placenta/Code_PAL/Preg_placenta_R_codes/Figures/upset_plot_results/figure4",
  output_prefix = "fig4_PTB_ContPTB_cell_movement_HBC",
   save_long_csv = FALSE,
)

# HBC - Microtuble Dynamics
mats = get_pathway_gene_heatmap_mats(
  path_of_interest = "Microtubule dynamics",
  list_files = list_files,
  df_loc = all_lin_loc,
  file_regex = "PTB",
  loc_of_IPA_gene_files = "/projects/Pregnancy_placenta/Results/20_Figures/Pathways/IPA_genes/",
  cell_types_keep = c("HBC"),
  comparison_regex = "Cont_PTB"
)

ht = plot_pathway_gene_heatmap(
  log2fc_mat = mats$log2fc_mat,
  ast_mat = mats$ast_mat,
  path_of_interest = "Microtubule dynamics",
  column_title = "Microtubule dynamics - HBC in PTB vs ContPTB",
  annotation_values = ifelse(grepl("TermPE", colnames(mats$log2fc_mat)), "TermPE", "PrePE"),
  annotation_colors = c("TermPE" = "#9ecae1", "PrePE" = "#fdae6b"),
  column_labels = sub("-.*", "", colnames(mats$log2fc_mat))
)


## draw Upset Plotf
upset_input_df <- make_complex_upset_input_genes(
  mats$pathway_gene_df,
  as_logical = FALSE
)
comparisons = unique(mats$pathway_gene_df$comparison)
upset(upset_input_df,comparisons,
      set_sizes=(
        upset_set_size()
        + geom_text(aes(label=..count..), hjust=1.1, stat='count')))

## save uniue comparison
csv_path_int = "/projects/Pregnancy_placenta/Code_PAL/Preg_placenta_R_codes/Figures/pathway_csv_ht_results/Fig4_PTB_vs_Cont_PTB_Microtubule_dynamics_HBC.csv"
write.csv(mats$pathway_gene_df ,file = csv_path_int)

make_gene_overlap_summary_csv(
  csv_file = csv_path_int,
  output_dir ="/projects/Pregnancy_placenta/Code_PAL/Preg_placenta_R_codes/Figures/upset_plot_results/figure4",
  output_prefix = "fig4_PTB_ContPTB_Microtubule_dynamics_HBC",
   save_long_csv = FALSE,
)






