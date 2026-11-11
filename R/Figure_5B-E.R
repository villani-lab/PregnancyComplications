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


# ======================================================================
# Figure 5C-E
# ======================================================================

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




# load functions ----
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



list_files = list.files("/projects/Pregnancy_placenta/Results/20_Figures/Pathways/IPA_genes/")

# Fig 5C - Cell Movement of Phagocytes
#===================================================================================================
## draw Heatmap
mats = get_pathway_gene_heatmap_mats(
  path_of_interest = "Cell movement of phagocytes",
  list_files = list_files,
  df_loc = all_lin_loc,
  file_regex = "PrePE|TermPE",
  loc_of_IPA_gene_files = "/projects/Pregnancy_placenta/Results/20_Figures/Pathways/IPA_genes/",
  cell_types_keep = c("HBC_2","HBC_3","CD14"),
  comparison_regex = "Cont_unlabored"
)
## Make heatmap
ht = plot_pathway_gene_heatmap(
  log2fc_mat = mats$log2fc_mat,
  ast_mat = mats$ast_mat,
  path_of_interest = "Cell movement of phagocytes",
  column_title = "Cell movement of phagocytes - HBC_2 and HBC_3 in TermPE vs Cont Unlabored",
  annotation_values = ifelse(grepl("TermPE", colnames(mats$log2fc_mat)), "TermPE", "PrePE"),
  annotation_colors = c("TermPE" = "#9ecae1", "PrePE" = "#fdae6b"),
  column_labels = sub("-.*", "", colnames(mats$log2fc_mat))
)

draw(ht)

## draw Upset Plotf
comp_keep<-c("HBC_3-TermPE_vs_Cont_unlabored",
             'HBC_2_CXCL2_CXCL3-TermPE_vs_Cont_unlabored',
             'HBC_2_CXCL2_CXCL3-PrePE_vs_Cont_unlabored',
             'CD14_monocyte-PrePE_vs_Cont_unlabored')
upset_input_df <- make_complex_upset_input_genes(
  subset(mats$pathway_gene_df,comparison %in% comp_keep),
  as_logical = FALSE
)
comparisons = unique(mats$pathway_gene_df$comparison)
upset(upset_input_df,comp_keep)

## save unique comparison
csv_path_int = "/projects/Pregnancy_placenta/Code_PAL/Preg_placenta_R_codes/Figures/pathway_csv_ht_results/Fig5_Cell_Movement_of_Phagocytes_HBC_CD14_V2.csv"
write.csv(subset(mats$pathway_gene_df,comparison %in% comp_keep) ,file = csv_path_int)

make_gene_overlap_summary_csv(
  csv_file = csv_path_int,
  output_dir ="/projects/Pregnancy_placenta/Code_PAL/Preg_placenta_R_codes/Figures/upset_plot_results/figure5",
  output_prefix = "fig5_movement_of_phagocytes_V2",
   save_long_csv = FALSE,
)


# Fig 5D - Interferon ALpha beta pathway
#===================================================================================================

## draw Heatmap
mats = get_pathway_gene_heatmap_mats(
  path_of_interest = "Interferon alpha/beta signaling",
  list_files = list_files,
  df_loc = all_lin_loc,
  file_regex = "PrePE|TermPE",
  loc_of_IPA_gene_files = "/projects/Pregnancy_placenta/Results/20_Figures/Pathways/IPA_genes/",
  cell_types_keep = c("HBC_2","HBC_3"),
  comparison_regex = "Cont_unlabored",
  diseases_or_pathway = 'PW'
)
## Make heatmap
ht = plot_pathway_gene_heatmap(
  log2fc_mat = mats$log2fc_mat,
  ast_mat = mats$ast_mat,
  path_of_interest = "Interferon alpha/beta signaling",
  column_title = "Interferon alpha/beta signaling - HBC_2 and HBC_3 in TermPE vs Cont Unlabored",
  annotation_values = ifelse(grepl("TermPE", colnames(mats$log2fc_mat)), "TermPE", "PrePE"),
  annotation_colors = c("TermPE" = "#9ecae1", "PrePE" = "#fdae6b"),
  column_labels = sub("-.*", "", colnames(mats$log2fc_mat))
)

draw(ht)

## draw Upset Plotf

upset_input_df <- make_complex_upset_input_genes(
 mats$pathway_gene_df,
  as_logical = FALSE
)
comparisons = unique(mats$pathway_gene_df$comparison)
upset(upset_input_df,comparisons)

## save unique comparison
csv_path_int = "/projects/Pregnancy_placenta/Code_PAL/Preg_placenta_R_codes/Figures/pathway_csv_ht_results/Fig5D_Interferron_signaling_HBC.csv"
write.csv(mats$pathway_gene_df ,file = csv_path_int)

make_gene_overlap_summary_csv(
  csv_file = csv_path_int,
  output_dir ="/projects/Pregnancy_placenta/Code_PAL/Preg_placenta_R_codes/Figures/upset_plot_results/figure5",
  output_prefix = "fig5_interferron_signaling_V2",
   save_long_csv = FALSE,
)



# Fig 5E - Migration of Vascular Endothelial
#===================================================================================================
## draw Heatmap
mats = get_pathway_gene_heatmap_mats(
  path_of_interest = "Migration of vascular endothelial cells",
  list_files = list_files,
  df_loc = all_lin_loc,
  file_regex = "PrePE",
  loc_of_IPA_gene_files = "/projects/Pregnancy_placenta/Results/20_Figures/Pathways/IPA_genes/",
  cell_types_keep = c("Endo"),
  comparison_regex = "Cont_unlabored"
)
## Make heatmap
ht = plot_pathway_gene_heatmap(
  log2fc_mat = mats$log2fc_mat,
  ast_mat = mats$ast_mat,
  path_of_interest = "Migration of vascular endothelial cells",
  column_title = "Migration of vascular endothelial cells - Endo PrePE vs Cont Unlabored",
  annotation_values = ifelse(grepl("TermPE", colnames(mats$log2fc_mat)), "TermPE", "PrePE"),
  annotation_colors = c("TermPE" = "#9ecae1", "PrePE" = "#fdae6b"),
  column_labels = sub("-.*", "", colnames(mats$log2fc_mat))
)

draw(ht)

## draw Upset Plotf
upset_input_df <- make_complex_upset_input_genes(
  mats$pathway_gene_df,
  as_logical = FALSE
)
comparisons = unique(mats$pathway_gene_df$comparison)
upset(upset_input_df,comparisons)

## save unique comparison
csv_path_int = "/projects/Pregnancy_placenta/Code_PAL/Preg_placenta_R_codes/Figures/pathway_csv_ht_results/Fig5_PrePE_Cont_unlabored_migration_vasc_endo.csv"
write.csv(mats$pathway_gene_df ,file = csv_path_int)

make_gene_overlap_summary_csv(
  csv_file = csv_path_int,
  output_dir ="/projects/Pregnancy_placenta/Code_PAL/Preg_placenta_R_codes/Figures/upset_plot_results/figure5",
  output_prefix = "fig5_migration_endo_prepe",
   save_long_csv = FALSE,
)


# Fig 5D - Inflammatory response in trophoblasts
#===================================================================================================
mats = get_pathway_gene_heatmap_mats(
  path_of_interest = "Inflammatory response",
  list_files = list_files,
  df_loc = all_lin_loc,
  file_regex = "PrePE",
  loc_of_IPA_gene_files = "/projects/Pregnancy_placenta/Results/20_Figures/Pathways/IPA_genes/",
  cell_types_keep = c("SCT","VCTfusion_2","VCTccc"),
  comparison_regex = "Cont_unlabored"
)
## Make heatmap
ht = plot_pathway_gene_heatmap(
  log2fc_mat = mats$log2fc_mat,
  ast_mat = mats$ast_mat,
  path_of_interest = "Inflammatory response",
  column_title = "Inflammatory response - Troph PrePE vs Cont Unlabored",
  annotation_values = ifelse(grepl("TermPE", colnames(mats$log2fc_mat)), "TermPE", "PrePE"),
  annotation_colors = c("TermPE" = "#9ecae1", "PrePE" = "#fdae6b"),
  column_labels = sub("-.*", "", colnames(mats$log2fc_mat))
)

draw(ht)

## draw Upset Plotf
upset_input_df <- make_complex_upset_input_genes(
  mats$pathway_gene_df,
  as_logical = FALSE
)
comparisons = unique(mats$pathway_gene_df$comparison)
upset(upset_input_df,comparisons)

## save unique comparison
csv_path_int = "/projects/Pregnancy_placenta/Code_PAL/Preg_placenta_R_codes/Figures/pathway_csv_ht_results/Fig5_PrePE_Cont_unlabored_inflammatory_response_Troph.csv"
write.csv(mats$pathway_gene_df ,file = csv_path_int)

make_gene_overlap_summary_csv(
  csv_file = csv_path_int,
  output_dir ="/projects/Pregnancy_placenta/Code_PAL/Preg_placenta_R_codes/Figures/upset_plot_results/figure5",
  output_prefix = "fig5_inflammatory_response_troph",
   save_long_csv = FALSE,
)

### MISC FUNCTIONS
make_complex_upset_input_genes <- function(pathway_gene_df,
                                           comparison_col = "comparison",
                                           gene_col = "gene_name",
                                           logfc_col = "log2fc",
                                           num_pos_col = "Num_Pos",
                                           as_logical = FALSE) {

  library(dplyr)
  library(tidyr)

  df <- pathway_gene_df %>%
    mutate(
      "{logfc_col}" := as.numeric(.data[[logfc_col]])
    )

  # Count number of comparisons where each gene has positive log2fc
  num_pos_df <- df %>%
    group_by(.data[[gene_col]]) %>%
    summarise(
      "{num_pos_col}" := sum(.data[[logfc_col]] > 0, na.rm = TRUE),
      .groups = "drop"
    )

  # Make gene x comparison presence/absence matrix
  gene_binary_df <- df %>%
    distinct(.data[[gene_col]], .data[[comparison_col]]) %>%
    mutate(present = ifelse(as_logical, TRUE, 1)) %>%
    pivot_wider(
      names_from = .data[[comparison_col]],
      values_from = present,
      values_fill = ifelse(as_logical, FALSE, 0)
    )

  # Add Num_Pos and place it after gene_name
  out_df <- gene_binary_df %>%
    left_join(num_pos_df, by = gene_col) %>%
    relocate(all_of(num_pos_col), .after = all_of(gene_col))

  return(out_df)
}
make_complex_upset_input <- function(pathway_gene_df,
                                     comparison_col = "comparison",
                                     gene_col = "gene_name",
                                     logfc_col = "log2fc",
                                     num_pos_col = "Num_Pos",
                                     as_logical = FALSE) {

  library(dplyr)
  library(tidyr)

  df <- pathway_gene_df %>%
    mutate(
      "{logfc_col}" := as.numeric(.data[[logfc_col]])
    )

  # Count number of positive genes per comparison
  num_pos_df <- df %>%
    group_by(.data[[comparison_col]]) %>%
    summarise(
      "{num_pos_col}" := sum(.data[[logfc_col]] > 0, na.rm = TRUE),
      .groups = "drop"
    )

  # Make comparison x gene presence/absence matrix
  gene_binary_df <- df %>%
    distinct(.data[[comparison_col]], .data[[gene_col]]) %>%
    mutate(present = ifelse(as_logical, TRUE, 1)) %>%
    pivot_wider(
      names_from = .data[[gene_col]],
      values_from = present,
      values_fill = ifelse(as_logical, FALSE, 0)
    )

  # Add Num_Pos and place it after comparison
  out_df <- gene_binary_df %>%
    left_join(num_pos_df, by = comparison_col) %>%
    relocate(all_of(num_pos_col), .after = all_of(comparison_col))

  return(out_df)
}

# preparing heatmap for UpSet Plot
mats_long_df = as.data.frame(mats$log2fc_mat)
mats_long_df$gene_name = rownames(mats$log2fc_mat)
rownames(mats_long_df) = NULL
mat_long_df = mats_long_df %>%
  pivot_longer(!gene_name,names_to = 'comparison' , values_to = 'log2fc')
mat_long_df = na.omit(mat_long_df)
write.csv(mat_long_df ,file = "/projects/Pregnancy_placenta/Code_PAL/Preg_placenta_R_codes/Figures/pathway_csv_ht_results/Fig4_PTB_vs_Cont_PTB_Troph_invasion_of_cells.csv")

## upset plot
csv_path_int = "/projects/Pregnancy_placenta/Code_PAL/Preg_placenta_R_codes/Figures/pathway_csv_ht_results/Fig4_PTB_vs_Cont_PTB_Troph_invasion_of_cells.csv"
make_gene_overlap_summary_csv(
  csv_file = csv_path_int,
  output_dir ="/projects/Pregnancy_placenta/Code_PAL/Preg_placenta_R_codes/Figures/upset_plot_results",
  output_prefix = "fig4_invasion_of_cells_upset_plot_summary_csv"
)

example_result <- make_gene_overlap_upset(
    csv_file = csv_path_int,
    comparisons = "all",
    positive_color = "black",
    negative_color = "red",
    save_pdf = FALSE,
    output_dir = getwd(),
    output_file = "inflammatory_response_gene_overlap_upset.pdf"
  )
