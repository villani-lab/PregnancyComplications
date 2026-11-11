
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
source("/projects/Pregnancy_placenta/Code_PAL/Preg_placenta_R_codes/Figures/funcs_gene_pathway.R")


# load directory of files ----
setwd("/projects/Pregnancy_placenta/Results/13_DE_PAL")

# load file going over where things are located ----
parent_lin = read.csv("/projects/Pregnancy_placenta/Results/13_DE_PAL/output_csv_IPA/redo_Parentlin/Parent_Test_3.csv")
sub_lin = read.csv("/projects/Pregnancy_placenta/Results/13_DE_PAL/output_csv_IPA/redo_sublin/sublin_Test_3.csv")

all_lin_loc = rbind(data.frame(parent_lin),data.frame(sub_lin))

# Load directory of variables----
IPA_pathway_genes_directory = "/projects/Pregnancy_placenta/Results/13_DE_PAL/output_csv_IPA/output_IPA_genes"



list_files = list.files("/projects/Pregnancy_placenta/Results/20_Figures/Pathways/IPA_genes/")

# Fig 5D - Engulfment of cells HBC

## draw Heatmap
mats = get_pathway_gene_heatmap_mats(
  path_of_interest = "Engulfment of cells",
  list_files = list_files,
  df_loc = all_lin_loc,
  file_regex = "PrePE|TermPE",
  loc_of_IPA_gene_files = "/projects/Pregnancy_placenta/Results/20_Figures/Pathways/IPA_genes/",
  cell_types_keep = c("HBC_2","HBC_3"),
  comparison_regex = "Cont_unlabored"
)
## Make heatmap
ht = plot_pathway_gene_heatmap(
  log2fc_mat = mats$log2fc_mat,
  ast_mat = mats$ast_mat,
  path_of_interest = "Engulfment of cells",
  column_title = "Engulfment of cells - HBC_2 and HBC_3 in TermPE vs Cont Unlabored",
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
csv_path_int = "/projects/Pregnancy_placenta/Code_PAL/Preg_placenta_R_codes/Figures/pathway_csv_ht_results/Fig5_TermPE_Cont_unlabored_engulfment_of_cells_HBC.csv"
write.csv(mats$pathway_gene_df ,file = csv_path_int)

make_gene_overlap_summary_csv(
  csv_file = csv_path_int,
  output_dir ="/projects/Pregnancy_placenta/Code_PAL/Preg_placenta_R_codes/Figures/upset_plot_results/figure5",
  output_prefix = "fig5_engulfment_of_cells_HBC",
   save_long_csv = FALSE,
)

# Fig 5D - Cell Movement of Phagocytes

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



# Fig 5F - Migration of Vascular Endothelial

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


## draw Heatmap
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

#########----IGNORE BELOW--------_######
# working directory (IGNORE BELOW)
df <- data.frame(
  full_name = list_files[
    grepl("PrePE|TermPE", list_files) & grepl("-diseases.xls", list_files)
  ]
)

df <- df[
  grepl("SCT|VCTccc|VCTfusion_2|VCT_1", df$full_name),
  ,
  drop = FALSE
]

cell_type <- c()
condition = c()
control = c()
for (i in 1:nrow(df)){
  temp_cell = str_split_fixed(df$full_name[i],"-",3)[1]


 temp_comparison = str_split_fixed(df$full_name[i],"-",3)[2]
  temp_condition = str_split_fixed(temp_comparison,"_vs_",2)[1]
  temp_control = str_split_fixed(temp_comparison,"_vs_",2)[2]

  cell_type[i] = temp_cell
  condition[i] = temp_condition
  control[i] = temp_control
}
df$cell_type = cell_type
df$condition = condition
df$control = control

## Focus on VCTccc vs other trophoblasts

sub_PTB = df
## focus on angiogenesis
Path_of_interest = "Inflammatory response"
loc_of_IPA_gene_files = "/projects/Pregnancy_placenta/Results/20_Figures/Pathways/IPA_genes/"

sub_PTB$Vasculogenesis_in = "N"
sub_PTB$Vasculogenesis_score = NA
sub_PTB_angiogenesis = list()
for (iter in 1:nrow(sub_PTB)){
    temp_all_pathways_cond_int = read_pathway_genes(cell_type = sub_PTB$cell_type[iter],
                            condition = sub_PTB$condition[iter],
                            control = sub_PTB$control[iter],
                            loc_of_IPA_gene_files = loc_of_IPA_gene_files)

  if (Path_of_interest %in% temp_all_pathways_cond_int$`Diseases or Functions Annotation`){
    sub_PTB$Vasculogenesis_in[iter] = "Y"
    sub_PTB$Vasculogenesis_score[iter] =temp_all_pathways_cond_int[temp_all_pathways_cond_int$`Diseases or Functions Annotation`==Path_of_interest,]$`Activation z-score`
    log2_df = get_log2fc_of_pathway_genes(path_int  = Path_of_interest,
                            df_pathway_genes = temp_all_pathways_cond_int,
                            cell_type = sub_PTB$cell_type[iter],
                            condition = sub_PTB$condition[iter],
                            control = sub_PTB$control[iter],
                            df_loc = all_lin_loc)

    name = gsub("-diseases.xls","",sub_PTB$full_name[iter])
    log2_df$comparison = name

    sub_PTB_angiogenesis[[name]] = log2_df


  }
}

sub_PTB_angiogenesis_df = bind_rows(sub_PTB_angiogenesis)
#
# test without PTB
sub_PTB_angiogenesis_df = sub_PTB_angiogenesis_df[grepl("_Cont_unlabored",sub_PTB_angiogenesis_df$comparison),]

# subset out the log2FC<0
#sub_PTB_angiogenesis_df = sub_PTB_angiogenesis_df[sub_PTB_angiogenesis_df$log2fc>0,]
log2fc_df = as.data.frame(sub_PTB_angiogenesis_df%>%
                             pivot_wider(id_cols = 'gene_name',
                                         names_from = 'comparison',
                                         values_from = "log2fc",values_fill = NA))


log2fc_mat = as.matrix(log2fc_df[,-1])
ast_df = as.data.frame(sub_PTB_angiogenesis_df%>%
                             pivot_wider(id_cols = 'gene_name',
                                         names_from = 'comparison',
                                         values_from = "ast",values_fill = NA))

ast_mat =as.matrix(ast_df[,-1])

# generate asterisks
cell_fun = function(j, i, x, y, w, h, fill) {
                 if (!is.na(ast_mat[i, j]) && ast_mat[i, j] != "n.s.") {
                   grid::grid.text(ast_mat[i, j], x, y - grid::unit(0.6, "mm"), just = c("center", "center"))
                 }}
rownames(log2fc_mat) = log2fc_df$gene_name

###-----HEATMAP ASTHETICS----###
cond <- ifelse(grepl("TermPE", colnames(log2fc_mat)),
               "TermPE",
               "PrePE")

top_anno <- HeatmapAnnotation(
    Cond = cond,
    col = list(
        Cond = c(
            "TermPE" = "#9ecae1",
            "PrePE" = "#fdae6b"
        )
    ),
    border = TRUE
    #gp = grid::gpar(col = "black", lwd = 1)   # border around annotation boxes
)

new_colnames <- sub("-.*", "", colnames(log2fc_mat))
#new_colnames[new_colnames == "HBC_2_CXCL2_CXCL3"] <- "HBC_2"
library(grid)

cell_w <- 4   # mm per column
cell_h <- 4   # mm per row



col_fun <- colorRamp2(
  c(
    min(log2fc_mat, na.rm = TRUE),
    0,
    max(log2fc_mat, na.rm = TRUE)
  ),
  c(
    "#2B3CFF",  # blue for negative
    "#FFFEFE",  # white near 0
    "#BB0103"   # red for positive
  )
)


all_angiogenesis = Heatmap(
    log2fc_mat,
    column_labels = new_colnames,
    name = "log2FC",
    col = col_fun,
  na_col = "grey90",
    show_row_names = TRUE,
    top_annotation = top_anno,
    row_names_side = "left",
    column_title = paste0(Path_of_interest,"- Endo Sub PrePE vs Cont Unlbr"),
    cluster_columns = FALSE,
    cluster_rows = FALSE,
   # col = colorRamp2(c(min(log2fc_mat, na.rm = TRUE),max(log2fc_mat, na.rm = TRUE)), c("#FFFEFE", "#BB0103")),
    rect_gp = grid::gpar(col = "black", lwd = 0.5),

    width  = ncol(log2fc_mat) * unit(cell_w, "mm"),
    height = nrow(log2fc_mat) * unit(cell_h, "mm"),

    cell_fun = function(j, i, x, y, w, h, fill) {
        if (!is.na(ast_mat[i, j]) && ast_mat[i, j] != "n.s.") {
            grid::grid.text(
                ast_mat[i, j],
                x,
                y - unit(0.6, "mm"),
                just = c("center", "center")
            )
        }
    }
)

# FOR SUB FILTERED MATRIX !

filt_log2mat <- function(mat,
                         cell_subtypes =  sub("-.*", "", colnames(mat)),
                         top_or_lowest_FC =rep("top",times = ncol(mat)),
                         top_n_genes = 5,
                         add_annotations = "Y"){
  top_gene <- list()

  for(i in 1:ncol(mat)){
    top_or_bottom = top_or_lowest_FC[i]
    if(top_or_bottom =="top"){
      gene_list = names(mat[,i][order(mat[,i],na.last = TRUE,decreasing=TRUE)][1:top_n_genes])
    }else{
      gene_list = names(mat[,i][order(mat[,i],na.last = TRUE)][1:top_n_genes])
    }
    top_gene[[cell_subtypes[i]]]=gene_list
  }

    gene_celltype_df <- do.call(
      rbind,
      lapply(names(top_gene), function(ct) {
        data.frame(
          gene =top_gene[[ct]],
          Cell_Type = ct,
          stringsAsFactors = FALSE
        )
      })
    )

  gene_celltype_df <- aggregate(
    Cell_Type ~ gene,
    data = gene_celltype_df,
    FUN = function(x) paste(unique(x), collapse = " + ")
  )
  gene_order <- unique(unlist(top_gene))

  gene_celltype_df <- gene_celltype_df[
    match(gene_order, gene_celltype_df$gene),
  ]

  rownames(gene_celltype_df) <- NULL

 # gene_celltype_df
  return(gene_celltype_df)
}


filt_log2mat_5 = filt_log2mat(log2fc_mat,cell_subtypes =  sub("_vs_.*", "", colnames(log2fc_mat)),
                              top_or_lowest_FC = c("top",'top','top',"top",'bottom','top'))
filt_log2mat_mat = log2fc_mat[rownames(log2fc_mat) %in% filt_log2mat_5$gene,]
 v <- log2fc_mat[
  filt_log2mat_5$gene,
  ,
  drop = FALSE
]
# Row annotation dataframe
celltype_cols <- c(
  "SCT-PrePE" = "#E41A1C",
  "SCT-TermPE" = "#E41A1C",
  "SCT-PrePE + SCT-TermPE" = "#E41A1C",

  "VCT_1-PrePE" = "#377EB8",
  "VCT_1-TermPE" = "#377EB8",
  "VCT_1-PrePE + VCT_1-TermPE" = "#377EB8",

  "VCTfusion_2-PrePE" = "#4DAF4A",

  "SCT-TermPE + VCTfusion_2-PrePE" = "#984EA3",
  "VCT_1-PrePE + VCTfusion_2-PrePE" = "#984EA3",

  "VCTccc-PrePE" = "#FF7F00"
)
# Rebuild row annotation directly in the order of the heatmap matrix
row_anno_df <- data.frame(
  gene = rownames(filt_log2mat_mat),
  Cell_Type = filt_log2mat_5$Cell_Type[
    match(rownames(filt_log2mat_mat), filt_log2mat_5$gene)
  ],
  stringsAsFactors = FALSE
)

rownames(row_anno_df) <- row_anno_df$gene

# Check mapping before plotting
row_anno_df

# Row annotation object
# Row annotation object
row_ha <- rowAnnotation(
  Cell_Type = row_anno_df$Cell_Type,
  col = list(
    Cell_Type = celltype_cols
  ),
  annotation_name_side = "top",
  show_annotation_name = TRUE
)

col_fun <- colorRamp2(
  c(
    min(filt_log2mat_mat, na.rm = TRUE),
    0,
    max(filt_log2mat_mat, na.rm = TRUE)
  ),
  c(
    "#2B3CFF",  # blue for negative
    "#FFFEFE",  # white near 0
    "#BB0103"   # red for positive
  )
)
rownames(ast_mat) = rownames(log2fc_mat)
filt_ast_mat <- ast_mat[
  filt_log2mat_5$gene,
  ,
  drop = FALSE
]
sub_all_angiogenesis = Heatmap(
   filt_log2mat_mat,
    column_labels = new_colnames,
    name = "log2FC",
    col = col_fun,
  na_col = "grey90",
    show_row_names = TRUE,
    right_annotation = row_ha,
    top_annotation = top_anno,
    row_names_side = "left",
    column_title = paste0(Path_of_interest,"-  PrePE vs Cont Unlbr"),
    cluster_columns = FALSE,
    cluster_rows = FALSE,
   # col = colorRamp2(c(min(log2fc_mat, na.rm = TRUE),max(log2fc_mat, na.rm = TRUE)), c("#FFFEFE", "#BB0103")),
    rect_gp = grid::gpar(col = "black", lwd = 0.5),

    width  = ncol(filt_log2mat_mat) * unit(cell_w, "mm"),
    height = nrow(filt_log2mat_mat) * unit(cell_h, "mm"),

    cell_fun = function(j, i, x, y, w, h, fill) {
        if (!is.na(filt_ast_mat [i, j]) && filt_ast_mat [i, j] != "n.s.") {
            grid::grid.text(
                filt_ast_mat [i, j],
                x,
                y - unit(0.6, "mm"),
                just = c("center", "center")
            )
        }
    }
)

# save heatmap

pdf(file="andy_pathway_gene_heatmaps/vasculogenesis_HBC_PTB_all_controls.pdf",height=22,width=4.8)
draw(all_angiogenesis)
dev.off()

pdf(file="andy_pathway_gene_heatmaps/vasculogenesis_HBC_PTB_ctrl-lbr_ctrl_unlbr.pdf",height=16,width=4.8)
draw(all_angiogenesis)
dev.off()


#------------------IGNORE BELOW
#loc_of_IPA_gene_files = "/projects/Pregnancy_placenta/Results/20_Figures/Pathways/IPA_genes/"
#Example, plotting the genes implicated in Angionesis across three condtions:
# HBC_1-PTB_vs_Cont_PTB
# HBC_10_cycling-PTB_vs_Cont_PTB
# VCT_p2-PTB_vs_Cont_PTB

## Condition 1----
path_int ="Angiogenesis"
cell_type = "HBC_1"
condition = "PTB"
control = "Cont_PTB"
loc_of_IPA_gene_files = "/projects/Pregnancy_placenta/Results/20_Figures/Pathways/IPA_genes/"
temp_all_pathways_cond_int = read_pathway_genes(cell_type = cell_type,
                            condition = condition,
                            control = control,
                            loc_of_IPA_gene_files = loc_of_IPA_gene_files)



log2fc_genes_hbc1_ptb_contptb = get_log2fc_of_pathway_genes(pathway = path_int,
                            df_pathway_genes = temp_all_pathways_cond_int,
                            cell_type = cell_type,
                            condition = condition,
                            control = control,
                            df_loc = all_lin_loc)
## Condition 2----
log2fc_genes_hbc1_ptb_contptb$condition = "HBC1"

path_int ="Angiogenesis"
cell_type = "HBC_10_cycling"
condition = "PTB"
control = "Cont_PTB"

temp_all_pathways_cond_int = read_pathway_genes(cell_type = cell_type,
                                                condition = condition, control = control,
                                                loc_of_IPA_gene_files = loc_of_IPA_gene_files)

log2fc_genes_hbc10_ptb_contlptb = get_log2fc_of_pathway_genes(pathway = path_int,
                            df_pathway_genes = temp_all_pathways_cond_int,
                            cell_type = cell_type,
                            condition = condition,
                            control = control,
                            df_loc = all_lin_loc)

log2fc_genes_hbc10_ptb_contlptb$condition = "HBC10"

## Condition 3----
path_int ="Angiogenesis"
cell_type = "VCTp_2"
condition = "PTB"
control = "Cont_PTB"

temp_all_pathways_cond_int = read_pathway_genes(cell_type = cell_type,
                                                condition = condition, control = control,
                                                loc_of_IPA_gene_files = loc_of_IPA_gene_files)

log2fc_genes_vctp2_ptb_contlptb = get_log2fc_of_pathway_genes(pathway = path_int,
                            df_pathway_genes = temp_all_pathways_cond_int,
                            cell_type = cell_type,
                            condition = condition,
                            control = control,
                            df_loc = all_lin_loc)

log2fc_genes_vctp2_ptb_contlptb$condition = "VCTp2"


## merge pathways and prep for Heatmap ----
merged_pathways = rbind(log2fc_genes_hbc1_ptb_contptb,log2fc_genes_hbc10_ptb_contlptb,log2fc_genes_vctp2_ptb_contlptb)
# make heatmap
log2fc_df = as.data.frame(merged_pathways%>%
                             pivot_wider(id_cols = 'gene_name',
                                         names_from = 'condition',
                                         values_from = "log2fc",values_fill = NA))

ast_df = as.data.frame(merged_pathways %>%
                             pivot_wider(id_cols = 'gene_name',
                                         names_from = 'condition',
                                         values_from = "ast",values_fill = NA))



log2fc_mat = as.matrix(log2fc_df[,-1])
rownames(log2fc_mat) = log2fc_df$gene_name
Heatmap(log2fc_mat,
        name = "log2FC",
        show_row_names = TRUE,
        row_names_side = "left",
        column_title = "Angiogenesis - PTB vs Cont PTB",
        cluster_columns = FALSE,
        cluster_rows = FALSE
)

