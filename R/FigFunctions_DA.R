suppressPackageStartupMessages({
  library(openxlsx)
  library(tidyr)
  library(stringr)
  library(ComplexHeatmap)
  library(colorRamp2)

})


prep_df_stat_for_plotting = function(input_df)
{
  # get the DF of the log2FC, qvalue and pvalue
  order_comp1 = c("PrePE","TermPE","FGR","DM1","PTB","Cont_labored")
  order_comp2 = c("Cont_labored","Cont_unlabored","Cont_PTB")

  df_temp = as.data.frame(input_df)
  df_temp$Comp1 = stringr::str_split_fixed(df_temp$Comparison,"-",2)[,1]
  df_temp$Comp2 = stringr::str_split_fixed(df_temp$Comparison,"-",2)[,2]

  df_temp$Comp2 = factor(df_temp$Comp2,levels=order_comp2)
  df_temp$Comp1 = factor(df_temp$Comp1, levels = order_comp1)

  df_temp = df_temp[order(df_temp$Comp2,df_temp$Comp1),]

  ## get the log2FC
  log2FC_temp = as.data.frame(df_temp %>%
                                pivot_wider(id_cols = 'Cluster',
                                            names_from = 'Comparison',
                                            values_from = "log.FC."))

  rownames(log2FC_temp) = log2FC_temp$Cluster
  log2FC_temp = log2FC_temp[,-c(1)]

  ## get qvalue
  qvalue_temp = df_temp
  qvalue_temp$qval_ast = stars.pval(qvalue_temp$qval)

  qvalue_temp = as.data.frame(qvalue_temp %>%
                                pivot_wider(id_cols = 'Cluster',
                                            names_from = 'Comparison',
                                            values_from = "qval_ast"))

  rownames(qvalue_temp) = qvalue_temp$Cluster
  qvalue_temp = qvalue_temp[,-c(1)]

  ## get pvalue
  pvalue_temp = df_temp
  pvalue_temp$pval_ast = stars.pval(pvalue_temp$P.val)

  pvalue_temp = as.data.frame(pvalue_temp %>%
                                pivot_wider(id_cols = 'Cluster',
                                            names_from = 'Comparison',
                                            values_from = "pval_ast"))

  rownames(pvalue_temp) = pvalue_temp$Cluster
  pvalue_temp = pvalue_temp[,-c(1)]

  mat_list = list(log2_df = log2FC_temp,
                  qvalue_df = qvalue_temp,
                  pvalue_df = pvalue_temp)

  return(mat_list)
}


stars.pval <- function(x){
  stars <- c("***", "**", "*", "n.s.")
  cutoffs <- c(0, 0.001, 0.01, 0.05, 1)
  i <- findInterval(x, cutoffs, left.open = TRUE, rightmost.closed = TRUE)
  stars[i]
}




DA_heatmaps <- function(lineages_of_interest,
                             log2FC_list,
                             qval_list,
                             pval_list,
                             plot_row_order = FALSE,
                             desired_row_order=NULL,  # New argument for custom order
                             order_comp2 = c("Cont_labored", "Cont_unlabored", "Cont_PTB"),
                             plot_type = c("qval_ast", "pval_ast", "both"),
                             save_pdf = TRUE,
                             output_prefix = "DA_heatmap",
                             output_dir = ".",
                             pdf_width = 9,
                             pdf_height = 6,
                             row_split_by=data.frame())
{
  plot_type <- match.arg(plot_type)

  # Merge log2FC
  merged_df_log2FC <- do.call(rbind, lapply(lineages_of_interest, function(name)
  {
    df <- log2FC_list[[name]]
    df$Lineage <- name
    return(df)
  }))

  # Merge qvals and pvals
  merged_df_qval <- do.call(rbind, lapply(lineages_of_interest, function(name)
  {
    qval_list[[name]]
  }))

  merged_df_pval <- do.call(rbind, lapply(lineages_of_interest, function(name)
  {
    pval_list[[name]]
  }))

  # Prepare heatmap matrix
  mat_log2FC_ht <- as.matrix(merged_df_log2FC[, -ncol(merged_df_log2FC)])
  mat_column_labels <- stringr::str_split_fixed(colnames(mat_log2FC_ht), "-", 2)[, 1]
  column_split_raw <- stringr::str_split_fixed(colnames(mat_log2FC_ht), "-", 2)[, 2]

  # Visually separate 'Cont_unlabored' column (adjust logic if needed)
  column_split_modified <- ifelse(colnames(mat_log2FC_ht) == "Cont_labored-Cont_unlabored",
                                  "Cont_unlabored_sep",
                                  column_split_raw)

  column_split_cols <- factor(column_split_modified, levels = c("Cont_labored", "Cont_unlabored", "Cont_PTB", "Cont_unlabored_sep"))

  cond2_blocks <- c("Cont_labored" = "#0073C2FF",
                    "Cont_unlabored" =  "#EFC000FF",
                    "Cont_PTB" = "#CD534CFF")

  # Remap for annotation coloring/legend
  annotation_split <- ifelse(column_split_cols == "Cont_unlabored_sep",
                             "Cont_unlabored",
                             as.character(column_split_cols))

  annotation_split <- factor(annotation_split, levels = c("Cont_labored", "Cont_unlabored", "Cont_PTB"))

  ## Need to reorder everythibng based on column split cols
  annotation_split = annotation_split[order(column_split_cols)]
  mat_log2FC_ht = mat_log2FC_ht[,order(column_split_cols)]
  merged_df_qval = merged_df_qval[,order(column_split_cols)]
  merged_df_pval = merged_df_pval[,order(column_split_cols)]

  Cond_annotation <- ComplexHeatmap::HeatmapAnnotation(Ctrl = annotation_split,
                                                       col = list(Ctrl = cond2_blocks),
                                                       annotation_legend_param = list(
                                                         Ctrl = list(nrow = 1, title_position = "topcenter", direction = "horizontal")),
                                                       show_annotation_name = FALSE,
                                                       border = TRUE)

  log2FC_col <- circlize::colorRamp2(c(min(mat_log2FC_ht, na.rm = TRUE), 0, max(mat_log2FC_ht, na.rm = TRUE)), c("#2B8F0D", "#FFFEFE", "#A62991")  )

  row_split_vec <- row_split_by[rownames(mat_log2FC_ht), "split", drop = TRUE]
  row_split_vec <- factor(row_split_vec, levels = levels(row_split_by[,"split"]), ordered = TRUE)
  stopifnot(all(rownames(mat_log2FC_ht) %in% rownames(row_split_by)))  # fail loudly if any rows are missing

  # Custom row order: filter to valid rows
  if(plot_row_order == TRUE)
  {
    valid_order <- desired_row_order[desired_row_order %in% rownames(mat_log2FC_ht)]
    # Reorder all matrices consistently
    mat_log2FC_ht <- mat_log2FC_ht[valid_order, , drop = FALSE]
    merged_df_qval <- merged_df_qval[valid_order, , drop = FALSE]
    merged_df_pval <- merged_df_pval[valid_order, , drop = FALSE]
    merged_df_log2FC <- merged_df_log2FC[match(valid_order, rownames(mat_log2FC_ht)), , drop = FALSE]
    row_split_vec <- row_split_by[rownames(mat_log2FC_ht), "split", drop = TRUE]
    row_split_vec <- factor(row_split_vec, levels = levels(row_split_by[,"split"]), ordered = TRUE)

  }  else{valid_order = rownames(mat_log2FC_ht)}

  cell_mm <- 4
  make_heatmap <- function(ast_matrix) {
    ComplexHeatmap::Heatmap(mat_log2FC_ht,
                            name = "log2FC - DA\nCond / Ctrl",
                            cluster_rows = FALSE,
                            row_order = valid_order,
                            cluster_columns = FALSE,
                            border_gp = grid::gpar(col = "black", lwd = 1),
                            rect_gp = grid::gpar(col = "black", lwd = 0.5),
                            top_annotation = Cond_annotation,
                            row_title_side = "right",
                            row_title_rot = 0,
                            column_title = NULL,
                            row_names_side = "left",
                            col = log2FC_col,
                            row_split = row_split_vec,
                            column_labels = mat_column_labels[order(column_split_cols)],
                            column_split = column_split_cols[order(column_split_cols)],
                            column_title_gp = grid::gpar(fontface = "bold"),
                            cell_fun = function(j, i, x, y, w, h, fill) {
                              if (!is.na(ast_matrix[i, j]) && ast_matrix[i, j] != "n.s.") {
                                grid::grid.text(ast_matrix[i, j], x, y = y - 0.1 * h,          # <- tweak: to make the asteriks int he middle of the cell
                                                just = c("centre", "centre"),)
                              }
                            },
                           width  = unit(ncol(mat_log2FC_ht) * cell_mm, "mm"),
                           height = unit(nrow(mat_log2FC_ht) * cell_mm, "mm"))
  }
  if (plot_type %in% c("qval_ast", "both"))
  {
    ht_q <- make_heatmap(as.matrix(merged_df_qval))
    if (save_pdf)
    {
      pdf(file = file.path(output_dir, paste0(output_prefix, "_qval_heatmap.pdf")),
          width = pdf_width, height = pdf_height)
      ComplexHeatmap::draw(ht_q, annotation_legend_side = "top",
                           column_title = "Conditions",
                           column_title_side = "bottom",
                           column_title_gp = grid::gpar(fontface = "bold"))
      dev.off()
      message("Saved qval heatmap to PDF.")
    } else
    {
      ComplexHeatmap::draw(ht_q, annotation_legend_side = "top",
                           column_title = "Conditions",
                           column_title_side = "bottom",
                           column_title_gp = grid::gpar(fontface = "bold"))
    }
  }

  if (plot_type %in% c("pval_ast", "both"))
  {
    ht_p <- make_heatmap(as.matrix(merged_df_pval))
    if (save_pdf)
    {
      pdf(file = file.path(output_dir, paste0(output_prefix, "_pval_heatmap.pdf")),
          width = pdf_width, height = pdf_height)
      ComplexHeatmap::draw(ht_p, annotation_legend_side = "top",
                           column_title = "Conditions",
                           column_title_side = "bottom",
                           column_title_gp = grid::gpar(fontface = "bold"))
      dev.off()
      message("Saved pval heatmap to PDF.")
    } else {
      ComplexHeatmap::draw(ht_p, annotation_legend_side = "top",
                           column_title = "Conditions",
                           column_title_side = "bottom",
                           column_title_gp = grid::gpar(fontface = "bold"))
    }
  }
}
