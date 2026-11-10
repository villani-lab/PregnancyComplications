
# Helper functions used by SupFig5_other_control_analysis.R and
# SupFig6_DM1_FGR.R.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
})

normalize_annotation_label <- function(x) {
  x <- as.character(x)
  x <- gsub("\\\\n", "\n", x)
  trimws(x)
}

heatmap_font_family <- "Arial"

open_heatmap_pdf_device <- function(file, width, height, family = heatmap_font_family) {
  if (capabilities("cairo")) {
    grDevices::cairo_pdf(
      filename = file,
      width = width,
      height = height,
      family = family
    )
  } else {
    warning(
      "Cairo graphics are not available; falling back to default PDF device and `sans` font."
    )
    grDevices::pdf(
      file = file,
      width = width,
      height = height,
      family = "sans"
    )
  }
}

prepare_pathway_row_split <- function(pathway_names,
                                      pathway_annotation_df = NULL,
                                      pathway_annotation_file = NULL,
                                      pathway_annotation_sheet = 1,
                                      pathway_col = "Pathway",
                                      annotation_col = "Annotation",
                                      annotation_order = NULL,
                                      unmatched_annotation = "Other pathways") {

  if (is.null(pathway_annotation_df) && is.null(pathway_annotation_file)) {
    return(list(
      row_order = pathway_names,
      row_split = NULL
    ))
  }

  if (is.null(pathway_annotation_df)) {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stop("Package `readxl` is required when `pathway_annotation_file` is provided.")
    }

    pathway_annotation_df <- readxl::read_excel(
      path = pathway_annotation_file,
      sheet = pathway_annotation_sheet
    )
  }

  if (!all(c(pathway_col, annotation_col) %in% colnames(pathway_annotation_df))) {
    stop(
      paste0(
        "`pathway_annotation_df` must contain columns `",
        pathway_col,
        "` and `",
        annotation_col,
        "`."
      )
    )
  }

  pathway_annotation_df <- pathway_annotation_df %>%
    dplyr::transmute(
      Pathway = trimws(as.character(.data[[pathway_col]])),
      Annotation = normalize_annotation_label(.data[[annotation_col]])
    ) %>%
    dplyr::filter(!is.na(Pathway), Pathway != "") %>%
    dplyr::mutate(
      Annotation = dplyr::if_else(
        is.na(Annotation) | Annotation == "",
        unmatched_annotation,
        Annotation
      )
    ) %>%
    dplyr::distinct(Pathway, .keep_all = TRUE)

  if (nrow(pathway_annotation_df) == 0) {
    return(list(
      row_order = pathway_names,
      row_split = NULL
    ))
  }

  pathway_names_clean <- trimws(as.character(pathway_names))
  matched_annotations <- stats::setNames(
    pathway_annotation_df$Annotation,
    pathway_annotation_df$Pathway
  )

  row_annotations <- matched_annotations[pathway_names_clean]
  if (!is.null(unmatched_annotation)) {
    row_annotations[is.na(row_annotations) | row_annotations == "NA"] <- unmatched_annotation
  }

  annotated_pathways <- pathway_annotation_df$Pathway[pathway_annotation_df$Pathway %in% pathway_names_clean]
  remaining_pathways <- pathway_names_clean[!pathway_names_clean %in% annotated_pathways]
  row_order <- c(annotated_pathways, remaining_pathways)

  if (is.null(annotation_order)) {
    annotation_order <- unique(pathway_annotation_df$Annotation)
  }

  if (!is.null(unmatched_annotation) && any(row_annotations == unmatched_annotation, na.rm = TRUE)) {
    annotation_order <- c(annotation_order[annotation_order != unmatched_annotation], unmatched_annotation)
  }

  row_split <- row_annotations[row_order]
  if (all(is.na(row_split))) {
    row_split <- NULL
  } else {
    row_split <- factor(row_split, levels = unique(annotation_order))
  }

  list(
    row_order = row_order,
    row_split = row_split
  )
}

filter_pathways_by_annotation <- function(pathway_names,
                                          pathway_annotation_df = NULL,
                                          pathway_annotation_file = NULL,
                                          pathway_annotation_sheet = 1,
                                          pathway_col = "Pathway",
                                          annotation_col = "Annotation",
                                          filter_annotation_group = "None",
                                          keep_annotation_group = NULL,
                                          filter_pathway = NULL,
                                          keep_pathway = NULL,
                                          omit_unannotated_pathways = FALSE) {

  pathway_names_clean <- trimws(as.character(pathway_names))

  if (!is.null(keep_pathway)) {
    keep_pathway <- trimws(as.character(keep_pathway))
    pathway_names_clean <- pathway_names_clean[pathway_names_clean %in% keep_pathway]
  }

  if (!is.null(filter_pathway)) {
    filter_pathway <- trimws(as.character(filter_pathway))
    pathway_names_clean <- pathway_names_clean[!pathway_names_clean %in% filter_pathway]
  }

  if (identical(filter_annotation_group, "None") && !isTRUE(omit_unannotated_pathways)) {
    return(pathway_names_clean)
  }

  if (is.null(pathway_annotation_df) && is.null(pathway_annotation_file)) {
    stop("A pathway annotation table or file is required when pathway annotation filtering is used.")
  }

  if (is.null(pathway_annotation_df)) {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stop("Package `readxl` is required when `pathway_annotation_file` is provided.")
    }

    pathway_annotation_df <- readxl::read_excel(
      path = pathway_annotation_file,
      sheet = pathway_annotation_sheet
    )
  }

  if (!all(c(pathway_col, annotation_col) %in% colnames(pathway_annotation_df))) {
    stop(
      paste0(
        "`pathway_annotation_df` must contain columns `",
        pathway_col,
        "` and `",
        annotation_col,
        "`."
      )
    )
  }

  pathway_annotation_df <- pathway_annotation_df %>%
    dplyr::transmute(
      Pathway = trimws(as.character(.data[[pathway_col]])),
      Annotation = normalize_annotation_label(.data[[annotation_col]])
    ) %>%
    dplyr::filter(!is.na(Pathway), Pathway != "") %>%
    dplyr::distinct(Pathway, .keep_all = TRUE)

  matched_annotations <- stats::setNames(
    pathway_annotation_df$Annotation,
    pathway_annotation_df$Pathway
  )
  row_annotations <- matched_annotations[pathway_names_clean]

  if (isTRUE(omit_unannotated_pathways)) {
    keep_idx <- !is.na(row_annotations) & row_annotations != "" & row_annotations != "NA"
    pathway_names_clean <- pathway_names_clean[keep_idx]
    row_annotations <- row_annotations[keep_idx]
  }

  if (!is.null(keep_annotation_group)) {
    keep_annotation_group <- trimws(as.character(keep_annotation_group))
    keep_idx <- !is.na(row_annotations) & row_annotations %in% keep_annotation_group
    pathway_names_clean <- pathway_names_clean[keep_idx]
    row_annotations <- row_annotations[keep_idx]
  }

  if (identical(filter_annotation_group, "None")) {
    return(pathway_names_clean)
  }

  if (length(filter_annotation_group) == 1 && is.na(filter_annotation_group)) {
    keep_idx <- is.na(row_annotations) | row_annotations == "" | row_annotations == "NA"
  } else {
    filter_annotation_group <- normalize_annotation_label(filter_annotation_group)
    keep_idx <- is.na(row_annotations) | !row_annotations %in% filter_annotation_group
  }

  pathway_names_clean[keep_idx]
}


prep_IPA_df <- function(dir_to_pathway,
                        pattern_int){
  
  comp_int = list.files(path = dir_to_pathway, pattern = pattern_int)
  
  path_zscore = paste0(dir_to_pathway,"/",comp_int[grepl("zscore",comp_int)])
  
  df_zscore = readxl::read_xls(path = path_zscore,skip=1)

  #test duplicates
  if(any(duplicated(df_zscore[,1]))){
    dups = as.vector(df_zscore[which(duplicated(df_zscore[,1])),][[1]])
    for (dup_gene in dups){
      df_zscore[df_zscore[,1]==dup_gene,]
      num_dup = nrow(df_zscore[df_zscore[,1]==dup_gene,])
      df_zscore[df_zscore[,1]==dup_gene,][,1][[1]] = paste0(df_zscore[df_zscore[,1]==dup_gene,][,1][[1]],"-DUP-",seq(1,num_dup))
    }
  }


  df_zscore = df_zscore %>% pivot_longer(
    cols =colnames(df_zscore)[-1],
    names_to = c("Clust_Comparison"),
    values_to = "Zscore"
  )
  
  #path_padj = paste0("IPA_results/Parent_lin_test3_Comp1/",canonical_pathway[grepl("BHpvalue",canonical_pathway)])
  path_padj = paste0(dir_to_pathway,"/",comp_int[grepl("BHpvalue",comp_int)])
  df_padj = readxl::read_xls(path = path_padj,skip=1)
  df_padj  = df_padj  %>% pivot_longer(
    cols =colnames(df_padj)[-1],
    names_to = c("Clust_Comparison"),
    values_to = "log10_padj"
  )
  
  #path_pval = paste0("IPA_results/Parent_lin_test3_Comp1/",canonical_pathway[grepl("_pvalue",canonical_pathway)])
  path_pval = paste0(dir_to_pathway,"/",comp_int[grepl("_pvalue",comp_int)])
  df_pval = readxl::read_xls(path = path_pval,skip=1)
  df_pval  = df_pval  %>% pivot_longer(
    cols =colnames(df_pval)[-1],
    names_to = c("Clust_Comparison"),
    values_to = "log10_pval"
  )
  
  if(pattern_int== "canonical_pathways"){
    colnames(df_zscore)[colnames(df_zscore) == "Canonical Pathways"] <- "Canonical_Pathways"
    colnames(df_padj)[colnames(df_padj) == "Canonical Pathways"] <- "Canonical_Pathways"
    colnames(df_pval)[colnames(df_pval) == "Canonical Pathways"] <- "Canonical_Pathways"
    merged_df = merge(
      df_zscore,
      df_padj,
      by = c("Canonical_Pathways","Clust_Comparison"),
      all = TRUE
    )
    
    merged_df = merge(
      merged_df,
      df_pval,
      by = c("Canonical_Pathways","Clust_Comparison"),
      all = TRUE
    )
  }else if(pattern_int =="diseases"){
    colnames(df_zscore)[colnames(df_zscore) == "Diseases and Bio Functions"] <- "Canonical_Pathways"
    colnames(df_padj)[colnames(df_padj) == "Diseases and Bio Functions"] <- "Canonical_Pathways"
    colnames(df_pval)[colnames(df_pval) == "Diseases and Bio Functions"] <- "Canonical_Pathways"
    
    merged_df = merge(
      df_zscore,
      df_padj,
      by = c("Canonical_Pathways","Clust_Comparison"),
      all = TRUE
    )
    
    merged_df = merge(
      merged_df,
      df_pval,
      by = c("Canonical_Pathways","Clust_Comparison"),
      all = TRUE
    )
  }
  else{
    colnames(df_zscore)[colnames(df_zscore) == "Upstream Regulators"] <- "Canonical_Pathways"
    colnames(df_padj)[colnames(df_padj) == "Upstream Regulators"] <- "Canonical_Pathways"
    colnames(df_pval)[colnames(df_pval) == "Upstream Regulators"] <- "Canonical_Pathways"
    
    merged_df = merge(
      df_zscore,
      df_padj,
      by = c("Canonical_Pathways","Clust_Comparison"),
      all = TRUE
    )
    
    merged_df = merge(
      merged_df,
      df_pval,
      by = c("Canonical_Pathways","Clust_Comparison"),
      all = TRUE
    )
  }
  
  if( sum(grepl("MT-CO1",merged_df$Clust_Comparison))>0){
    indices = grepl("MT-CO1",merged_df$Clust_Comparison)
    merged_df$Clust_Comparison[indices] = gsub("MT-CO1","MTCO1",merged_df$Clust_Comparison[indices])
  }
 
  
  merged_df$Zscore = as.numeric(merged_df$Zscore)
  merged_df$Cluster = str_split_fixed(merged_df$Clust_Comparison,"-",2)[,1]
  merged_df$Comparison = str_split_fixed(merged_df$Clust_Comparison,"-",2)[,2]
  merged_df$Condition = str_split_fixed(merged_df$Comparison, "_vs_",2)[,1]
  merged_df$Control= str_split_fixed(merged_df$Comparison, "_vs_",2)[,2]
  
  return(merged_df)
}

stars.pval <- function(x){
  stars <- c("***", "**", "*", "n.s.")
  cutoffs <- c(0, 0.001, 0.01, 0.05, 1)
  i <- findInterval(x, cutoffs, left.open = TRUE, rightmost.closed = TRUE)
  stars[i]
}

plot_heatmap_single_control_mainlin <- function(df_comp_merged,
                                                control_int,
                                                padj_or_pval_ast = "padj",
                                                filter_all = "Y",
                                                cond_int = NULL,
                                                lineage_column = NULL,
                                                lineage_colors = NULL,
                                                main_lin_order = NULL,
                                                condition_order = NULL,
                                                zscore_min = 0,
                                                padj_thresh = NULL,
                                                pvalue_thresh = NULL,
                                                zscore_disease = NULL,
                                                padj_disease = NULL,
                                                pvalue_disease = NULL,
                                                zscore_canonical = NULL,
                                                padj_canonical = NULL,
                                                pvalue_canonical = NULL,
                                                pathway_annotation_df = NULL,
                                                pathway_annotation_file = NULL,
                                                pathway_annotation_sheet = 1,
                                                pathway_col = "Pathway",
                                                annotation_col = "Annotation",
                                                annotation_order = NULL,
                                                unmatched_annotation = "Other pathways",
                                                filter_annotation_group = "None",
                                                keep_annotation_group = NULL,
                                                filter_pathway = NULL,
                                                keep_pathway = NULL,
                                                omit_unannotated_pathways = FALSE,
                                                save_pdf = TRUE,
                                                output_prefix = "IPA_single_control_mainlin_heatmap",
                                                output_dir = ".",
                                                pdf_width = 12,
                                                pdf_height = 8) {

  df_comp_merged = na.omit(df_comp_merged)

  if (!"Main_lin" %in% colnames(df_comp_merged)) {
    stop("`df_comp_merged` must contain a `Main_lin` column.")
  }

  if (!"Database" %in% colnames(df_comp_merged)) {
    stop("`df_comp_merged` must contain a `Database` column with values like `disease` and `canonical`.")
  }

  if (is.null(lineage_column)) {
    lineage_column <- if ("Main_lin_fixed" %in% colnames(df_comp_merged)) "Main_lin_fixed" else "Main_lin"
  }

  if (!lineage_column %in% colnames(df_comp_merged)) {
    stop(paste0("`df_comp_merged` must contain the lineage column `", lineage_column, "`."))
  }

  cluster_column <- if ("Cluster_fixed" %in% colnames(df_comp_merged)) "Cluster_fixed" else "Cluster"

  df_comp_merged$pval_ast = stars.pval(10^(-1 * df_comp_merged$log10_pval))
  df_comp_merged$qval_ast = stars.pval(10^(-1 * df_comp_merged$log10_padj))

  filtered_df <- df_comp_merged %>%
    dplyr::filter(Control == control_int)

  if (!is.null(cond_int)) {
    filtered_df <- filtered_df %>% dplyr::filter(Condition %in% cond_int)
  }

  if (nrow(filtered_df) == 0) {
    stop("No rows remain after filtering by `control_int` and `cond_int`.")
  }

  filtered_df$Database <- as.character(filtered_df$Database)

  if (is.null(zscore_disease)) {
    zscore_disease <- zscore_min
  }
  if (is.null(zscore_canonical)) {
    zscore_canonical <- zscore_min
  }
  if (is.null(padj_disease)) {
    padj_disease <- padj_thresh
  }
  if (is.null(padj_canonical)) {
    padj_canonical <- padj_thresh
  }
  if (is.null(pvalue_disease)) {
    pvalue_disease <- pvalue_thresh
  }
  if (is.null(pvalue_canonical)) {
    pvalue_canonical <- pvalue_thresh
  }

  keep_by_database <- function(df_in, use_any = FALSE) {
    df_disease <- df_in %>% dplyr::filter(Database == "disease")
    df_canonical <- df_in %>% dplyr::filter(Database == "canonical")
    df_other <- df_in %>% dplyr::filter(!Database %in% c("disease", "canonical"))

    filter_group <- function(df_sub, z_thr, padj_thr, pval_thr) {
      if (nrow(df_sub) == 0) {
        return(df_sub)
      }

      if (use_any) {
        if (padj_or_pval_ast == "padj") {
          if (is.null(padj_thr)) {
            stop("`padj` threshold cannot be NULL when `filter_all = \"N\"` and `padj_or_pval_ast = \"padj\"`.")
          }
          df_sub <- df_sub %>%
            dplyr::group_by(Canonical_Pathways) %>%
            dplyr::filter(any(log10_padj >= -log10(padj_thr) & abs(Zscore) > z_thr)) %>%
            dplyr::ungroup()
        } else {
          if (is.null(pval_thr)) {
            stop("`pvalue` threshold cannot be NULL when `filter_all = \"N\"` and `padj_or_pval_ast != \"padj\"`.")
          }
          df_sub <- df_sub %>%
            dplyr::group_by(Canonical_Pathways) %>%
            dplyr::filter(any(log10_pval >= -log10(pval_thr) & abs(Zscore) > z_thr)) %>%
            dplyr::ungroup()
        }
      } else {
        df_sub <- df_sub %>% dplyr::filter(abs(Zscore) > z_thr)

        if (!is.null(padj_thr)) {
          df_sub <- df_sub %>% dplyr::filter(log10_padj >= -log10(padj_thr))
        }

        if (!is.null(pval_thr)) {
          df_sub <- df_sub %>% dplyr::filter(log10_pval >= -log10(pval_thr))
        }
      }

      df_sub
    }

    dplyr::bind_rows(
      filter_group(df_disease, zscore_disease, padj_disease, pvalue_disease),
      filter_group(df_canonical, zscore_canonical, padj_canonical, pvalue_canonical),
      filter_group(df_other, zscore_min, padj_thresh, pvalue_thresh)
    )
  }

  if (filter_all == "N") {
    filtered_df <- keep_by_database(filtered_df, use_any = TRUE)
  } else {
    filtered_df <- keep_by_database(filtered_df, use_any = FALSE)
  }

  if (nrow(filtered_df) == 0) {
    stop("No rows remain after significance filtering.")
  }

  if (is.null(main_lin_order)) {
    main_lin_order <- unique(filtered_df[[lineage_column]])
  }

  if (is.null(condition_order)) {
    condition_order <- unique(filtered_df$Condition)
  }

  filtered_df[[lineage_column]] <- factor(filtered_df[[lineage_column]], levels = main_lin_order)
  filtered_df$Condition <- factor(filtered_df$Condition, levels = condition_order)

  filtered_df$plot_column <- paste(filtered_df[[cluster_column]], filtered_df$Condition, sep = "__")
  filtered_df <- filtered_df %>%
    dplyr::arrange(.data[[lineage_column]], Condition, CompNum, .data[[cluster_column]])

  col_meta <- filtered_df %>%
    dplyr::select(plot_column, dplyr::all_of(cluster_column), Condition, CompNum, dplyr::all_of(lineage_column)) %>%
    dplyr::distinct() %>%
    dplyr::arrange(.data[[lineage_column]], Condition, CompNum, .data[[cluster_column]])

  zscore_df <- filtered_df %>%
    tidyr::pivot_wider(
      id_cols = "Canonical_Pathways",
      names_from = "plot_column",
      values_from = "Zscore",
      values_fill = NA
    ) %>%
    as.data.frame()

  if (padj_or_pval_ast == "padj") {
    ast_df <- filtered_df %>%
      tidyr::pivot_wider(
        id_cols = "Canonical_Pathways",
        names_from = "plot_column",
        values_from = "qval_ast",
        values_fill = NA
      ) %>%
      as.data.frame()
  } else {
    ast_df <- filtered_df %>%
      tidyr::pivot_wider(
        id_cols = "Canonical_Pathways",
        names_from = "plot_column",
        values_from = "pval_ast",
        values_fill = NA
      ) %>%
      as.data.frame()
  }

  zscore_mat <- as.matrix(zscore_df[, -1, drop = FALSE])
  rownames(zscore_mat) <- zscore_df$Canonical_Pathways

  ast_mat <- as.matrix(ast_df[, -1, drop = FALSE])
  rownames(ast_mat) <- ast_df$Canonical_Pathways

  ordered_columns <- col_meta$plot_column
  zscore_mat <- zscore_mat[, ordered_columns, drop = FALSE]
  ast_mat <- ast_mat[, ordered_columns, drop = FALSE]

  pathways_to_keep <- filter_pathways_by_annotation(
    pathway_names = rownames(zscore_mat),
    pathway_annotation_df = pathway_annotation_df,
    pathway_annotation_file = pathway_annotation_file,
    pathway_annotation_sheet = pathway_annotation_sheet,
    pathway_col = pathway_col,
    annotation_col = annotation_col,
    filter_annotation_group = filter_annotation_group,
    keep_annotation_group = keep_annotation_group,
    filter_pathway = filter_pathway,
    keep_pathway = keep_pathway,
    omit_unannotated_pathways = omit_unannotated_pathways
  )

  zscore_mat <- zscore_mat[rownames(zscore_mat) %in% pathways_to_keep, , drop = FALSE]
  ast_mat <- ast_mat[rownames(ast_mat) %in% pathways_to_keep, , drop = FALSE]

  if (nrow(zscore_mat) == 0) {
    stop("No pathways remain after applying annotation-group and/or pathway filters.")
  }

  row_split_info <- prepare_pathway_row_split(
    pathway_names = rownames(zscore_mat),
    pathway_annotation_df = pathway_annotation_df,
    pathway_annotation_file = pathway_annotation_file,
    pathway_annotation_sheet = pathway_annotation_sheet,
    pathway_col = pathway_col,
    annotation_col = annotation_col,
    annotation_order = annotation_order,
    unmatched_annotation = unmatched_annotation
  )

  zscore_mat <- zscore_mat[row_split_info$row_order, , drop = FALSE]
  ast_mat <- ast_mat[row_split_info$row_order, , drop = FALSE]

  lineages <- as.character(col_meta[[lineage_column]])
  conditions <- as.character(col_meta$Condition)
  cluster_labels <- as.character(col_meta[[cluster_column]])

  condition_colors <- c(
    "TermPE" = "#111111",
    "PrePE" = "#BDBDBD",
    "FGR" = "#D55E00",
    "PTB" = "#CC79A7",
    "DM1" = "#0072B2"
  )
  missing_conditions <- setdiff(unique(conditions), names(condition_colors))
  if (length(missing_conditions) > 0) {
    extra_cols <- grDevices::hcl.colors(length(missing_conditions), "Dark 3")
    names(extra_cols) <- missing_conditions
    condition_colors <- c(condition_colors, extra_cols)
  }

  if (is.null(lineage_colors)) {
    lineage_colors <- c(
      "Myeloid" = "#590902",
      "HBC" = "#BF3604",
      "Trophoblasts" = "#60B364",
      "Stromal" = "#022873",
      "Endothelial" = "#6B0AC9",
      "Monocytes" = "#590902",
      "PAMM" = "#590902",
      "DC" = "#590902",
      "Troph" = "#60B364",
      "EndoF" = "#6B0AC9",
      "EndoM" = "#6B0AC9",
      "dS" = "#022873",
      "fFB" = "#022873",
      "fPV" = "#022873",
      "NK" = "#7570B3",
      "Tcells" = "#1F78B4",
      "Bcells" = "#33A02C",
      "-" = "#999999"
    )
  }
  missing_lineages <- setdiff(unique(lineages), names(lineage_colors))
  if (length(missing_lineages) > 0) {
    extra_cols <- grDevices::hcl.colors(length(missing_lineages), "Set 2")
    names(extra_cols) <- missing_lineages
    lineage_colors <- c(lineage_colors, extra_cols)
  }

  if (length(unique(stats::na.omit(conditions))) <= 1) {
    top_annot <- ComplexHeatmap::HeatmapAnnotation(
      Lineage = lineages,
      col = list(
        Lineage = lineage_colors
      ),
      show_annotation_name = TRUE,
      annotation_name_gp = grid::gpar(fontface = "bold", fontfamily = heatmap_font_family),
      annotation_legend_param = list(
        Lineage = list(title_position = "topcenter", nrow = 1)
      ),
      border = TRUE
    )
  } else {
    top_annot <- ComplexHeatmap::HeatmapAnnotation(
      Lineage = lineages,
      Complication = conditions,
      col = list(
        Lineage = lineage_colors,
        Complication = condition_colors
      ),
      show_annotation_name = TRUE,
      annotation_name_gp = grid::gpar(fontface = "bold", fontfamily = heatmap_font_family),
      annotation_legend_param = list(
        Lineage = list(title_position = "topcenter", nrow = 1),
        Complication = list(title_position = "topcenter", nrow = 1)
      ),
      border = TRUE
    )
  }

  zscore_limits <- range(zscore_mat, na.rm = TRUE)
  zscore_col <- circlize::colorRamp2(
    c(zscore_limits[1], 0, zscore_limits[2]),
    c("#313695", "#FFFFFF", "#A50026")
  )

  ht <- ComplexHeatmap::Heatmap(
    zscore_mat,
    name = "Zscore",
    heatmap_legend_param = list(direction = "horizontal", title_position = "topcenter"),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    column_split = factor(lineages, levels = main_lin_order),
    column_labels = cluster_labels,
    row_names_gp = grid::gpar(fontfamily = heatmap_font_family),
    column_names_gp = grid::gpar(fontfamily = heatmap_font_family),
    top_annotation = top_annot,
    border_gp = grid::gpar(col = "black", lwd = 1),
    rect_gp = grid::gpar(col = "grey30", lwd = 0.5),
    width = ncol(zscore_mat) * grid::unit(5, "mm"),
    height = nrow(zscore_mat) * grid::unit(5, "mm"),
    row_names_side = "left",
    row_split = row_split_info$row_split,
    row_title_side = "right",
    row_title_rot = 0,
    row_title_gp = grid::gpar(fontface = "bold", fontfamily = heatmap_font_family),
    column_title = NULL,
    col = zscore_col,
    na_col = "white",
    cell_fun = function(j, i, x, y, w, h, fill) {
      if (!is.na(ast_mat[i, j]) && ast_mat[i, j] != "n.s.") {
        grid::grid.text(ast_mat[i, j], x, y, gp = grid::gpar(fontfamily = heatmap_font_family))
      }
    }
  )

  if (save_pdf) {
    open_heatmap_pdf_device(
      file = file.path(output_dir, paste0(output_prefix, ".pdf")),
      width = pdf_width,
      height = pdf_height
    )
    ComplexHeatmap::draw(
      ht,
      annotation_legend_side = "bottom",
      heatmap_legend_side = "bottom",
      merge_legend = TRUE,
      column_title = paste0("Cell subset - complication (Control: ", control_int, ")"),
      column_title_side = "bottom",
      column_title_gp = grid::gpar(fontface = "bold", fontfamily = heatmap_font_family)
    )
    dev.off()
    message("Saved heatmap to PDF.")
  } else {
    ComplexHeatmap::draw(
      ht,
      annotation_legend_side = "bottom",
      heatmap_legend_side = "bottom",
      merge_legend = TRUE,
      column_title = paste0("Cell subset - complication (Control: ", control_int, ")"),
      column_title_side = "bottom",
      column_title_gp = grid::gpar(fontface = "bold", fontfamily = heatmap_font_family)
    )
  }
}


plot_heatmap_single_condition_mainlin <- function(df_comp_merged,
                                                  cond_int,
                                                  control_int = NULL,
                                                  padj_or_pval_ast = "padj",
                                                  filter_all = "Y",
                                                  lineage_column = NULL,
                                                  lineage_colors = NULL,
                                                  main_lin_order = NULL,
                                                  control_order = NULL,
                                                  zscore_min = 0,
                                                  padj_thresh = NULL,
                                                  pvalue_thresh = NULL,
                                                  zscore_disease = NULL,
                                                  padj_disease = NULL,
                                                  pvalue_disease = NULL,
                                                  zscore_canonical = NULL,
                                                  padj_canonical = NULL,
                                                  pvalue_canonical = NULL,
                                                  pathway_annotation_df = NULL,
                                                  pathway_annotation_file = NULL,
                                                  pathway_annotation_sheet = 1,
                                                  pathway_col = "Pathway",
                                                  annotation_col = "Annotation",
                                                  annotation_order = NULL,
                                                  unmatched_annotation = "Other pathways",
                                                  filter_annotation_group = "None",
                                                  keep_annotation_group = NULL,
                                                  filter_pathway = NULL,
                                                  keep_pathway = NULL,
                                                  omit_unannotated_pathways = FALSE,
                                                  save_pdf = TRUE,
                                                  output_prefix = "IPA_single_condition_mainlin_heatmap",
                                                  output_dir = ".",
                                                  pdf_width = 12,
                                                  pdf_height = 8) {

  df_comp_merged = na.omit(df_comp_merged)

  if (!"Main_lin" %in% colnames(df_comp_merged)) {
    stop("`df_comp_merged` must contain a `Main_lin` column.")
  }

  if (!"Database" %in% colnames(df_comp_merged)) {
    stop("`df_comp_merged` must contain a `Database` column with values like `disease` and `canonical`.")
  }

  if (is.null(lineage_column)) {
    lineage_column <- if ("Main_lin_fixed" %in% colnames(df_comp_merged)) "Main_lin_fixed" else "Main_lin"
  }

  if (!lineage_column %in% colnames(df_comp_merged)) {
    stop(paste0("`df_comp_merged` must contain the lineage column `", lineage_column, "`."))
  }

  cluster_column <- if ("Cluster_fixed" %in% colnames(df_comp_merged)) "Cluster_fixed" else "Cluster"

  df_comp_merged$pval_ast = stars.pval(10^(-1 * df_comp_merged$log10_pval))
  df_comp_merged$qval_ast = stars.pval(10^(-1 * df_comp_merged$log10_padj))

  filtered_df <- df_comp_merged %>%
    dplyr::filter(Condition %in% cond_int)

  if (!is.null(control_int)) {
    filtered_df <- filtered_df %>% dplyr::filter(Control %in% control_int)
  }

  if (nrow(filtered_df) == 0) {
    stop("No rows remain after filtering by `cond_int` and `control_int`.")
  }

  filtered_df$Database <- as.character(filtered_df$Database)

  if (is.null(zscore_disease)) {
    zscore_disease <- zscore_min
  }
  if (is.null(zscore_canonical)) {
    zscore_canonical <- zscore_min
  }
  if (is.null(padj_disease)) {
    padj_disease <- padj_thresh
  }
  if (is.null(padj_canonical)) {
    padj_canonical <- padj_thresh
  }
  if (is.null(pvalue_disease)) {
    pvalue_disease <- pvalue_thresh
  }
  if (is.null(pvalue_canonical)) {
    pvalue_canonical <- pvalue_thresh
  }

  keep_by_database <- function(df_in, use_any = FALSE) {
    df_disease <- df_in %>% dplyr::filter(Database == "disease")
    df_canonical <- df_in %>% dplyr::filter(Database == "canonical")
    df_other <- df_in %>% dplyr::filter(!Database %in% c("disease", "canonical"))

    filter_group <- function(df_sub, z_thr, padj_thr, pval_thr) {
      if (nrow(df_sub) == 0) {
        return(df_sub)
      }

      if (use_any) {
        if (padj_or_pval_ast == "padj") {
          if (is.null(padj_thr)) {
            stop("`padj` threshold cannot be NULL when `filter_all = \"N\"` and `padj_or_pval_ast = \"padj\"`.")
          }
          df_sub <- df_sub %>%
            dplyr::group_by(Canonical_Pathways) %>%
            dplyr::filter(any(log10_padj >= -log10(padj_thr) & abs(Zscore) > z_thr)) %>%
            dplyr::ungroup()
        } else {
          if (is.null(pval_thr)) {
            stop("`pvalue` threshold cannot be NULL when `filter_all = \"N\"` and `padj_or_pval_ast != \"padj\"`.")
          }
          df_sub <- df_sub %>%
            dplyr::group_by(Canonical_Pathways) %>%
            dplyr::filter(any(log10_pval >= -log10(pval_thr) & abs(Zscore) > z_thr)) %>%
            dplyr::ungroup()
        }
      } else {
        df_sub <- df_sub %>% dplyr::filter(abs(Zscore) > z_thr)

        if (!is.null(padj_thr)) {
          df_sub <- df_sub %>% dplyr::filter(log10_padj >= -log10(padj_thr))
        }

        if (!is.null(pval_thr)) {
          df_sub <- df_sub %>% dplyr::filter(log10_pval >= -log10(pval_thr))
        }
      }

      df_sub
    }

    dplyr::bind_rows(
      filter_group(df_disease, zscore_disease, padj_disease, pvalue_disease),
      filter_group(df_canonical, zscore_canonical, padj_canonical, pvalue_canonical),
      filter_group(df_other, zscore_min, padj_thresh, pvalue_thresh)
    )
  }

  if (filter_all == "N") {
    filtered_df <- keep_by_database(filtered_df, use_any = TRUE)
  } else {
    filtered_df <- keep_by_database(filtered_df, use_any = FALSE)
  }

  if (nrow(filtered_df) == 0) {
    stop("No rows remain after significance filtering.")
  }

  if (is.null(main_lin_order)) {
    main_lin_order <- unique(filtered_df[[lineage_column]])
  }

  if (is.null(control_order)) {
    control_order <- unique(filtered_df$Control)
  }

  filtered_df[[lineage_column]] <- factor(filtered_df[[lineage_column]], levels = main_lin_order)
  filtered_df$Control <- factor(filtered_df$Control, levels = control_order)

  filtered_df$plot_column <- paste(filtered_df[[cluster_column]], filtered_df$Control, sep = "__")
  filtered_df <- filtered_df %>%
    dplyr::arrange(.data[[lineage_column]], Control, CompNum, .data[[cluster_column]])

  col_meta <- filtered_df %>%
    dplyr::select(plot_column, dplyr::all_of(cluster_column), Control, CompNum, dplyr::all_of(lineage_column)) %>%
    dplyr::distinct() %>%
    dplyr::arrange(.data[[lineage_column]], Control, CompNum, .data[[cluster_column]])

  zscore_df <- filtered_df %>%
    tidyr::pivot_wider(
      id_cols = "Canonical_Pathways",
      names_from = "plot_column",
      values_from = "Zscore",
      values_fill = NA
    ) %>%
    as.data.frame()

  if (padj_or_pval_ast == "padj") {
    ast_df <- filtered_df %>%
      tidyr::pivot_wider(
        id_cols = "Canonical_Pathways",
        names_from = "plot_column",
        values_from = "qval_ast",
        values_fill = NA
      ) %>%
      as.data.frame()
  } else {
    ast_df <- filtered_df %>%
      tidyr::pivot_wider(
        id_cols = "Canonical_Pathways",
        names_from = "plot_column",
        values_from = "pval_ast",
        values_fill = NA
      ) %>%
      as.data.frame()
  }

  zscore_mat <- as.matrix(zscore_df[, -1, drop = FALSE])
  rownames(zscore_mat) <- zscore_df$Canonical_Pathways

  ast_mat <- as.matrix(ast_df[, -1, drop = FALSE])
  rownames(ast_mat) <- ast_df$Canonical_Pathways

  ordered_columns <- col_meta$plot_column
  zscore_mat <- zscore_mat[, ordered_columns, drop = FALSE]
  ast_mat <- ast_mat[, ordered_columns, drop = FALSE]

  pathways_to_keep <- filter_pathways_by_annotation(
    pathway_names = rownames(zscore_mat),
    pathway_annotation_df = pathway_annotation_df,
    pathway_annotation_file = pathway_annotation_file,
    pathway_annotation_sheet = pathway_annotation_sheet,
    pathway_col = pathway_col,
    annotation_col = annotation_col,
    filter_annotation_group = filter_annotation_group,
    keep_annotation_group = keep_annotation_group,
    filter_pathway = filter_pathway,
    keep_pathway = keep_pathway,
    omit_unannotated_pathways = omit_unannotated_pathways
  )

  zscore_mat <- zscore_mat[rownames(zscore_mat) %in% pathways_to_keep, , drop = FALSE]
  ast_mat <- ast_mat[rownames(ast_mat) %in% pathways_to_keep, , drop = FALSE]

  if (nrow(zscore_mat) == 0) {
    stop("No pathways remain after applying annotation-group and/or pathway filters.")
  }

  row_split_info <- prepare_pathway_row_split(
    pathway_names = rownames(zscore_mat),
    pathway_annotation_df = pathway_annotation_df,
    pathway_annotation_file = pathway_annotation_file,
    pathway_annotation_sheet = pathway_annotation_sheet,
    pathway_col = pathway_col,
    annotation_col = annotation_col,
    annotation_order = annotation_order,
    unmatched_annotation = unmatched_annotation
  )

  zscore_mat <- zscore_mat[row_split_info$row_order, , drop = FALSE]
  ast_mat <- ast_mat[row_split_info$row_order, , drop = FALSE]

  lineages <- as.character(col_meta[[lineage_column]])
  controls <- as.character(col_meta$Control)
  cluster_labels <- as.character(col_meta[[cluster_column]])

  control_colors <- c(
    "Cont_labored" = "#2F2F2F",
    "Cont_unlabored" = "#4A5A6A",
    "Cont_PTB" = "#CD534CFF"
  )
  missing_controls <- setdiff(unique(controls), names(control_colors))
  if (length(missing_controls) > 0) {
    extra_cols <- grDevices::hcl.colors(length(missing_controls), "Dark 3")
    names(extra_cols) <- missing_controls
    control_colors <- c(control_colors, extra_cols)
  }

  if (is.null(lineage_colors)) {
    lineage_colors <- c(
      "Myeloid" = "#590902",
      "HBC" = "#BF3604",
      "Trophoblasts" = "#60B364",
      "Stromal" = "#022873",
      "Endothelial" = "#6B0AC9",
      "Monocytes" = "#590902",
      "PAMM" = "#590902",
      "DC" = "#590902",
      "Troph" = "#60B364",
      "EndoF" = "#6B0AC9",
      "EndoM" = "#6B0AC9",
      "dS" = "#022873",
      "fFB" = "#022873",
      "fPV" = "#022873",
      "NK" = "#7570B3",
      "Tcells" = "#1F78B4",
      "Bcells" = "#33A02C",
      "-" = "#999999"
    )
  }
  missing_lineages <- setdiff(unique(lineages), names(lineage_colors))
  if (length(missing_lineages) > 0) {
    extra_cols <- grDevices::hcl.colors(length(missing_lineages), "Set 2")
    names(extra_cols) <- missing_lineages
    lineage_colors <- c(lineage_colors, extra_cols)
  }

  if (length(unique(stats::na.omit(controls))) <= 1) {
    top_annot <- ComplexHeatmap::HeatmapAnnotation(
      Lineage = lineages,
      col = list(
        Lineage = lineage_colors
      ),
      show_annotation_name = TRUE,
      annotation_name_gp = grid::gpar(fontface = "bold", fontfamily = heatmap_font_family),
      annotation_legend_param = list(
        Lineage = list(title_position = "topcenter", nrow = 1)
      ),
      border = TRUE
    )
  } else {
    top_annot <- ComplexHeatmap::HeatmapAnnotation(
      Lineage = lineages,
      Control = controls,
      col = list(
        Lineage = lineage_colors,
        Control = control_colors
      ),
      show_annotation_name = TRUE,
      annotation_name_gp = grid::gpar(fontface = "bold", fontfamily = heatmap_font_family),
      annotation_legend_param = list(
        Lineage = list(title_position = "topcenter", nrow = 1),
        Control = list(title_position = "topcenter", nrow = 1)
      ),
      border = TRUE
    )
  }

  zscore_limits <- range(zscore_mat, na.rm = TRUE)
  zscore_col <- circlize::colorRamp2(
    c(zscore_limits[1], 0, zscore_limits[2]),
    c("#313695", "#FFFFFF", "#A50026")
  )

  ht <- ComplexHeatmap::Heatmap(
    zscore_mat,
    name = "Zscore",
    heatmap_legend_param = list(direction = "horizontal", title_position = "topcenter"),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    column_split = factor(lineages, levels = main_lin_order),
    column_labels = cluster_labels,
    row_names_gp = grid::gpar(fontfamily = heatmap_font_family),
    column_names_gp = grid::gpar(fontfamily = heatmap_font_family),
    top_annotation = top_annot,
    border_gp = grid::gpar(col = "black", lwd = 1),
    rect_gp = grid::gpar(col = "grey30", lwd = 0.5),
    width = ncol(zscore_mat) * grid::unit(5, "mm"),
    height = nrow(zscore_mat) * grid::unit(5, "mm"),
    row_names_side = "left",
    row_split = row_split_info$row_split,
    row_title_side = "right",
    row_title_rot = 0,
    row_title_gp = grid::gpar(fontface = "bold", fontfamily = heatmap_font_family),
    column_title = NULL,
    col = zscore_col,
    na_col = "white",
    cell_fun = function(j, i, x, y, w, h, fill) {
      if (!is.na(ast_mat[i, j]) && ast_mat[i, j] != "n.s.") {
        grid::grid.text(ast_mat[i, j], x, y, gp = grid::gpar(fontfamily = heatmap_font_family))
      }
    }
  )

  cond_label <- paste(unique(as.character(filtered_df$Condition)), collapse = ", ")

  if (save_pdf) {
    open_heatmap_pdf_device(
      file = file.path(output_dir, paste0(output_prefix, ".pdf")),
      width = pdf_width,
      height = pdf_height
    )
    ComplexHeatmap::draw(
      ht,
      annotation_legend_side = "bottom",
      heatmap_legend_side = "bottom",
      merge_legend = TRUE,
      column_title = paste0("Cell subset - control (Condition: ", cond_label, ")"),
      column_title_side = "bottom",
      column_title_gp = grid::gpar(fontface = "bold", fontfamily = heatmap_font_family)
    )
    dev.off()
    message("Saved heatmap to PDF.")
  } else {
    ComplexHeatmap::draw(
      ht,
      annotation_legend_side = "bottom",
      heatmap_legend_side = "bottom",
      merge_legend = TRUE,
      column_title = paste0("Cell subset - control (Condition: ", cond_label, ")"),
      column_title_side = "bottom",
      column_title_gp = grid::gpar(fontface = "bold", fontfamily = heatmap_font_family)
    )
  }
}



  
