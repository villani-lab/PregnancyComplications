suppressPackageStartupMessages({
  library(openxlsx)
  library(readxl)
  library(tidyr)
  library(stringr)
  library(ComplexHeatmap)
  library(colorRamp2)
  library(dplyr)
  library(circlize)
  library(grid) # for adjusting matrix box size
  })
# Functions----- (still optimizing)
read_pathway_genes = function(cell_type,
                              condition,
                              control,
                              loc_of_IPA_gene_files,
                              diseases_or_pathway = "diseases"){

  file_name = paste0(cell_type,"-",condition,"_vs_",
                            control,"-",
                     diseases_or_pathway,".xls")

  file_int = read_excel(col_names = TRUE,
                        skip = 1,
                        path = paste0(loc_of_IPA_gene_files,"/",file_name))


  return(file_int)
}

# NOTE: this version takes an already-loaded IPA results data frame (df_pathway_genes, "diseases" or canonical
# pathway table selected with diseases_or_pathway). A different function with a similar purpose that reads the
# IPA export files from disk itself is get_log2fc_of_pathway_genes_from_IPA_files() in figure_functions.R;
# the two are intentionally kept separate and have different arguments.
get_log2fc_of_pathway_genes_from_df = function(path_int, #pathway of interest
                                       df_pathway_genes, #IPA results, has pathway info and molecules
                                       cell_type,
                                       condition,
                                       control,
                                       df_loc, #contains the info as to where the files are
                                        diseases_or_pathway = "diseases"
                                       ){
  if(diseases_or_pathway=="diseases"){
    df_path_int = subset(df_pathway_genes,`Diseases or Functions Annotation`== path_int)
  }
  else{
    df_path_int = subset(df_pathway_genes,`Ingenuity Canonical Pathways`== path_int)
  }


  # get the DEGs of the comparison of interest of the cell
  mask = (df_loc$Cell_type==cell_type & df_loc$name == paste0(condition,"_vs_",control))
  loc_int = df_loc[mask,]

  info_on_loc = str_split(loc_int$name_ipa,"_",4)[[1]]
  if(info_on_loc[[1]][1]=="Parent"){
    csv_file = paste0("/projects/Pregnancy_placenta/Results/13_DE_PAL/output_csv_IPA/redo_",
              info_on_loc[[1]][1],"lin/",loc_int$name_ipa,".csv")
  }
  else{csv_file = paste0("/projects/Pregnancy_placenta/Results/13_DE_PAL/output_csv_IPA/redo_",
              info_on_loc[[1]][1],"/",loc_int$name_ipa,".csv")}


  log2fc_gene_df = read.csv(csv_file)
  vals_int = log2fc_gene_df[,grepl(paste0(cell_type,".",condition,"_vs_",control),colnames(log2fc_gene_df))]
  rownames(vals_int) = log2fc_gene_df$GeneName

  #only keep the genes that are in the pathway
  vals_int = vals_int[which(rownames(vals_int) %in% str_split(df_path_int$Molecules,",")[[1]]),]
  vals_int = data.frame(vals_int)
  colnames(vals_int) = c("log2fc","padj")
  vals_int$ast = stars.pval(1*vals_int$padj)

  vals_int$gene_name = rownames(vals_int)
  rownames(vals_int) = NULL

  return(vals_int)
}


stars.pval <- function(x){
  stars <- c("***", "**", "*", "n.s.")
  cutoffs <- c(0, 0.001, 0.01, 0.05, 1)
  i <- findInterval(x, cutoffs, left.open = TRUE, rightmost.closed = TRUE)
  stars[i]
}


get_pathway_gene_heatmap_mats = function(path_of_interest,
                                         list_files,
                                         df_loc,
                                         loc_of_IPA_gene_files,
                                         file_regex = "PrePE|TermPE",
                                         disease_regex = "-diseases.xls",
                                         cell_types_keep = c("SCT", "VCTccc", "VCTfusion_2", "VCT_1"),
                                         comparison_regex = "_Cont_unlabored",
                                         log2fc_filter = NULL,
                                         diseases_or_pathway = "diseases") {

  df = data.frame(
    full_name = list_files[
      grepl(file_regex, list_files) & grepl(disease_regex, list_files)
    ],
    stringsAsFactors = FALSE
  )

  if (!is.null(cell_types_keep)) {
    cell_pattern = paste(cell_types_keep, collapse = "|")
    df = df[grepl(cell_pattern, df$full_name), , drop = FALSE]
  }

  if (nrow(df) == 0) {
    stop("No IPA files matched the requested filters.")
  }

  cell_type = character(nrow(df))
  condition = character(nrow(df))
  control = character(nrow(df))

  for (i in seq_len(nrow(df))) {
    temp_cell = str_split_fixed(df$full_name[i], "-", 3)[1]
    temp_comparison = str_split_fixed(df$full_name[i], "-", 3)[2]
    temp_condition = str_split_fixed(temp_comparison, "_vs_", 2)[1]
    temp_control = str_split_fixed(temp_comparison, "_vs_", 2)[2]

    cell_type[i] = temp_cell
    condition[i] = temp_condition
    control[i] = temp_control
  }

  df$cell_type = cell_type
  df$condition = condition
  df$control = control
  df$pathway_in = "N"
  df$pathway_score = NA_real_

  pathway_gene_list = list()

  for (iter in seq_len(nrow(df))) {
    temp_all_pathways_cond_int = read_pathway_genes(
      cell_type = df$cell_type[iter],
      condition = df$condition[iter],
      control = df$control[iter],
      loc_of_IPA_gene_files = loc_of_IPA_gene_files,
      diseases_or_pathway = diseases_or_pathway
    )
    if(diseases_or_pathway=="diseases"){
      pathway_hits = temp_all_pathways_cond_int$`Diseases or Functions Annotation` == path_of_interest
    }
    else{ # assuming canonical pathway
      pathway_hits = temp_all_pathways_cond_int$`Ingenuity Canonical Pathways` == path_of_interest
    }


    if (any(pathway_hits, na.rm = TRUE)) {
      df$pathway_in[iter] = "Y"
      if(diseases_or_pathway=="diseases"){
        df$pathway_score[iter] = temp_all_pathways_cond_int[pathway_hits, ]$`Activation z-score`[1]
      }
      else{
        df$pathway_score[iter] = temp_all_pathways_cond_int[pathway_hits, ]$`z-score`[1]
      }
      log2_df = get_log2fc_of_pathway_genes_from_df(
        path_int = path_of_interest,
        df_pathway_genes = temp_all_pathways_cond_int,
        cell_type = df$cell_type[iter],
        condition = df$condition[iter],
        control = df$control[iter],
        df_loc = df_loc,
        diseases_or_pathway = diseases_or_pathway
      )

      comparison_name = gsub("-diseases.xls", "", df$full_name[iter])
      log2_df$comparison = comparison_name
      pathway_gene_list[[comparison_name]] = log2_df
    }
  }

  pathway_gene_df = bind_rows(pathway_gene_list)

  if (nrow(pathway_gene_df) == 0) {
    stop(paste0("Pathway not found in matched files: ", path_of_interest))
  }

  if (!is.null(comparison_regex)) {
    pathway_gene_df = pathway_gene_df[
      grepl(comparison_regex, pathway_gene_df$comparison),
      ,
      drop = FALSE
    ]
  }

  if (!is.null(log2fc_filter)) {
    if (identical(log2fc_filter, "positive")) {
      pathway_gene_df = pathway_gene_df[pathway_gene_df$log2fc > 0, , drop = FALSE]
    } else if (identical(log2fc_filter, "negative")) {
      pathway_gene_df = pathway_gene_df[pathway_gene_df$log2fc < 0, , drop = FALSE]
    } else {
      stop("log2fc_filter must be NULL, 'positive', or 'negative'.")
    }
  }

  if (nrow(pathway_gene_df) == 0) {
    stop("No rows remained after applying comparison/log2fc filters.")
  }

  log2fc_df = as.data.frame(
    pathway_gene_df %>%
      pivot_wider(
        id_cols = "gene_name",
        names_from = "comparison",
        values_from = "log2fc",
        values_fill = NA
      )
  )

  ast_df = as.data.frame(
    pathway_gene_df %>%
      pivot_wider(
        id_cols = "gene_name",
        names_from = "comparison",
        values_from = "ast",
        values_fill = NA
      )
  )

  log2fc_mat = as.matrix(log2fc_df[, -1, drop = FALSE])
  ast_mat = as.matrix(ast_df[, -1, drop = FALSE])
  rownames(log2fc_mat) = log2fc_df$gene_name
  rownames(ast_mat) = ast_df$gene_name

  return(list(
    log2fc_mat = log2fc_mat,
    ast_mat = ast_mat,
    pathway_gene_df = pathway_gene_df,
    metadata = df
  ))
}


plot_pathway_gene_heatmap = function(log2fc_mat,
                                     ast_mat,
                                     path_of_interest,
                                     column_title = NULL,
                                     annotation_name = "Cond",
                                     annotation_values = NULL,
                                     annotation_colors = c(
                                       "TermPE" = "#9ecae1",
                                       "PrePE" = "#fdae6b"
                                     ),
                                     column_labels = NULL,
                                     column_label_regex = NULL,
                                     na_col = "grey90",
                                     cluster_columns = FALSE,
                                     cluster_rows = FALSE,
                                     cell_w = 4,
                                     cell_h = 4) {

  if (!all(dim(log2fc_mat) == dim(ast_mat))) {
    stop("log2fc_mat and ast_mat must have the same dimensions.")
  }

  if (!identical(rownames(log2fc_mat), rownames(ast_mat))) {
    stop("log2fc_mat and ast_mat must have the same rownames.")
  }

  if (!identical(colnames(log2fc_mat), colnames(ast_mat))) {
    stop("log2fc_mat and ast_mat must have the same colnames.")
  }

  if (is.null(column_labels)) {
    column_labels = colnames(log2fc_mat)
    if (!is.null(column_label_regex)) {
      column_labels = sub(column_label_regex, "", column_labels)
    }
  }

  if (is.null(annotation_values)) {
    annotation_values = ifelse(grepl("TermPE", colnames(log2fc_mat)), "TermPE", "PrePE")
  }

  if (length(annotation_values) != ncol(log2fc_mat)) {
    stop("annotation_values must have one entry per heatmap column.")
  }

  annotation_df = data.frame(annotation = annotation_values, stringsAsFactors = FALSE)
  colnames(annotation_df) = annotation_name

  top_anno = HeatmapAnnotation(
    df = annotation_df,
    col = setNames(list(annotation_colors), annotation_name),
    border = TRUE
  )

  min_val = min(log2fc_mat, na.rm = TRUE)
  max_val = max(log2fc_mat, na.rm = TRUE)
  col_fun = colorRamp2(
    c(min_val, 0, max_val),
    c("#2B3CFF", "#FFFEFE", "#BB0103")
  )

  if (is.null(column_title)) {
    column_title = path_of_interest
  }

  Heatmap(
    log2fc_mat,
    column_labels = column_labels,
    name = "log2FC",
    col = col_fun,
    na_col = na_col,
    show_row_names = TRUE,
    top_annotation = top_anno,
    row_names_side = "left",
    column_title = column_title,
    cluster_columns = cluster_columns,
    cluster_rows = cluster_rows,
    rect_gp = grid::gpar(col = "black", lwd = 0.5),
    width = ncol(log2fc_mat) * unit(cell_w, "mm"),
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
}

# upset plot _ TEMP

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(cowplot)
})

pretty_comparison <- function(x) {
  gsub("_vs_", " vs ", x, fixed = TRUE)
}

match_requested_comparisons <- function(available_comparisons, requested_comparisons) {
  if (length(requested_comparisons) == 1 && identical(tolower(requested_comparisons), "all")) {
    return(available_comparisons)
  }

  available_labels <- pretty_comparison(available_comparisons)
  names(available_labels) <- available_comparisons

  matched <- vapply(
    requested_comparisons,
    FUN.VALUE = character(1),
    FUN = function(x) {
      if (x %in% available_comparisons) {
        return(x)
      }

      label_match <- available_comparisons[available_labels == x]
      if (length(label_match) > 0) {
        return(label_match[1])
      }

      stop(
        "Comparison not found: ", x,
        "\nAvailable comparisons are:\n",
        paste(c(available_comparisons, available_labels), collapse = "\n")
      )
    }
  )

  unique(matched)
}

make_gene_overlap_upset <- function(
  csv_file,
  comparisons = "all",
  comparison_col = "comparison",
  gene_col = "gene_name",
  logfc_col = "log2fc",
  positive_color = "black",
  negative_color = "red",
  save_pdf = TRUE,
  output_dir = getwd(),
  output_file = NULL,
  top_n_intersections = 20,
  width = 16,
  height = 8,
  title = "Gene overlap across comparison groups",
  subtitle = NULL,
  base_size = 11,
  point_color = "#12344D",
  line_color = "#2F6C8F",
  absent_point_color = "grey85",
  return_data = FALSE
) {
  if (!file.exists(csv_file)) {
    stop("Input file not found: ", csv_file)
  }

  dat <- read.csv(csv_file, stringsAsFactors = FALSE, check.names = FALSE)
  dat <- dat[, !is.na(names(dat)) & names(dat) != "", drop = FALSE]

  required_cols <- c(comparison_col, gene_col, logfc_col)
  missing_cols <- setdiff(required_cols, names(dat))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  dat <- dat %>%
    transmute(
      comparison = trimws(.data[[comparison_col]]),
      gene_name = trimws(.data[[gene_col]]),
      log2fc = as.numeric(.data[[logfc_col]])
    ) %>%
    filter(comparison != "", gene_name != "")

  selected_comparisons <- match_requested_comparisons(
    available_comparisons = unique(dat$comparison),
    requested_comparisons = comparisons
  )

  dat <- dat %>%
    filter(comparison %in% selected_comparisons)

  if (nrow(dat) == 0) {
    stop("No rows remain after filtering the requested comparisons.")
  }

  gene_direction <- dat %>%
    group_by(comparison, gene_name) %>%
    summarise(log2fc = dplyr::first(log2fc), .groups = "drop") %>%
    mutate(
      direction = ifelse(log2fc >= 0, "Positive log2FC", "Negative log2FC")
    )

  membership_long <- gene_direction %>%
    select(comparison, gene_name) %>%
    distinct(comparison, gene_name)

  set_sizes_direction <- gene_direction %>%
    count(comparison, direction, name = "set_size") %>%
    complete(
      comparison,
      direction = c("Positive log2FC", "Negative log2FC"),
      fill = list(set_size = 0L)
    )

  set_sizes <- set_sizes_direction %>%
    group_by(comparison) %>%
    summarise(set_size = sum(set_size), .groups = "drop") %>%
    arrange(desc(set_size), comparison)

  comparison_levels <- set_sizes$comparison
  comparison_labels <- pretty_comparison(comparison_levels)
  names(comparison_labels) <- comparison_levels

  comparison_factor_levels <- rev(comparison_levels)
  comparison_label_levels <- rev(unname(comparison_labels))

  set_sizes_direction <- set_sizes_direction %>%
    mutate(
      comparison = factor(comparison, levels = comparison_factor_levels),
      comparison_label = factor(
        comparison_labels[as.character(comparison)],
        levels = comparison_label_levels
      ),
      direction = factor(
        direction,
        levels = c("Negative log2FC", "Positive log2FC")
      )
    )

  membership_wide <- membership_long %>%
    mutate(value = 1L) %>%
    pivot_wider(
      names_from = comparison,
      values_from = value,
      values_fill = 0L
    )

  intersection_counts <- membership_wide %>%
    rowwise() %>%
    mutate(
      intersection = paste(
        comparison_levels[which(c_across(all_of(comparison_levels)) == 1)],
        collapse = " & "
      )
    ) %>%
    ungroup() %>%
    count(intersection, name = "intersection_size") %>%
    filter(intersection != "") %>%
    arrange(desc(intersection_size), intersection) %>%
    slice_head(n = top_n_intersections) %>%
    mutate(intersection_id = paste0("I", seq_len(n())))

  if (nrow(intersection_counts) == 0) {
    stop("No intersections were found for the selected comparisons.")
  }

  if (is.null(subtitle)) {
    subtitle <- paste0(
      "Top ", nrow(intersection_counts),
      " intersections based on shared gene values"
    )
  }

  matrix_df <- intersection_counts %>%
    mutate(comparison = strsplit(intersection, " & ")) %>%
    unnest(comparison) %>%
    mutate(present = 1L) %>%
    right_join(
      expand.grid(
        intersection_id = intersection_counts$intersection_id,
        comparison = comparison_levels,
        stringsAsFactors = FALSE
      ),
      by = c("intersection_id", "comparison")
    ) %>%
    left_join(
      intersection_counts %>% select(intersection_id, intersection_size),
      by = "intersection_id"
    ) %>%
    mutate(
      present = ifelse(is.na(present), 0L, present),
      intersection_id = factor(
        intersection_id,
        levels = intersection_counts$intersection_id
      ),
      comparison = factor(comparison, levels = comparison_factor_levels),
      comparison_label = factor(
        comparison_labels[as.character(comparison)],
        levels = comparison_label_levels
      )
    )

  label_df <- tibble(
    comparison = factor(comparison_factor_levels, levels = comparison_factor_levels),
    comparison_label = factor(comparison_label_levels, levels = comparison_label_levels),
    x = 1
  )

  intersection_bar <- ggplot(
    intersection_counts,
    aes(
      x = factor(intersection_id, levels = intersection_id),
      y = intersection_size
    )
  ) +
    geom_col(fill = line_color, width = 0.75) +
    geom_text(aes(label = intersection_size), vjust = -0.3, size = 3.2) +
    labs(
      title = title,
      subtitle = subtitle,
      x = NULL,
      y = "Shared genes"
    ) +
    theme_minimal(base_size = base_size) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major.x = element_blank(),
      plot.title = element_text(face = "bold"),
      plot.margin = margin(t = 8, r = 6, b = 4, l = 6)
    ) +
    expand_limits(y = max(intersection_counts$intersection_size) * 1.12)

  label_plot <- ggplot(
    label_df,
    aes(x = x, y = comparison_label, label = comparison_label)
  ) +
    geom_text(hjust = 1, size = 3.2, color = "grey30") +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(x = NULL, y = NULL) +
    theme_void(base_size = base_size) +
    theme(
      plot.margin = margin(t = 6, r = 4, b = 6, l = 8)
    )

  intersection_matrix <- ggplot(
    matrix_df,
    aes(x = intersection_id, y = comparison_label)
  ) +
    geom_point(aes(size = factor(present), color = factor(present)), stroke = 0) +
    geom_line(
      data = matrix_df %>% filter(present == 1L),
      aes(group = intersection_id),
      color = line_color,
      linewidth = 0.5
    ) +
    scale_size_manual(values = c("0" = 1.6, "1" = 4.2), guide = "none") +
    scale_color_manual(
      values = c("0" = absent_point_color, "1" = point_color),
      guide = "none"
    ) +
    labs(x = "Intersection", y = NULL) +
    theme_minimal(base_size = base_size) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 9),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.x = element_blank(),
      plot.margin = margin(t = 6, r = 6, b = 6, l = 6)
    )

  set_size_bar <- ggplot(
    set_sizes_direction,
    aes(x = set_size, y = comparison_label, fill = direction)
  ) +
    geom_col(width = 0.7) +
    scale_fill_manual(
      values = c(
        "Negative log2FC" = negative_color,
        "Positive log2FC" = positive_color
      ),
      name = NULL
    ) +
    labs(x = "Genes per comparison", y = NULL) +
    theme_minimal(base_size = base_size) +
    theme(
      panel.grid.major.y = element_blank(),
      plot.margin = margin(t = 6, r = 10, b = 6, l = 0),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = "top"
    )

  center_stack <- plot_grid(
    intersection_bar,
    intersection_matrix,
    ncol = 1,
    rel_heights = c(2.2, 1.8),
    align = "v",
    axis = "lr"
  )

  right_stack <- plot_grid(
    NULL,
    set_size_bar,
    ncol = 1,
    rel_heights = c(2.2, 1.8)
  )

  left_stack <- plot_grid(
    NULL,
    label_plot,
    ncol = 1,
    rel_heights = c(2.2, 1.8)
  )

  final_plot <- plot_grid(
    left_stack,
    center_stack,
    right_stack,
    nrow = 1,
    rel_widths = c(1.8, 4.8, 2.6),
    align = "h",
    axis = "tb"
  )

  saved_pdf <- NULL
  if (isTRUE(save_pdf)) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    }

    if (is.null(output_file)) {
      stem <- tools::file_path_sans_ext(basename(csv_file))
      output_file <- paste0(stem, "_gene_overlap_upset.pdf")
    }

    saved_pdf <- file.path(output_dir, output_file)
    ggsave(saved_pdf, final_plot, width = width, height = height, bg = "white")
  }

  result <- list(
    plot = final_plot,
    saved_pdf = saved_pdf,
    selected_comparisons = comparison_levels,
    intersection_counts = intersection_counts,
    set_sizes = set_sizes_direction
  )

  if (isTRUE(return_data)) {
    return(result)
  }

  invisible(result)
}

make_gene_overlap_summary_csv <- function(
  csv_file,
  comparisons = "all",
  comparison_col = "comparison",
  gene_col = "gene_name",
  logfc_col = "log2fc",
  padj_col = "padj",
  ast_col = "ast",
  focal_comparison = NULL,
  output_dir = getwd(),
  output_prefix = NULL,
  save_wide_csv = TRUE,
  save_long_csv = TRUE,
  return_data = FALSE
) {
  if (!file.exists(csv_file)) {
    stop("Input file not found: ", csv_file)
  }

  dat <- read.csv(csv_file, stringsAsFactors = FALSE, check.names = FALSE)
  dat <- dat[, !is.na(names(dat)) & names(dat) != "", drop = FALSE]

  required_cols <- c(comparison_col, gene_col, logfc_col)
  missing_cols <- setdiff(required_cols, names(dat))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  optional_cols <- intersect(c(padj_col, ast_col), names(dat))

  dat <- dat %>%
    transmute(
      comparison = trimws(.data[[comparison_col]]),
      gene_name = trimws(.data[[gene_col]]),
      log2fc = as.numeric(.data[[logfc_col]]),
      padj = if (padj_col %in% optional_cols) as.numeric(.data[[padj_col]]) else NA_real_,
      ast = if (ast_col %in% optional_cols) as.character(.data[[ast_col]]) else NA_character_
    ) %>%
    filter(comparison != "", gene_name != "")

  selected_comparisons <- match_requested_comparisons(
    available_comparisons = unique(dat$comparison),
    requested_comparisons = comparisons
  )

  dat <- dat %>%
    filter(comparison %in% selected_comparisons) %>%
    distinct(comparison, gene_name, .keep_all = TRUE) %>%
    mutate(
      comparison_label = pretty_comparison(comparison)
    )

  if (nrow(dat) == 0) {
    stop("No rows remain after filtering the requested comparisons.")
  }

  if (is.null(focal_comparison)) {
    focal_comparison <- selected_comparisons[1]
  } else {
    focal_comparison <- match_requested_comparisons(
      available_comparisons = selected_comparisons,
      requested_comparisons = focal_comparison
    )[1]
  }

  make_safe_name <- function(x) {
    gsub("[^A-Za-z0-9]+", "_", x)
  }

  selected_safe <- make_safe_name(selected_comparisons)
  names(selected_safe) <- selected_comparisons
  focal_safe <- selected_safe[[focal_comparison]]

  wide <- reshape(
    dat[, c("gene_name", "comparison", "log2fc", "padj", "ast")],
    idvar = "gene_name",
    timevar = "comparison",
    direction = "wide"
  )

  for (comp in selected_comparisons) {
    safe_tag <- selected_safe[[comp]]
    logfc_name <- paste0("log2fc.", comp)
    padj_name <- paste0("padj.", comp)
    ast_name <- paste0("ast.", comp)

    wide[[paste0("in_", safe_tag)]] <- !is.na(wide[[logfc_name]])
    names(wide)[names(wide) == logfc_name] <- paste0("log2fc_", safe_tag)
    names(wide)[names(wide) == padj_name] <- paste0("padj_", safe_tag)
    names(wide)[names(wide) == ast_name] <- paste0("ast_", safe_tag)
  }

  in_cols <- paste0("in_", selected_safe)
  wide[["n_comparisons_present"]] <- rowSums(wide[, in_cols, drop = FALSE])
  wide[["present_in"]] <- apply(
    wide[, in_cols, drop = FALSE],
    1,
    function(z) {
      present_names <- names(z)[z]
      if (length(present_names) == 0) {
        return("None")
      }
      paste(gsub("^in_", "", present_names), collapse = " | ")
    }
  )

  focal_in_col <- paste0("in_", focal_safe)
  other_in_cols <- setdiff(in_cols, focal_in_col)

  wide[["overlap_category"]] <- ifelse(
    wide[["n_comparisons_present"]] == length(selected_comparisons),
    paste0("Shared across all ", length(selected_comparisons)),
    ifelse(
      wide[["n_comparisons_present"]] > 1,
      paste0("Shared across ", wide[["n_comparisons_present"]]),
      ifelse(
        wide[[focal_in_col]],
        paste0(focal_safe, "-only"),
        "Unique to other selected comparison"
      )
    )
  )

  for (comp in setdiff(selected_comparisons, focal_comparison)) {
    other_safe <- selected_safe[[comp]]
    wide[[paste0("overlaps_with_", other_safe, "_from_", focal_safe)]] <-
      wide[[focal_in_col]] & wide[[paste0("in_", other_safe)]]
  }

  if (length(other_in_cols) > 1) {
    wide[[paste0("overlaps_with_all_other_selected_from_", focal_safe)]] <-
      wide[[focal_in_col]] &
      apply(wide[, other_in_cols, drop = FALSE], 1, all)
  }

  wide <- wide[order(-wide[["n_comparisons_present"]], wide[["gene_name"]]), ]
  long_out <- dat[order(dat[["comparison"]], dat[["gene_name"]]), ]

  if (is.null(output_prefix)) {
    stem <- tools::file_path_sans_ext(basename(csv_file))
    output_prefix <- paste0(stem, "_", focal_safe, "_overlap_summary")
  }

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  saved_wide_csv <- NULL
  saved_long_csv <- NULL

  if (isTRUE(save_wide_csv)) {
    saved_wide_csv <- file.path(output_dir, paste0(output_prefix, ".csv"))
    write.csv(wide, saved_wide_csv, row.names = FALSE, na = "")
  }

  if (isTRUE(save_long_csv)) {
    saved_long_csv <- file.path(output_dir, paste0(output_prefix, "_long.csv"))
    write.csv(long_out, saved_long_csv, row.names = FALSE, na = "")
  }

  result <- list(
    wide_summary = wide,
    long_summary = long_out,
    focal_comparison = focal_comparison,
    selected_comparisons = selected_comparisons,
    saved_wide_csv = saved_wide_csv,
    saved_long_csv = saved_long_csv
  )

  if (isTRUE(return_data)) {
    return(result)
  }

  invisible(result)
}


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