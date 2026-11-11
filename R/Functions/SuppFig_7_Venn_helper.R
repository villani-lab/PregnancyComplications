# Helper Functions for Supp Figure 7 (Venn diagram and shared-pathway scatter) ----
# Documentation: reusable input checks and figure formatting for R/SuppFig_7_Venn_scatter.R, which
# loads the inputs, applies |Z-score| > 2 and adjusted P < 0.25 to the labored/unlabored control
# comparisons and saves the Venn diagram, Z-score scatter and source tables (see that script).
# (Labelled "Figure S6" in the original code; Supp Fig 7 in the revised paper.)
# Sourced by SuppFig_7_Venn_scatter.R; this file does not run the analysis itself.
# Sections: input names, TEST3 source checks, table exports, figure formatting.

# Input Names ----
clean_ipa_names <- function(x) {
  x <- enc2utf8(trimws(as.character(x)))
  gsub("\u2212|\u2013|\u2011", "-", x)
}

cell_display_names <- function(x) {
  # Translate saved object names to the names used in the manuscript.
  aliases <- c(
    HBC_2_CXCL2_CXCL3 = "HBC_2", HBC_4_IL1B_MMP19 = "HBC_4",
    HBC_7_lnRNA = "HBC_7", HBC_8_neuro = "HBC_8",
    HBC_9_CXCL10_IFIT44 = "HBC_9", HBC_10_cycling = "HBC_10",
    mixed_HBC_13_S100A8 = "HBC_13",
    VCTp_1 = "VCT_cycling_1", VCTp_2 = "VCT_cycling_2",
    Endo_f1 = "fEndo_1", Endo_f2 = "fEndo_2", Endo_f3 = "fEndo_3",
    Endo_m1 = "mEndo_1", fFB1 = "fFB_1", fPV2 = "fPV_2",
    fPV3 = "fPV_3", dS1 = "dS_1", dS3 = "dS_3", mMac_upper = "mMac_4_6_7"
  )

  labels <- sub("_L[12]$", "", x)
  matched <- labels %in% names(aliases)
  labels[matched] <- unname(aliases[labels[matched]])
  labels
}

prepare_ipa_results <- function(canonical_pathways, disease_functions) {
  # Keep the source identifiers needed to locate the original Excel entries.
  columns <- c("Canonical_Pathways", "Cluster", "Clust_Comparison", "Comparison",
               "Condition", "Control", "CompNum", "Zscore", "log10_padj")
  canonical <- as.data.frame(canonical_pathways)[, columns]
  diseases <- as.data.frame(disease_functions)[, columns]
  canonical$Database <- "canonical"
  diseases$Database <- "disease"
  results <- rbind(canonical, diseases)

  name_columns <- c("Canonical_Pathways", "Cluster", "Clust_Comparison",
                    "Comparison", "Condition", "Control")
  results[name_columns] <- lapply(results[name_columns], clean_ipa_names)
  results$Pathway <- results$Canonical_Pathways
  results$Pathway_id <- paste(results$Database, results$Pathway, sep = "::")
  results$Cell_type <- results$Cluster
  results$Cell_label <- cell_display_names(results$Cluster)
  results$Source_level <- ifelse(endsWith(results$Cluster, "_L1"), "parent", "sublineage")

  expected_comparison <- paste(results$Condition, "vs", results$Control, sep = "_")
  expected_header <- paste(sub("_L[12]$", "", results$Cluster),
                           expected_comparison, sep = "-")
  stopifnot(all(grepl("_L[12]$", results$Cluster)),
            all(results$Comparison == expected_comparison),
            all(results$Clust_Comparison == expected_header))
  results
}

# Check the Original TEST3 Statistics ----
check_ipa_statistics <- function(results, project_dir) {
  # The RDS objects define which comparisons belong in this analysis.
  # Check each selected Z and BH value against the corresponding TEST3 matrix.
  # Keep genuine missing values; never fill them from another IPA report.
  result_keys <- results[c("Pathway_id", "Cell_type", "Condition", "Control")]
  if (anyDuplicated(result_keys)) {
    stop("Duplicate pathway/cell/complication/control results in the saved inputs.")
  }

  files <- character()
  checks <- list()
  results$Zscore_file <- results$log10_padj_file <- NA_character_
  results$Zscore_Excel_row <- results$log10_padj_Excel_row <- NA_integer_

  batch_groups <- split(seq_len(nrow(results)),
                        paste(results$Source_level, results$CompNum, results$Database))

  for (rows in batch_groups) {
    batch <- results[rows[1], ]

    if (batch$Source_level == "parent") {
      stopifnot(batch$CompNum %in% 1:4)
      batch_folder <- paste0("Parent_lin_test3_CompNum_", batch$CompNum)
      if (batch$CompNum == 4) batch_folder <- "Parent_lin_test3b_CompNum_1"
      folder <- file.path(project_dir, "IPA_results/Parent_lin", batch_folder)
    } else {
      stopifnot(batch$CompNum %in% 1:5)
      batch_folder <- paste0("Sublin_test3_CompNum_", batch$CompNum)
      folder <- file.path(project_dir, "IPA_results/sub_lin", batch_folder)
    }

    database_name <- if (batch$Database == "canonical") "canonical_pathways" else "diseases"
    matrix_suffixes <- c(Zscore = "_zscore.xls", log10_padj = "_BHpvalue.xls")

    for (metric in names(matrix_suffixes)) {
      suffix <- paste0(database_name, matrix_suffixes[[metric]])
      candidates <- list.files(folder, full.names = TRUE)
      source_file <- candidates[endsWith(candidates, suffix)]
      if (length(source_file) != 1) {
        stop("Expected one file ending in ", suffix, " in ", folder)
      }

      raw_matrix <- as.data.frame(readxl::read_xls(
        source_file, skip = 1, col_types = "text", .name_repair = "minimal"
      ))
      pathway_names <- clean_ipa_names(raw_matrix[[1]])
      comparison_names <- clean_ipa_names(names(raw_matrix)[-1])
      stopifnot(!anyNA(pathway_names), !anyDuplicated(pathway_names),
                !anyDuplicated(comparison_names))

      pathway_row <- match(results$Pathway[rows], pathway_names)
      comparison_column <- match(results$Clust_Comparison[rows], comparison_names)
      if (anyNA(pathway_row) || anyNA(comparison_column)) {
        stop("A saved pathway/comparison is missing from the original matrix: ", source_file)
      }

      raw_text <- as.matrix(raw_matrix[-1])[cbind(pathway_row, comparison_column)]
      raw_values <- suppressWarnings(as.numeric(raw_text))
      unexpected <- is.na(raw_values) & !is.na(raw_text) &
        !raw_text %in% c("", "NA", "N/A", "NaN", "-")
      if (any(unexpected)) stop("Unexpected nonnumeric value in ", source_file)
      if (metric == "log10_padj" && any(raw_values < 0, na.rm = TRUE)) {
        stop("Negative -log10(adjusted P) in ", source_file)
      }

      saved_values <- results[[metric]][rows]
      matches <- is.na(saved_values) & is.na(raw_values)
      numeric_pairs <- !is.na(saved_values) & !is.na(raw_values)
      matches[numeric_pairs] <- saved_values[numeric_pairs] == raw_values[numeric_pairs] |
        (is.finite(saved_values[numeric_pairs]) & is.finite(raw_values[numeric_pairs]) &
           abs(saved_values[numeric_pairs] - raw_values[numeric_pairs]) < 1e-10)
      mismatches <- sum(!matches)
      if (mismatches > 0) {
        stop(mismatches, " saved/raw statistic mismatches in ", source_file)
      }

      relative_file <- sub(paste0(normalizePath(project_dir), "/"), "",
                            normalizePath(source_file), fixed = TRUE)
      results[[metric]][rows] <- raw_values
      results[[paste0(metric, "_file")]][rows] <- relative_file
      results[[paste0(metric, "_Excel_row")]][rows] <- pathway_row + 2
      files <- union(files, source_file)
      checks[[length(checks) + 1]] <- data.frame(
        Source_level = batch$Source_level, Batch = batch$CompNum,
        Database = batch$Database, Metric = metric,
        Checked_values = length(rows), Mismatches = mismatches,
        Source_file = relative_file
      )
    }
  }

  list(data = results, checks = do.call(rbind, checks), files = files)
}

# Save Tables ----
write_tsv <- function(data, filename) {
  # Preserve Greek letters when R starts with an ASCII-only locale.
  previous_locale <- Sys.getlocale("LC_CTYPE")
  if (!l10n_info()[["UTF-8"]]) {
    Sys.setlocale("LC_CTYPE", "en_US.UTF-8")
    on.exit(Sys.setlocale("LC_CTYPE", previous_locale))
  }
  write.table(data, filename, sep = "\t", quote = TRUE, qmethod = "double",
               row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
}

save_complication_tables <- function(results, pathways, populations, complications, output_dir) {
  table_dir <- file.path(output_dir, "by_complication")
  dir.create(table_dir, showWarnings = FALSE)
  pathway_summaries <- list()

  for (condition in names(complications)) {
    condition_results <- results[results$Condition == condition, ]
    significance_matrix <- data.frame(Pathway = pathways$Pathway, check.names = FALSE)

    # Each population gets two columns: unlabored X/- and labored X/-.
    for (i in seq_len(nrow(populations))) {
      cell_results <- condition_results[condition_results$Cell_type == populations$Cell_type[i], ]
      pathway_rows <- match(pathways$Pathway_id, cell_results$Pathway_id)

      for (control in c("Unlabored", "Labored")) {
        column_name <- paste(populations$Cell_label[i], control, sep = " | ")
        significance_matrix[[column_name]] <- cell_results[[paste0(control, "_flag")]][pathway_rows]
      }
    }

    counts <- as.data.frame(pathways)
    for (category in c("Unlabored only", "Both", "Labored only")) {
      counts[[gsub(" ", "_", category)]] <- vapply(pathways$Pathway_id, function(pathway) {
        sum(condition_results$Pathway_id == pathway & condition_results$Category == category)
      }, integer(1))
    }
    counts$Total <- counts$Unlabored_only + counts$Both + counts$Labored_only
    pathway_summaries[[condition]] <- counts

    write_tsv(significance_matrix, file.path(table_dir, paste0(condition, "_significance_matrix.tsv")))
    write_tsv(counts, file.path(table_dir, paste0(condition, "_pathway_counts.tsv")))
  }

  total_counts <- pathway_summaries[[1]]
  for (column in c("Unlabored_only", "Both", "Labored_only", "Total")) {
    total_counts[[column]] <- Reduce(`+`, lapply(pathway_summaries, function(x) x[[column]]))
  }
  write_tsv(total_counts, file.path(output_dir, "Consolidated_pathway_counts.tsv"))
}

# Check Figure Inputs ----
check_plot_inputs <- function(results, labels, colors, family) {
  if (!capabilities("cairo")) stop("Cairo PDF support is needed to save the figures.")
  if (!family %in% systemfonts::system_fonts()$family) stop("Install the font: ", family)
  if (anyNA(labels$Label_x) || anyNA(labels$Label_y)) {
    stop("Add positions for the new shared results to inputs/label_positions.tsv.")
  }
  if (!all(labels$Cell_label %in% names(colors))) {
    stop("Add colors for the new shared populations to inputs/cell_colors.tsv.")
  }
  if (any(abs(c(results$Unlabored_Zscore, results$Labored_Zscore)) > 5.5)) {
    stop("A Z-score exceeds the figure axes. Update the axis limits and label layout.")
  }
}

# Scatter Labels ----
make_label_lines <- function(labels, family) {
  # End each leader line beside its text, leaving space around the point.
  # The publication panel is 340 points wide and spans 11 Z-score units.
  points_per_z <- 340 / 11
  labels$Label <- sub("Interferon alpha/beta signaling",
                      "Interferon alpha/beta\nsignaling", labels$Pathway, fixed = TRUE)
  labels$Start_x <- labels$Start_y <- labels$End_x <- labels$End_y <- NA_real_

  for (i in seq_len(nrow(labels))) {
    text_lines <- strsplit(labels$Label[i], "\n", fixed = TRUE)[[1]]
    text_width <- max(systemfonts::string_width(text_lines, family = family, size = 8))
    text_height <- 8 + (length(text_lines) - 1) * 9.6
    dx <- labels$Unlabored_Zscore[i] - labels$Label_x[i]
    dy <- labels$Labored_Zscore[i] - labels$Label_y[i]
    horizontal_edge <- if (dx == 0) Inf else (text_width / 2 + 2) / points_per_z / abs(dx)
    vertical_edge <- if (dy == 0) Inf else (text_height / 2 + 2) / points_per_z / abs(dy)
    fraction <- min(horizontal_edge, vertical_edge)

    end_x <- labels$Label_x[i] + fraction * dx
    end_y <- labels$Label_y[i] + fraction * dy
    distance <- sqrt((end_x - labels$Unlabored_Zscore[i])^2 +
                       (end_y - labels$Labored_Zscore[i])^2)
    if (!is.finite(distance) || distance <= 3.5 / points_per_z) {
      stop("Move the label further from its point: ", labels$Pathway[i])
    }
    point_gap <- 3.5 / points_per_z / distance
    labels$Start_x[i] <- labels$Unlabored_Zscore[i] +
      point_gap * (end_x - labels$Unlabored_Zscore[i])
    labels$Start_y[i] <- labels$Labored_Zscore[i] +
      point_gap * (end_y - labels$Labored_Zscore[i])
    labels$End_x[i] <- end_x
    labels$End_y[i] <- end_y
  }
  labels
}

# Venn Formatting ----
plot_pathway_venn <- function(counts, colors, title, family) {
  # Coordinates below reproduce the approved equal-circle layout (280 x 246 pt).
  # Areas are schematic; the numbers come from the threshold counts.
  text_label <- function(label, x, y, size = 8.5, bold = FALSE) {
    grid::textGrob(label, x = x / 280, y = y / 246,
                   gp = grid::gpar(fontfamily = family, fontsize = size,
                                   fontface = if (bold) "bold" else "plain"))
  }
  grid::grobTree(
    grid::circleGrob(x = 97 / 280, y = 105 / 246, r = grid::unit(66, "bigpts"),
                     gp = grid::gpar(fill = grDevices::adjustcolor(colors[1], 0.5), col = NA)),
    grid::circleGrob(x = 183 / 280, y = 105 / 246, r = grid::unit(66, "bigpts"),
                     gp = grid::gpar(fill = grDevices::adjustcolor(colors[2], 0.5), col = NA)),
    text_label(title, 140, 220, size = 10.5, bold = TRUE),
    text_label("Unlabored control", 86, 184),
    text_label("Labored control", 194, 184),
    text_label(counts$Unlabored_only, 71, 105, size = 15, bold = TRUE),
    text_label(counts$Both, 140, 105, size = 15, bold = TRUE),
    text_label(counts$Labored_only, 209, 105, size = 15, bold = TRUE),
    text_label(paste0("Total: ", counts$Total_unlabored), 86, 28),
    text_label(paste0("Total: ", counts$Total_labored), 194, 28)
  )
}

# Center the Scatter on Its PDF Page ----
save_scatter_pdf <- function(plot, filename, family) {
  # Keep the existing 440 x 476 pt page and 340 pt square plotting panel.
  # Draw ggplot's panel and axes at fixed positions so unequal axis-label widths
  # do not shift the panel. The title and bottom legend are centered separately.
  grDevices::cairo_pdf(filename, width = 440 / 72, height = 476 / 72, family = family)
  on.exit(grDevices::dev.off())
  plot_grobs <- ggplot2::ggplotGrob(plot)
  get_part <- function(name) plot_grobs$grobs[[match(name, plot_grobs$layout$name)]]

  draw_part <- function(part, x, y, width, height) {
    grid::pushViewport(grid::viewport(x = x / 440, y = y / 476,
                                      width = grid::unit(width, "bigpts"),
                                      height = grid::unit(height, "bigpts")))
    grid::grid.draw(part)
    grid::popViewport()
  }

  grid::grid.newpage()
  draw_part(get_part("panel"), 220, 256, 340, 340)

  bottom_axis <- get_part("axis-b")
  axis_height <- grid::convertHeight(grid::grobHeight(bottom_axis), "bigpts", valueOnly = TRUE)
  draw_part(bottom_axis, 220, 86 - axis_height / 2, 340, axis_height)

  left_axis <- get_part("axis-l")
  axis_width <- grid::convertWidth(grid::grobWidth(left_axis), "bigpts", valueOnly = TRUE)
  draw_part(left_axis, 50 - axis_width / 2, 256, axis_width, 340)

  grid::grid.text(plot$labels$title, x = 220 / 440, y = 461 / 476,
                  gp = grid::gpar(fontfamily = family, fontsize = 10.5, fontface = "bold"))
  grid::grid.text(plot$labels$x, x = 220 / 440, y = 62 / 476,
                  gp = grid::gpar(fontfamily = family, fontsize = 8.5))
  grid::grid.text(plot$labels$y, x = 16 / 440, y = 256 / 476, rot = 90,
                  gp = grid::gpar(fontfamily = family, fontsize = 8.5))

  legend <- get_part("guide-box-bottom")
  legend_width <- grid::convertWidth(grid::grobWidth(legend), "bigpts", valueOnly = TRUE)
  legend_height <- grid::convertHeight(grid::grobHeight(legend), "bigpts", valueOnly = TRUE)
  draw_part(legend, 220, 28, legend_width, legend_height)
  invisible(filename)
}
