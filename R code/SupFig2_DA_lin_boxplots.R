suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(grid)
  library(readr)
  library(tidyr)
})

stars_padj <- function(x) {
  case_when(
    is.na(x) ~ NA_character_,
    x < 1e-4 ~ "****",
    x < 1e-3 ~ "***",
    x < 1e-2 ~ "**",
    x < 5e-2 ~ "*",
    TRUE ~ "ns"
  )
}

compute_pairwise_stats <- function(plot_df, lineage_name, control_groups, complication_groups) {
  comparison_grid <- expand.grid(
    control = control_groups,
    complication = complication_groups,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::arrange(
      match(control, control_groups),
      match(complication, complication_groups)
    )

  stats_df <- lapply(seq_len(nrow(comparison_grid)), function(i) {
    control_name <- comparison_grid$control[i]
    complication_name <- comparison_grid$complication[i]

    control_values <- plot_df %>%
      dplyr::filter(Condition_group == control_name) %>%
      dplyr::pull(Abundance)

    complication_values <- plot_df %>%
      dplyr::filter(Condition_group == complication_name) %>%
      dplyr::pull(Abundance)

    p_value <- if (length(control_values) > 0 && length(complication_values) > 0) {
      wilcox.test(complication_values, control_values, exact = FALSE)$p.value
    } else {
      NA_real_
    }

    data.frame(
      lineage = lineage_name,
      control = control_name,
      complication = complication_name,
      p_value = p_value,
      stringsAsFactors = FALSE
    )
  }) %>%
    dplyr::bind_rows()

  stats_df %>%
    dplyr::mutate(
      p_adj = p.adjust(p_value, method = "BH"),
      label = stars_padj(p_adj)
    )
}

build_lineage_plot <- function(plot_df, lineage_name, stats_df, condition_order, fill_values) {
  condition_labels <- c(
    Cont_labored = "Labored Control",
    Cont_unlabored = "Unlabored Control",
    Cont_PTB = "Preterm Control",
    PrePE = "Preterm PE"
  )
  x_positions <- stats::setNames(seq_along(condition_order), condition_order)
  data_y_max <- max(plot_df$Abundance, na.rm = TRUE)
  annotation_df <- stats_df %>%
    dplyr::filter(!is.na(label), label != "ns")

  if (nrow(annotation_df) > 0) {
    bracket_span <- max(2, data_y_max * 0.22)
    y_upper <- min(100, max(data_y_max * 1.15, data_y_max + bracket_span))
  } else {
    y_upper <- min(100, max(0.1, data_y_max * 1.15))
  }

  y_range <- max(y_upper, 0.1)

  if (nrow(annotation_df) > 0) {
    annotation_top <- y_upper * 0.88
    annotation_start <- min(
      annotation_top - max(0.3, y_range * 0.03),
      max(data_y_max + max(0.4, y_range * 0.05), y_upper * 0.70)
    )

    if (nrow(annotation_df) == 1) {
      y_positions <- annotation_start
    } else {
      y_positions <- seq(annotation_start, annotation_top, length.out = nrow(annotation_df))
    }

    step_height <- if (nrow(annotation_df) > 1) min(diff(y_positions)) else max(0.2, y_range * 0.03)
    bracket_height <- min(max(0.15, y_range * 0.02), step_height * 0.35)

    annotation_df <- annotation_df %>%
      dplyr::mutate(
        x_start = x_positions[control],
        x_end = x_positions[complication],
        y_position = y_positions
      )
  } else {
    bracket_height <- 1
  }

  base_plot <- ggplot(plot_df, aes(x = Condition_group, y = Abundance, fill = Group_type)) +
    geom_violin(
      width = 0.85,
      trim = FALSE,
      alpha = 0.5,
      color = "black",
      linewidth = 0.3,
      draw_quantiles = c(0.25, 0.5, 0.75)
    ) +
    geom_jitter(
      aes(color = Group_type),
      width = 0.16,
      size = 1.5,
      alpha = 0.95
    ) +
    scale_fill_manual(values = fill_values, guide = "none") +
    scale_color_manual(values = fill_values, guide = "none") +
    scale_x_discrete(drop = FALSE, labels = function(x) dplyr::coalesce(condition_labels[x], x)) +
    scale_y_continuous(
      limits = c(0, y_upper),
      expand = expansion(mult = c(0, 0.04))
    ) +
    labs(
      title = paste0("Distribution of ", lineage_name, " Abundance\nby Condition"),
      x = NULL,
      y = paste0("Abundance of ", lineage_name, " (%)")
    ) +
    theme_classic(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      plot.title = element_text(hjust = 0.5, size = 10.5, face = "bold"),
      plot.margin = margin(8, 8, 8, 8)
    )

  if (nrow(annotation_df) > 0) {
    base_plot <- base_plot +
      geom_segment(
        data = annotation_df,
        aes(x = x_start, xend = x_end, y = y_position, yend = y_position),
        inherit.aes = FALSE,
        linewidth = 0.35
      ) +
      geom_segment(
        data = annotation_df,
        aes(x = x_start, xend = x_start, y = y_position - bracket_height, yend = y_position),
        inherit.aes = FALSE,
        linewidth = 0.35
      ) +
      geom_segment(
        data = annotation_df,
        aes(x = x_end, xend = x_end, y = y_position - bracket_height, yend = y_position),
        inherit.aes = FALSE,
        linewidth = 0.35
      ) +
      geom_text(
        data = annotation_df,
        aes(x = (x_start + x_end) / 2, y = y_position + bracket_height * 0.15, label = label),
        inherit.aes = FALSE,
        size = 3
      )
  }

  base_plot
}

sample_condition_path <- "/projects/Pregnancy_placenta/Results/13_DE/sample_condition_list.csv"
sample_values_path <- "/projects/Pregnancy_placenta/Results/20_Figures/DA_lineage_level/AllCells_v3_frequencies_per_pt.csv"
plot_output_dir <- "/projects/Pregnancy_placenta/Results/20_Figures/DA_lineage_level/SupFig_DA_plots"
dir.create(plot_output_dir, recursive = TRUE, showWarnings = FALSE)

control_groups <- c("Cont_labored", "Cont_unlabored", "Cont_PTB")
complication_groups <- c("PTB", "DM1", "FGR", "PrePE", "TermPE")
condition_order <- c(control_groups, complication_groups)
fill_values <- c(Control = "#868686FF", Complication = "#CD534CFF")

samp_condition_file <- readr::read_csv(sample_condition_path, show_col_types = FALSE) %>%
  dplyr::select(Sample, Condition_group) %>%
  dplyr::distinct()

samp_values_df <- readr::read_csv(sample_values_path, show_col_types = FALSE)
lineages_to_plot <- setdiff(colnames(samp_values_df), "Sample")

multi_condition_samples <- samp_condition_file %>%
  dplyr::count(Sample, name = "n_conditions") %>%
  dplyr::filter(n_conditions > 1)

plotting_df <- samp_values_df %>%
  tidyr::pivot_longer(cols = -Sample, names_to = "Lineage", values_to = "Abundance") %>%
  dplyr::inner_join(
    samp_condition_file,
    by = "Sample",
    relationship = "many-to-many"
  ) %>%
  dplyr::filter(Condition_group %in% condition_order) %>%
  dplyr::mutate(
    Condition_group = factor(Condition_group, levels = condition_order, ordered = TRUE),
    Group_type = ifelse(Condition_group %in% control_groups, "Control", "Complication"),
    Group_type = factor(Group_type, levels = c("Control", "Complication"))
  )

all_stats <- lapply(lineages_to_plot, function(lineage_name) {
  compute_pairwise_stats(
    plot_df = dplyr::filter(plotting_df, Lineage == lineage_name),
    lineage_name = lineage_name,
    control_groups = control_groups,
    complication_groups = complication_groups
  )
}) %>%
  dplyr::bind_rows()

all_plots <- lapply(lineages_to_plot, function(lineage_name) {
  build_lineage_plot(
    plot_df = dplyr::filter(plotting_df, Lineage == lineage_name),
    lineage_name = lineage_name,
    stats_df = dplyr::filter(all_stats, lineage == lineage_name),
    condition_order = condition_order,
    fill_values = fill_values
  )
})

names(all_plots) <- lineages_to_plot

combined_pdf_path <- file.path(plot_output_dir, "SupFig2_DA_all_lineages_V1.pdf")
pdf(combined_pdf_path, width = 8.5, height = 11, onefile = TRUE)
grid::grid.newpage()
grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow = 5, ncol = 3)))

combined_plots <- all_plots[names(all_plots) != "Other"]

combined_plots <- lapply(combined_plots, function(curr_plot) {
  curr_plot +
    theme(
      plot.title = element_text(hjust = 0.5, size = 6.5, lineheight = 0.9, face = "bold"),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 6.5),
      axis.text.x = element_text(size = 5.5, angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y = element_text(size = 5.5),
      plot.margin = margin(2, 2, 2, 2)
    )
})

plot_index <- 1
for (row_idx in seq_len(5)) {
  for (col_idx in seq_len(3)) {
    if (plot_index <= length(combined_plots)) {
      print(
        combined_plots[[plot_index]],
        vp = grid::viewport(layout.pos.row = row_idx, layout.pos.col = col_idx)
      )
    }
    plot_index <- plot_index + 1
  }
}

dev.off()

readr::write_csv(
  all_stats,
  file.path(plot_output_dir, "SupFig2_DA_lineage_boxplot_stats.csv")
)

invisible(all_plots)
