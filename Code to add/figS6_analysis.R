# Purpose of Script ----
# Figure S6: compare pathway results using unlabored and labored controls.
# Run from the FigS6_code folder, either all at once or section by section.

# Preparing Packages and Functions ----
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readxl)
  library(systemfonts)
})

source("figS6_helper.R")

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
                 "inputs/label_positions.tsv", "figS6_analysis.R", "figS6_helper.R")
write_tsv(data.frame(File = input_files, MD5 = unname(tools::md5sum(input_files))),
          file.path(output_dir, "Input_checksums.tsv"))
write_tsv(data.frame(Adjusted_p_cutoff = padj_cutoff, Absolute_Z_cutoff = zscore_cutoff,
                     Pathways = nrow(selected_pathways), Populations = nrow(populations),
                     Complications = length(complications), Shared_results = nrow(shared_results)),
          file.path(output_dir, "Analysis_settings.tsv"))
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))

message("Figure S6 PDFs and tables saved in: ", normalizePath(output_dir))
