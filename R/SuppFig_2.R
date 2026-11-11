# Supp Fig 2) CellTypist comparison
# ------ dir org
setwd("/projects/Pregnancy_placenta/Results/22_Final_figures/SuppFigure_2_celltypist/")

projdir="/projects/home/rpragya/pregnancy_plac/label_comparisons"
inputdir=file.path(projdir,"data")
savedir=file.path("./results")
preddir=file.path("/projects/home/rpragya/pregnancy_plac/label_comparisons/results/predictions/")

# ------ functions
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(gtools)
})

celltypist_dotplot <- function(predicted_labels_csv,
                               probability_matrix_csv,
                               reference_csv,
                               lineage_map,                 # named vector/list: reference label -> lineage
                               lineage_colors,              # named vector: lineage -> colour
                               lineage_order      = names(lineage_colors),
                               annotation_bars    = TRUE,
                               reference_col      = "Cell_subset",
                               use_as_prediction  = c("predicted_labels", "majority_voting"),
                               prediction_order   = NULL,   # x axis, left to right (subset ok)
                               reference_order    = NULL,   # y axis, top to bottom (subset ok)
                               filter_prediction  = 0,      # drop predicted labels whose max fraction < this
                               dot_fill_color = 'RdBu',
                               title              = "CellTypist label transfer",
                               move_last = NULL, # labels to move to bottom of a lineage block on y axis
                               base_size          = 12,
                               max_size           = 6,
                               out_prefix="") {
  ## Reproduce celltypist.dotplot() in R from the CSVs written by
  ## predictions.to_table(folder=..., prefix=...)
  ##
  ## Inputs
  ##   predicted_labels.csv   : rows = cells; columns predicted_labels, majority_voting, conf_score, ...
  ##   probability_matrix.csv : rows = cells; one column per model class (probabilities)
  ##   reference labels       : csv with the query cells' own annotation
  ##                            (export from python: adata.obs[['annotation']].to_csv('reference.csv'))
  ##   lineage_map            : named vector / list, REFERENCE label -> lineage
  ##
  ## Dot size   = fraction of cells in each reference group assigned a given predicted label
  ## Dot colour = mean probability of the assigned label, averaged over the cells of the
  ##              reference group that received that label (same as celltypist.dotplot)
  ##
  ## Axes: predicted labels on x, reference labels on y.
  ##   y (reference): lineages in `lineage_order`, subsets by cell count within a lineage,
  ##                  unless `reference_order` is given (top to bottom).
  ##   x (predicted): each predicted label is "owned" by the reference subset in which it is
  ##                  most frequent; it inherits that subset's lineage (colour) and is placed
  ##                  under it, ordered by fraction, unless `prediction_order` is given.
  ## Lines separate lineages; axis text is coloured by lineage.

  use_as_prediction <- match.arg(use_as_prediction)
  lineage_map <- unlist(lineage_map)

  ## ---- read ---------------------------------------------------------------
  lab  <- fread(predicted_labels_csv)
  prob <- fread(probability_matrix_csv)
  ref  <- fread(reference_csv)
  setnames(lab,  1, "cell"); setnames(prob, 1, "cell"); setnames(ref, 1, "cell")
  setnames(ref, reference_col, "reference")

  cells <- Reduce(intersect, list(lab$cell, prob$cell, ref$cell))
  lab  <- lab[match(cells, cell)]
  prob <- prob[match(cells, cell)]
  ref  <- ref[match(cells, cell)]

  classes   <- setdiff(names(prob), "cell")
  reference <- as.character(ref$reference)
  predicted <- as.character(lab[[use_as_prediction]])

  ## ---- fraction of cells (dot size) ----------------------------------------
  ct <- table(factor(reference), factor(predicted, levels = classes))
  frac <- as.data.table(as.data.frame(sweep(ct, 1, rowSums(ct), "/")))
  setnames(frac, c("reference", "prediction", "fraction"))

  ## ---- mean probability of the assigned label (dot colour) -----------------
  pmat  <- as.matrix(prob[, ..classes])
  hit   <- match(predicted, classes)
  score <- ifelse(is.na(hit), apply(pmat, 1, max), pmat[cbind(seq_len(nrow(pmat)), hit)])
  meanp <- data.table(reference = reference, prediction = predicted, score = score)[
    , .(mean_prob = mean(score)), by = .(reference, prediction)]
  df <- merge(frac, meanp, by = c("reference", "prediction"), all.x = TRUE)
  df[is.na(mean_prob), mean_prob := 0]
  df[, `:=`(reference = as.character(reference), prediction = as.character(prediction))]

  ## ---- reference (y) order --------------------------------------------------
  ## lineage_map covers the REFERENCE labels only
  missing <- setdiff(unique(df$reference), names(lineage_map))
  if (length(missing)) stop("Reference labels missing from lineage_map: ", paste(missing, collapse = ", "))

  if (is.null(reference_order)) {
    ## lineages in lineage_order; within a lineage, subsets in natural alphanumeric order
    n_ref <- as.data.table(table(reference))[, .(reference = as.character(reference), n = N)]
    n_ref[, lin := factor(lineage_map[reference], levels = lineage_order)]
    n_ref <- n_ref[order(lin, reference %in% move_last, gtools::mixedorder(reference))]
    reference_order <- n_ref$reference
  } else {
    reference_order <- intersect(reference_order, unique(df$reference))
  }
  df <- df[reference %in% reference_order]

  ## ---- predicted (x) lineage + order ------------------------------------------
  ## each predicted label is "owned" by the reference subset in which it is most
  ## frequent; it inherits that subset's lineage and is placed under it
  keep  <- df[, .(mx = max(fraction)), by = prediction][mx >= filter_prediction, prediction]
  df    <- df[prediction %in% keep]
  owner <- df[, .SD[which.max(fraction)], by = prediction][, .(prediction, owner = reference, fraction)]
  owner[, pos := match(owner, reference_order)]
  setorder(owner, pos, -fraction)
  pred_lineage <- setNames(lineage_map[owner$owner], owner$prediction)

  if (is.null(prediction_order)) {
    prediction_order <- owner$prediction
  } else {
    prediction_order <- intersect(prediction_order, owner$prediction)
  }
  df <- df[prediction %in% prediction_order]

  x_levels <- prediction_order
  y_levels <- rev(reference_order)                        # first label on top
  df[, prediction := factor(prediction, levels = x_levels)]
  df[, reference  := factor(reference,  levels = y_levels)]
  write.table(df, file=paste0(out_prefix,"_cell_typist_data.csv"), sep = ",", row.names = FALSE)

  ## lineage boundaries (k + 0.5) and axis text colours
  breaks_of <- function(l) which(l[-1] != l[-length(l)]) + 0.5
  x_lin <- unname(pred_lineage[x_levels])
  y_lin <- unname(lineage_map[y_levels])
  x_breaks   <- breaks_of(x_lin)
  y_breaks   <- breaks_of(y_lin)
  x_bar <- data.table(x = seq_along(x_levels), y = 0.2, fill = unname(lineage_colors[x_lin]))
  y_bar <- data.table(x = 0.2, y = seq_along(y_levels), fill = unname(lineage_colors[y_lin]))

  ## invisible points to get a lineage legend
  lin_present <- intersect(lineage_order, unique(c(x_lin, y_lin)))
  lin_df <- data.table(lineage = factor(lin_present, levels = lineage_order))

  ## ---- plot -----------------------------------------------------------------
  ggplot(df, aes(x = prediction, y = reference)) +
    annotate("segment", x = x_breaks, xend = x_breaks, y = 0.5, yend = length(y_levels) + 0.5,
            colour = 'grey40', linewidth = 0.4) +
    annotate("segment", y = y_breaks, yend = y_breaks, x = 0.5, xend = length(x_levels) + 0.5,
          colour = 'grey40', linewidth = 0.4) +
    geom_point(aes(size = fraction * 100, fill = mean_prob),
               shape = 21, colour = "black", stroke = 0.3) +
    { if (annotation_bars) list(
        geom_tile(data = x_bar, aes(x = x, y = y), fill = x_bar$fill, width = 1,   height = 0.4, inherit.aes = FALSE),
        geom_tile(data = y_bar, aes(x = x, y = y), fill = y_bar$fill, width = 0.4, height = 1,   inherit.aes = FALSE)) } +
    scale_x_discrete(expand = expansion(add = if (annotation_bars) c(0.9, 0.6) else 0.6)) +
    scale_y_discrete(expand = expansion(add = if (annotation_bars) c(0.9, 0.6) else 0.6)) +
    geom_point(data = lin_df, aes(x = -Inf, y = -Inf, colour = lineage),
               inherit.aes = FALSE, alpha = 0, show.legend = TRUE) +
    scale_colour_manual(values = lineage_colors, name = "Lineage",
                        guide = guide_legend(override.aes = list(alpha = 1, size = 4, shape = 15))) +
    scale_fill_distiller(palette = dot_fill_color, direction = -1, limits = c(0, 1),
                         name = "Mean probability") +
    scale_size_area(max_size = max_size, limits = c(0, 100),
                    breaks = c(20, 40, 60, 80, 100), name = "Fraction of cells (%)") +
    labs(title = title, x = NULL, y = NULL) +
    theme_bw(base_size = base_size) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = base_size, colour = "black"),
          axis.text.y = element_text(size = base_size, colour = "black"),
          plot.title   = element_text(size = base_size * 1.3, hjust = 0.5),
          legend.title = element_text(size = base_size),
          legend.text  = element_text(size = base_size * 0.9),
          panel.grid   = element_blank())
}


## ---- plot ---------------------------------------------------------------
ref_df <- read.csv(file.path(preddir, "normand_pred_using_vt/reference.csv"))
lineage_map <- unique(ref_df[, c("Cell_subset", "Lineage")])
cell_to_lineage <- setNames(lineage_map$Lineage, lineage_map$Cell_subset)


palette<- c('Trophoblasts'= '#6aa96d',
"EVT"= "#b2df8a",
"fFibroblasts"="#6cb1db",
"fPerivascular"="#6388c4",
"dStromal"= "#28315d",
"Endothelial"="#967ed9",
"HBC"= "#6f102d",
"mMac"= "#d64e3b",
"Monocytes"= "#efaabb",
"DC"= "#e3c087",
"T"= "#e59154",
"NK"= "#edd188",
"B_Plasma"= "#e80d7a",
"HSC"= "#808080",
"Other"="#d3d3d3")

# vt
vt_x_order<- c(
  "VCT",
  "SCT",
  "EVT",
  "fFB1",
  "fFB2",
  "dS1",
  "dS2",
  "dS3",
  "dP1",
  "dP2",
  "Endo (f)",
  "Endo (m)",
  "Endo L",
  "HB",
  "dM2",
  "dM3",
  "MO",
  "dM1",
  "DC1",
  "DC2",
  "Plasma",
  "Tcells",
  "ILC3",
  "dNK3",
  "NK CD16+",
  "dNK2",
  "dNK p",
  "NK CD16-",
  "Granulocytes"
)


p1 <- celltypist_dotplot(
  predicted_labels_csv   = file.path(preddir, "normand_pred_using_vt/predicted_labels.csv"),
  probability_matrix_csv = file.path(preddir, "normand_pred_using_vt/probability_matrix.csv"),
  reference_csv          = file.path(preddir, "normand_pred_using_vt/reference.csv"),   # columns: cell, annotation
  lineage_map            = cell_to_lineage,
  lineage_colors         = palette,
  use_as_prediction      = "predicted_labels",
  prediction_order   = vt_x_order,
  move_last = 'SCT',
  filter_prediction      = 0.1,
  out_prefix="vt")
ggsave(file.path(savedir, "figures_final/vt_dotplot.pdf"), p1, width = 10, height = 19, useDingbats = FALSE)

# shook
shook_x_axis = c(
"HBCs_0",
"HBCs_1",
"HBCs_2",
"HBCs_3",
"HBCs_4",
"HBCs_5",
"HBCs_6",
"HBCs_7",
"Monocytes",
"PAMs"
)
p2 <- celltypist_dotplot(
  predicted_labels_csv   = file.path(preddir, "normand_pred_using_shook/predicted_labels.csv"),
  probability_matrix_csv = file.path(preddir, "normand_pred_using_shook/probability_matrix.csv"),
  reference_csv          = file.path(preddir, "normand_pred_using_shook/reference.csv"),   # columns: cell, annotation
  lineage_map            = cell_to_lineage,
  lineage_colors         = palette,
  annotation_bars    = FALSE,
  dot_fill_color = 'RdBu',
  prediction_order   = shook_x_axis,
  use_as_prediction      = "predicted_labels",
  filter_prediction      = 0.1,
  out_prefix="shook")
ggsave(file.path(savedir, "figures_final/shook_dotplot_rdbu.pdf"), p2, width = 6, height = 7, useDingbats = FALSE)





## Map cells - mutual best match
# two cell subsets x from Data_x and y from Data_y will be considered "matching if each has the highest probability to map to the other out of all cell subsets
# A mutual best match is a cell (x, y) that is simultaneously the row maximum and the column maximum.

# Shook data
df = read.table("shook_cell_typist_data.csv", sep=",", header = TRUE)

#df_table = reshape2::dcast(df[,c("reference", "prediction", "mean_prob")], reference ~ prediction, value.var = "mean_prob")
df_table = reshape2::dcast(df[,c("reference", "prediction", "fraction")], reference ~ prediction, value.var = "fraction")

rownames(df_table) = df_table[,"reference"]
df_table = df_table[,-1]

#' Unique mutual match between reference and predicted labels
#'
#' Takes a long-format dataframe of label-transfer results (one row per
#' reference/prediction pair) and keeps pairs that pass both thresholds
#' (fraction >= frac_min and mean_prob >= prob_min). A pair is returned only if
#' it is the sole passing pair for its reference and the sole passing pair for
#' its prediction; if any other cluster also passes both thresholds, the match
#' is considered ambiguous and dropped.
#'
#' @param df Dataframe with one row per reference/prediction pair.
#' @param row_col,col_col Names of the columns holding reference and prediction labels.
#' @param frac_col,prob_col Names of the columns holding the fraction of cells matched and the mean probability.
#' @param frac_min,prob_min Minimum fraction and mean probability (default 0.5).
#' @return The subset of rows of df that are unique mutual matches.
mutual_best_df <- function(df,
                           row_col  = "reference",
                           col_col  = "prediction",
                           frac_col = "fraction",
                           prob_col = "mean_prob",
                           frac_min = 0.5,
                           prob_min = 0.5)
{
  # 1. keep only pairs passing both thresholds
  hits <- df[df[[frac_col]] >= frac_min & df[[prob_col]] >= prob_min, , drop = FALSE]

  # 2. count passing partners for each reference and each prediction
  n_per_ref  <- ave(seq_len(nrow(hits)), hits[[row_col]], FUN = length)
  n_per_pred <- ave(seq_len(nrow(hits)), hits[[col_col]], FUN = length)

  # 3. keep pairs with exactly one passing partner on both sides
  hits[n_per_ref == 1 & n_per_pred == 1, , drop = FALSE]
}

mutual_best_df(df)
mutual_best_df(df, frac_min = 0.3)
mutual_best_df(df, frac_min = 0.3, prob_min=0.3)