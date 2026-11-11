# Figure 3B) DA heatmap across conditions
library(openxlsx)
library(tidyr)
library(stringr)
library(ComplexHeatmap)
library(colorRamp2)

code_dir="/projects/Pregnancy_placenta/Code_Publication/R/"
source(paste0(code_dir, "Functions/DE_DA_functions.R"))
setwd("/projects/Pregnancy_placenta/Results/22_Final_figures/Figure_3/")

all_da_stat = read.csv("/projects/Pregnancy_placenta/Results/22_Final_figures/DE_DA_results/all_DA_stats.csv")
myeloid_da_stat = all_da_stat[all_da_stat$Lineage=="Myeloid",]
myeloid_list= prep_df_stat_for_plotting(myeloid_da_stat)

log2fc_int_list = list(Myeloid = myeloid_list$log2_df)
qvalue_int_list = list(Myeloid = myeloid_list$qvalue_df)
pvalue_int_list = list(Myeloid = myeloid_list$pvalue_df)
myeloid_order=c("HBC_1","HBC_2","HBC_3","HBC_4","HBC_5","HBC_6","HBC_7","HBC_8","HBC_9","HBC_10","HBC_11","HBC_12","HBC_13","Mon_1_CD14","Mon_2_CD14","mMon_4_MS1","fMon_3_CD16","mMac_1","mMac_2","mMac_3","mMac_4","mMac_5","mMac_6","mMac_7","fDC2","mDC1","fDC3","mDC2/4","mDC_mature","fDC","pDC")

row_split_by=data.frame(c("HBC","HBC","HBC","HBC","HBC","HBC","HBC","HBC","HBC","HBC","HBC","HBC","HBC","Mon","Mon","Mon","Mon","mMac","mMac","mMac","mMac","mMac","mMac","mMac","DCs","DCs","DCs","DCs","DCs","DCs","DCs"), row.names = c("HBC_1","HBC_2","HBC_3","HBC_4","HBC_5","HBC_6","HBC_7","HBC_8","HBC_9","HBC_10","HBC_11","HBC_12","HBC_13","Mon_1_CD14","Mon_2_CD14","mMon_4_MS1","fMon_3_CD16","mMac_1","mMac_2","mMac_3","mMac_4","mMac_5","mMac_6","mMac_7","fDC2","mDC1","fDC3","mDC2/4","mDC_mature","fDC","pDC"))
colnames(row_split_by)="split"
row_split_by[,"split"] = factor(row_split_by[,"split"], levels=c("HBC","Mon","mMac","DCs"), ordered = TRUE)

DA_heatmaps(lineages_of_interest="Myeloid",
                 log2FC_list=log2fc_int_list,
                 qval_list=qvalue_int_list,
                 pval_list=pvalue_int_list,
                 plot_row_order = TRUE,
                 desired_row_order=myeloid_order,  # New argument for custom order
                 #order_comp2 = c("Cont_labored", "Cont_unlabored", "Cont_PTB"),
                 plot_type = "both",
                 save_pdf = TRUE,
                 output_prefix = "3B_myeloid_DA_heatmap",
                 output_dir = ".",
                 pdf_width = 7,
                 pdf_height = 8,
                 row_split_by=row_split_by)


# ======================================================================
# Figure 3D
# ======================================================================

# Figure 3D) Pathway enrichment across conditions (shared pathways dotplot)
# ------------------------------------------------------------
# Recreate a dot-plot from summary_parent_lin_P0.25_Z2.xlsx  with explicit pathway selection (matches the figure)
# Group separators are based on the Annotation column.
# ------------------------------------------------------------

library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(patchwork)
library(grid)
library(gtable)
library(cowplot)

setwd("/projects/Pregnancy_placenta/Results/22_Final_figures/Figure_3/")
code_path = "/projects/Pregnancy_placenta/Code_R/Figures/"
summary_canon_pw <- read.table(paste0(code_path, "AUX_files/IPA_res_parent_lin_canonical_pathways.txt"))
summary_diseases <- read.table(paste0(code_path, "AUX_files/IPA_res_parent_lin_diseases.txt"))

pw_zscores <- rbind(summary_canon_pw, summary_diseases)
df <- read_excel("summary_parent_lin_P0.25_Z2.xlsx" , sheet = 1)

# ---- EXPLICIT pathway selection (the "fix") ----
# These are the pathways that appear in the figure.
selected_pathways <- c(
  # Tissue repair
    "Wound Healing Signaling Pathway",
    "Extracellular matrix organization",
  # Metabolic
    "Warburg Effect Signaling Pathway",
    "Synthesis of lipid",
    "Fatty acid metabolism",
  # Angiogenesis
    "Angiogenesis",
    "Cell surface interactions at the vascular wall",
    "Vasculogenesis",
  # Stress
    "Cellular response to heat stress",
    "Protein Ubiquitination Pathway",
    "Stress Granule Signaling Pathway",
    "HSP90 chaperone cycle for steroid hormone receptors in the presence of ligand",
    "HIF1α Signaling",
  # Pro-inflammatory
    "Viral Infection",
    "Pathogen Induced Cytokine Storm Signaling Pathway",
    "Antimicrobial response",
    "Immune response of cells",
    "Activation of leukocytes",
    "S100 Family Signaling Pathway",
  # Migration
    "Invasion of cells",
    "Migration of cells",
    "Cell movement of phagocytes",
    "Chemotaxis of phagocytes",
    "Leukocyte migration",
  # Cell survival
    "Organismal death",
    "Apoptosis",
    "Cell survival",
    "Cell viability")


# Summary info
#---------------
df <- df %>%  filter(Pathway %in% selected_pathways)
# Keep pathway order exactly as in selected_pathways (bottom->top will be handled by factor order)
conditions_order <- c("PTB", "PrePE", "TermPE", "FGR","DM1")
df <- df %>%
  mutate(Condition = factor(Condition, levels = conditions_order),
         Pathway   = factor(Pathway, levels = selected_pathways))

# ---- Cluster definitions ----
all_clusters  <- c("VCT", "SCT","HBC","PAMM", "Monocytes", "DC", "Tcells", "NK","fFB", "Endo_f")
color_map = c(VCT = "#166131", SCT="#94cd99", PAMM="#c54527",HBC="#6f2723", Monocytes="#f49582", DC="#e2bf86", Tcells="#ef8e1b", NK="#efcd3d", fFB="#4090ce", Endo_f="#7c52a2")

# ---- Prep plotting dataframe ----
# creating a row per cluster/ comparison/ pw only for relevant pathways and clusters
plot_df <- df %>%
  mutate(Cluster = str_split(Clusters, ",")) %>%
  unnest(Cluster) %>%
  mutate(Cluster = str_trim(Cluster)) %>%
  filter(Cluster %in% c(all_clusters)) %>%
  mutate(Lineage = case_when(
                  Cluster %in% all_clusters ~ "all_clusters",
                  TRUE ~ NA_character_),
        Z_sign = factor(Z_sign, levels = c("Up", "Down")),
        Cluster = factor(Cluster, levels = names(color_map))) %>%
  filter(!is.na(Lineage))

# Z-score info
#---------------
pw_zscores["Cluster"] <- sapply(pw_zscores[,"Clust_Comparison"], FUN =  function(x) strsplit(x, "-", fixed = TRUE)[[1]][1])
zscores_filt <- pw_zscores[(pw_zscores$Cluster %in% all_clusters) & (pw_zscores$Canonical_Pathways %in% selected_pathways),]

# Is the Z-score the same direction in the same cell type and complication across controls?
# lets have a circle if the pw is significant in different directions in different controls, and a triangle towards its directionality if it is consistent across comprisons in which it is significant
# Any condition that it is not significant (Z<2, adj_p<0.05)
pick_z <- function(z) {
  z <- z[!is.na(z)]
  if (length(z) == 0) return(NA_real_)
  s <- sign(z)
  s <- s[s != 0]                 # ignore exact zeros (optional)
  if (length(s) == 0) return(0)  # all zeros
  if (length(unique(s)) > 1) return(0)
    z[which.max(abs(z))]           # most extreme (largest |Z|)
}

#padj_thresh = 0.2
Z_thresh = 2
all_zscore = data.frame()
for (clust in all_clusters)
  for (pw in selected_pathways)
    for (cond in conditions_order)
    {
      curr = zscores_filt %>% filter(Cluster == clust) %>% filter(Canonical_Pathways == pw) %>% filter(Condition == cond)
      #sig_curr = curr %>% filter(log10_padj > -log10(padj_thresh)) %>% filter(abs(Zscore) > Z_thresh)
      sig_curr = curr %>%  filter(abs(Zscore) > Z_thresh)
      if (dim(sig_curr)[1]==1)  # if there is just one sig result - take it's Zscore
        all_zscore = rbind(all_zscore, data.frame(PW = pw, Condition = cond, Cluster = clust, Z = sig_curr[,"Zscore"]))
      else
        if (dim(sig_curr)[1]>1) # if there is more than one sig result - check if they all agree
        {
          #print(paste0(clust," ",pw," ",cond))
          #print(sig_curr)
          z_out <- pick_z(curr$Zscore)
          all_zscore = rbind(all_zscore, data.frame(PW = pw, Condition = cond, Cluster = clust, Z = z_out))
        }
    }

# Add Z to plot_df and check missing combos
all_zscore2 <- all_zscore %>% rename(Pathway = PW)
plot_df <- plot_df %>% left_join(all_zscore2, by = c("Pathway", "Condition", "Cluster"))
write.csv(plot_df, "fig3D_allZscores.csv")

plot_df <- plot_df %>%
  mutate(
    Condition = factor(Condition, levels = conditions_order),
    Pathway   = factor(Pathway, levels = selected_pathways),
    Cluster   = factor(Cluster, levels = names(color_map))
  )

missing_z <- plot_df %>%
  distinct(Pathway, Condition, Cluster, .keep_all = TRUE) %>%
  filter(is.na(Z))

if (nrow(missing_z) > 0)
{
  message("Missing Z-scores for ", nrow(missing_z), " Pathway×Condition×Cluster combinations. Examples:")
  print(head(missing_z %>% select(Pathway, Condition, Cluster), 25))
} else {message("All Pathway×Condition×Cluster combinations in plot_df have a Z-score in all_zscore.")}

# Seeting triangle direction + size based on |Z|
plot_df <- plot_df %>%
  mutate(
    Z_dir = case_when(is.na(Z) ~ "missing", Z > 0 ~ "up", Z < 0 ~ "down", TRUE ~ "zero"),
    Z_abs = ifelse(is.na(Z), NA_real_, abs(Z))
  )

# ---- Compute separator lines based ONLY on Annotation ----
# A dashed line is drawn between pathways whenever Annotation changes
annotation_df <- plot_df %>%
  distinct(Pathway, Annotation) %>%
  arrange(Pathway) %>%
  mutate(y = as.numeric(Pathway))

sep_y <- annotation_df %>%
  mutate(annotation_next = lead(Annotation)) %>%
  filter(!is.na(annotation_next) & Annotation != annotation_next) %>%
  transmute(yintercept = y + 0.5) %>%
  pull(yintercept)


# ---- Panel function (manual horizontal spread) ----
# cond_gap: when >1 increases space between conditions
# within_scale: when <1 squeezes dots within condition
make_panel <- function(dat, title, show_y = TRUE, spread = 0.18,cond_gap = 1.4, within_scale = 1 )
{
  dat2 <- dat %>%
    mutate(x_base = as.numeric(Condition) * cond_gap) %>%
    group_by(Condition, Pathway) %>%
    arrange(Cluster, .by_group = TRUE) %>%
    mutate(
      n_in_cell = n(),
      idx_in_cell = row_number(),
      x_offset = (idx_in_cell - (n_in_cell + 1) / 2) * spread * within_scale,
      x_pos = x_base + x_offset
    ) %>%
    ungroup()

  cond_levels <- levels(dat2$Condition)
  ncond <- length(cond_levels)

  p <- ggplot(dat2, aes(x = x_pos, y = Pathway)) +
    geom_hline(yintercept = sep_y, linetype = "dashed", linewidth = 1, alpha = 0.7) +
    geom_point(data = dat2 %>% filter(!is.na(Z) & Z != 0), aes(color = Cluster, fill= Cluster, shape = Z_dir, size = Z_abs), stroke = 1.1) +  # triangles for Z (up/down), size ~ |Z|; keep missing as circles
    geom_point(data = dat2 %>% filter(is.na(Z) | Z == 0), aes(color = Cluster), fill = NA, size = 3.5, stroke = 1.1) +
    scale_shape_manual(values = c(up = 24, down = 25)) +  # upward facing tringle for Z>0 downward facing triangle for Z<0
    scale_size(range = c(1.5, 8), name = "|Z-score|", breaks = c(2, 4, 6, 8),  limits = c(2, NA)) +
    scale_x_continuous(
      breaks = (seq_len(ncond) * cond_gap), labels = cond_levels,
      limits = c(0.5 * cond_gap, (ncond + 0.5) * cond_gap)) +
    scale_y_discrete(drop = FALSE) +
    scale_color_manual(values = color_map, drop = FALSE) +
    scale_fill_manual(values = color_map, drop = FALSE, guide = "none")+
    guides(color = guide_legend(title = "Cell type"),
    size = guide_legend(override.aes = list(shape = 24, fill = "grey70", color = "grey30"))) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, face = "bold", size = 16),
      axis.text.y = element_text(size = 12),   # pathway font size
      panel.grid.minor = element_blank(),
      legend.position = "top",
      legend.box = "vertical",
      plot.margin = margin(0, 30, 0, 0, unit = "pt")
    )

  if (!show_y) {
    p <- p + theme(
      axis.title.y = element_blank(),
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank()
    )
  }
  p
}

make_panel(plot_df, "All_cells", show_y = TRUE, spread = 0.18)






# ---- Save ----
ggsave("Figure_3D_triangles.pdf", make_panel(plot_df, "All_cells", TRUE, 0.18), width = 15, height = 10)

