library(ggplot2)
library(reshape2)
library(xlsx)
library(ComplexHeatmap)
library(tidyr)
library(writexl)
library(stringr)
library(dplyr)

code_path="/projects/Pregnancy_placenta/Code_R/Figures/"
code_path_PAL="/projects/Pregnancy_placenta/Code_PAL/Preg_placenta_R_codes/Figures"
code_dir="/projects/Pregnancy_placenta/Code_Publication/R/"
source(paste0(code_dir, "Functions/figure_functions.R"))
setwd("/projects/Pregnancy_placenta/Results/22_Final_figures_PAL/SupFig_DM1_FGR/")



#========================================================
# modify PW Files
#========================================================
pw_annot_file_PE <- read_excel(paste0(code_path,"/AUX_files/PE_merged_PW_annotation_v6.xlsx")) %>%
  select(Pathway, Annotation) %>%
  mutate(Annotation = na_if(Annotation, "NA"))
pw_annot_file_PE = pw_annot_file_PE[,c("Pathway","Annotation")]
pw_annot_file_PE$Annotation=="NA"

pw_annot_file_PTB = read_excel(paste0(code_path,"/AUX_files/PTB_merged_Pathway_annotation_v5.xlsx")) %>%
  select(Pathway, Annotation) %>%
  mutate(Annotation = na_if(Annotation, "NA"))
pw_annot_file_PTB = pw_annot_file_PTB[,c("Pathway","Annotation")]

# Rename annotation columns to track origin
pw_PE  <- pw_annot_file_PE  %>% rename(Annotation_PE  = Annotation)
pw_PTB <- pw_annot_file_PTB %>% rename(Annotation_PTB = Annotation)

# Full join on Pathway
merged_pw <- full_join(pw_PE, pw_PTB, by = "Pathway") %>%
  arrange(Pathway)

## LOADING MY OWN MERGED DM1 results
merged_pw = read_excel(paste0(code_path_PAL,"/AUX_files/DM1_merged_PW_annotation_v3.xlsx"))%>%
  select(Pathway, Annotation) %>%
  mutate(Annotation = na_if(Annotation, "NA"))
merged_pw= merged_pw[,c("Pathway","Annotation")]

#========================================================
# SupFig X- pathway enrichment in DM1
#========================================================

df_canonical = read.table(paste0(code_path, "AUX_files/IPA_res_sublin_canonical_pathways.txt"), sep="\t", header =1)
df_disease = read.table(paste0(code_path, "AUX_files/IPA_res_sublin_diseases.txt"), sep="\t", header =1)

df = rbind(cbind(df_canonical, Type="Canonical_Pathways"), cbind(df_disease, Type="Diseases_and_biofunctions"))
DM1_df <- df[df$Comparison %in% "DM1_vs_Cont_unlabored",]
DM1_df[DM1_df["Lineage"]=="HBC","Lineage"] <- "Myeloid"


# get unique DM1 pathways
dm1_paths <- unique(DM1_df$Canonical_Pathways)

# add in_DM1 column to merged_pw
merged_pw <- merged_pw %>%
  mutate(in_DM1 = ifelse(Pathway %in% dm1_paths, "Y", "N"))

# find DM1 pathways not already present in merged_pw
missing_dm1 <- setdiff(dm1_paths, merged_pw$Pathway)

# make rows to append
missing_rows <- data.frame(
  Pathway = missing_dm1,
  #Annotation_PE = NA,
  #Annotation_PTB = NA,
  Annotation = NA,
  in_DM1 = "Y",
  stringsAsFactors = FALSE
)

# append missing rows
merged_pw <- bind_rows(merged_pw, missing_rows) %>%
  arrange(Pathway)

# merge columns, prioritie PE vs PTB

## ignoring this for now as we want to see which pathways appear but not annoted
#merged_pw_sub <- merged_pw %>%
#  filter(
#    in_DM1 == "Y",
#    !(is.na(Annotation_PE) & is.na(Annotation_PTB))
#  )
merged_pw_sub = merged_pw

#merged_pw_sub <- merged_pw_sub %>%
#  mutate(
    #Annotation = coalesce(Annotation_PE, Annotation_PTB)
#  )

#adding this to see whats not annotated
merged_pw_sub[is.na(merged_pw_sub$Annotation),]$Annotation = "Not_Annotated"

merged_pw_sub = merged_pw_sub[,c("Pathway","Annotation")]


write_xlsx(merged_pw_sub,paste0(code_path_PAL,"/AUX_files/DM1_merged_PW_annotation_v2_PLUS_NOT_ASSIGNED.xlsx") )

# ----- RUN pathway heatmap -----
plot_heatmap_PW(df_comp_merged=DM1_df,
                        zscore_min_cp =1,
                        padj_thresh_cp = 0.25,
                        zscore_min_db = 2,
                        padj_thresh_db = 0.05,
                        save_pdf = TRUE,
                        #output_prefix = paste0("SupFig_DM1_both_adjPval_0.25_v1"),
                output_prefix = paste0("SupFig_DM1_both_adjPval_0.05_ALL_PATHS_V5_AE"),
                        orderColBy = "Cells",
                       # order_vector=c("HBC_2_TermPE","HBC_3_TermPE","CD14_monocyte_TermPE","HBC_4_TermPE","HBC_13_TermPE","HBC_2_PrePE","HBC_3_PrePE","CD14_monocyte_PrePE","HBC_1_PrePE","SCT_TermPE","VCT_1_TermPE","VCT_3_TermPE","VCTfusion_1_TermPE","VCTfusion_2_TermPE","VCT_cycling_1_TermPE","VCT_cycling_2_TermPE","SCT_PrePE","VCT_1_PrePE","VCT_3_PrePE","VCTfusion_1_PrePE","VCTfusion_2_PrePE","VCTccc_PrePE","fEndo_1_PrePE","fEndo_2_PrePE","fEndo_3_PrePE"),
                        pdf_width = 15,
                      split_cols = c(DM1="black"),
                        pdf_height = 28,
                        cluster_names_col="Cluster",
                        filter_more="N",
                row_display_mode = "pathway_annotation_group",
                        #plot_title="PE vs unlabored_cont",
                       pw_annot_file = paste0(code_path_PAL,"/AUX_files/DM1_merged_PW_annotation_v4.xlsx"))


#========================================================
# SupFig X- pathway enrichment in FGR
#========================================================
merged_pw = read_excel(paste0(code_path_PAL,"/AUX_files/DM1_merged_PW_annotation_v3.xlsx"))%>%
  select(Pathway, Annotation) %>%
  mutate(Annotation = na_if(Annotation, "NA"))
merged_pw= merged_pw[,c("Pathway","Annotation")]

df_canonical = read.table(paste0(code_path, "AUX_files/IPA_res_sublin_canonical_pathways.txt"), sep="\t", header =1)
df_disease = read.table(paste0(code_path, "AUX_files/IPA_res_sublin_diseases.txt"), sep="\t", header =1)

df = rbind(cbind(df_canonical, Type="Canonical_Pathways"), cbind(df_disease, Type="Diseases_and_biofunctions"))
FGR_df <- df[df$Comparison %in% c("FGR_vs_Cont_unlabored","FGR_vs_Cont_labored"),]
FGR_df [FGR_df ["Lineage"]=="HBC","Lineage"] <- "Myeloid"


# get unique DM1 pathways
fgr_paths <- unique(FGR_df$Canonical_Pathways)

# add in_DM1 column to merged_pw
merged_pw_FGR <- merged_pw %>%
  mutate(in_FGR = ifelse(Pathway %in% fgr_paths, "Y", "N"))

# find FGR pathways not already present in merged_pw
missing_fgr <- setdiff(fgr_paths, merged_pw_FGR$Pathway)

# make rows to append
missing_rows <- data.frame(
  Pathway = missing_fgr,
  #Annotation_PE = NA,
  #Annotation_PTB = NA,
  in_FGR= "Y",
  Annotation = NA,
  stringsAsFactors = FALSE
)

# append missing rows
merged_pw_FGR <- bind_rows(merged_pw_FGR, missing_rows) %>%
  arrange(Pathway)

# merge columns, prioritie PE vs PTB
#merged_pw_FGR_sub <- merged_pw_FGR %>%
#  filter(
#    in_FGR == "Y",
#    !(is.na(Annotation_PE) & is.na(Annotation_PTB))
#  )
merged_pw_FGR_sub = merged_pw_FGR
#merged_pw_FGR_sub <- merged_pw_FGR_sub %>%
#  mutate(
#    Annotation = coalesce(Annotation_PE, Annotation_PTB)
#  )

#merged_pw_FGR_sub = merged_pw_FGR_sub[,c("Pathway","Annotation")]
merged_pw_FGR_sub[is.na(merged_pw_FGR_sub$Annotation),]$Annotation = "Not_Annotated"

merged_pw_FGR_sub = merged_pw_FGR_sub[,c("Pathway","Annotation")]

write_xlsx(merged_pw_sub,paste0(code_path_PAL,"/AUX_files/FGR_merged_PW_annotation_v3.xlsx") )


# ----- RUN pathway heatmap -----
plot_heatmap_PW(df_comp_merged=FGR_df,
                        zscore_min_cp = 2,
                        padj_thresh_cp = 0.25,
                        zscore_min_db = 3,
                        padj_thresh_db = 0.001,
                        save_pdf = TRUE,
                        #output_prefix = paste0("SupFig_FGR_CP_adjPval_0.05_DB_adjPval_0.001_test4"),
                      output_prefix = paste0("SupFig_FGR_SMFM_thresh_V1"),
                        orderColBy = "Cells",
                       # order_vector=c("HBC_2_TermPE","HBC_3_TermPE","CD14_monocyte_TermPE","HBC_4_TermPE","HBC_13_TermPE","HBC_2_PrePE","HBC_3_PrePE","CD14_monocyte_PrePE","HBC_1_PrePE","SCT_TermPE","VCT_1_TermPE","VCT_3_TermPE","VCTfusion_1_TermPE","VCTfusion_2_TermPE","VCT_cycling_1_TermPE","VCT_cycling_2_TermPE","SCT_PrePE","VCT_1_PrePE","VCT_3_PrePE","VCTfusion_1_PrePE","VCTfusion_2_PrePE","VCTccc_PrePE","fEndo_1_PrePE","fEndo_2_PrePE","fEndo_3_PrePE"),
                        pdf_width = 15,
                      split_cols = c(Cont_labored="Grey", Cont_unlabored="black"),
                        pdf_height = 28,
                        cluster_names_col="Cluster",
                        filter_more="N",
                row_display_mode = "pathway_annotation_group",
                split_by_cond_OR_ctrl = "Ctrl",
                        #plot_title="PE vs unlabored_cont",
                       pw_annot_file = paste0(code_path_PAL,"/AUX_files/FGR_merged_PW_annotation_v3.xlsx"))

