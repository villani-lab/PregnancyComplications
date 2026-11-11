# Summarizing pathways regardless of cell subsets
library(ggplot2)
library(reshape2)
library(xlsx)
library(ComplexHeatmap)
library(tidyr)
library(stringr)
library(dplyr)
source("/projects/Pregnancy_placenta/Code_R/Functions_IPA.R")

setwd("/projects/Pregnancy_placenta/Results/20_Figures/Pathways/")

#========================================================================
# 1) Summarize parent sublineages - Canonical Pathways / Diseases
#========================================================================
pw_dir = "/projects/Pregnancy_placenta/Code_PAL/DE_IPA_Code_Rmarkdown_PAL/IPA_results/Parent_lin/"
for (ipa_type in c("canonical_pathways", "diseases"))
{
  print(ipa_type)
  df_parent = data.frame()
  for(comp_no in c("1","2","3"))
  {
    df_curr = prep_IPA_df(dir_to_pathway = paste0(pw_dir,"/Parent_lin_test3_CompNum_",comp_no,"/"), pattern_int = ipa_type)
    df_curr$CompNum = comp_no
    df_parent = rbind(df_parent, df_curr)
  }
  write.table(df_parent, paste0("IPA_res_parent_lin_",ipa_type,".txt"))
}

# subset populations that don't have a "sub_lin" title
parent_sub = subset(df_parent,Cluster %in% c("SCT","VCTccc","Endo_mL","Endo_fProlif"))


# Summarize sublineages - Canonical Pathways / Diseases
#---------------------------------
pw_dir = "/projects/Pregnancy_placenta/Code_PAL/DE_IPA_Code_Rmarkdown_PAL/IPA_results/sub_lin/"

df_sublin = data.frame()
for(comp_no in c("1","2","3","4","5"))
{
  df_sublin_curr = prep_IPA_df(dir_to_pathway = paste0(pw_dir,"/Sublin_test3_CompNum_",comp_no,"/"), pattern_int = ipa_type)
  df_sublin_curr$CompNum = comp_no
  df_sublin = rbind(df_sublin, df_sublin_curr)
}

# Add cell subsets from parent  that are missing a "sub_lin"
df_sublin = rbind(df_sublin, parent_sub)


# Add lineage info to the shared_dist df
df_sublin["Lineage"] = "-"
df_sublin[grep("VCT|SCT",df_sublin$Cluster), "Lineage"] = "Trophoblasts"
df_sublin[grep("HBC",df_sublin$Cluster), "Lineage"] = "HBC"
df_sublin[grep("monocyte|cDC|mMac",df_sublin$Cluster), "Lineage"] = "Myeloid"
df_sublin[grep("Endo",df_sublin$Cluster), "Lineage"] = "Endothelial"
df_sublin[grep("dS|fPV|fFB",df_sublin$Cluster), "Lineage"] = "Stromal"

write.table(df_sublin, paste0("IPA_res_sublin_",ipa_type,".txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)




##############
# Summarize mid-level for endothelial cells and T cells
##############
ipa_type="diseases"
pw_dir = "/projects/Pregnancy_placenta/Code_PAL/DE_IPA_Code_Rmarkdown_PAL/IPA_results/Parent_lin/"

df_parent = data.frame()
for(comp_no in c("1","2","3"))
{
  df_curr = prep_IPA_df(dir_to_pathway = paste0(pw_dir,"/Parent_lin_test3_CompNum_",comp_no,"/"), pattern_int = ipa_type)
  df_curr$CompNum = comp_no
  df_parent = rbind(df_parent, df_curr)
}

# subset populations that don't have a "sub_lin" title
parent_sub = subset(df_parent,Cluster %in% c("SCT","VCTccc","Endo_mL","Endo_fProlif"))



