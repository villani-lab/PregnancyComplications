code_location = "/projects/Pregnancy_placenta/Code_R/Figures/AUX_files/"


####################################################################################################################
#   Summarizing pathways regardless of cell subsets 
####################################################################################################################
  # 1) General counts - per pathway / condition / sign (up/ down) - how many cell subsets have significant enrichemnt
  #-------------------------------------------------------------------------------------------------------------------
  # summary_IPA_logP - for every condition / pathway / sign (up/ down) - counts number subsets-comparison that go are significant.
  # Note that every comparison (condition against any control) is counted separately. So if n=3 and there are only 2 cell types in the column "Cluster" this means that one of these clusters was significant in more than one comparison
  #-------------------------------------------------------------------------------------------------------------------
setwd("/projects/Pregnancy_placenta/Results/15_IPA/summarize_PWs_Andrea/init")
ipa_types=c("parent_lin_canonical_pathways", "parent_lin_diseases","sublin_canonical_pathways", "sublin_diseases")

p_thresh=0.25
for (ipa_type in ipa_types) # ipa_type="sublin_canonical_pathways"
{
  print(ipa_type)
  df_sublin = read.table(paste0(code_location,"/IPA_res_",ipa_type,".txt"), sep = "\t", header = TRUE, quote = "", comment.char = "")
  df = na.omit(df_sublin)
  for (p_z_pair in list(c(p_thresh, 2)))
  {
    log_p_thresh = -log10(p_z_pair[1])
    z_thresh = p_z_pair[2]
    df_filt = df[(df$log10_padj > log_p_thresh) & (abs(df$Zscore) > z_thresh),]
    df_filt$Z_sign = ifelse(df_filt$Zscore > 0, "Up", "Down")
    cell_df = count(df_filt, Canonical_Pathways, Condition, Z_sign, Cluster)
    count_df = count(cell_df, Canonical_Pathways, Condition, Z_sign) # how many cell subsets per PW+condition+sign
    # group by Canonical_Pathways, Condition and Z_sign and make a list of all the Clusters
    cells_list_df = cell_df %>% group_by(Canonical_Pathways, Condition, Z_sign) %>% summarise(Clusters = paste(Cluster, collapse = ","))
    final_df = merge(count_df, cells_list_df, by = c("Canonical_Pathways", "Condition", "Z_sign"))
    names(final_df)[4] = "Subset_number"
    write.xlsx(final_df, paste0("summary_",ipa_type,"_P",p_thresh,"_Z",z_thresh,".xlsx"), row.names = FALSE)
  }
}



# Merge canonical pathwyas and diseases file + adding annotations for pathways
#=================================================================================
setwd("/projects/Pregnancy_placenta/Results/15_IPA/summarize_PWs_Andrea/")

add_annotations = function(df)
{
  pw_ann = read.xlsx("all_pw_annotation.xlsx", sheetIndex = 1)
  df_new = merge(df, pw_ann[,c("Pathway", "Annotation")], by ="Pathway", all.x=TRUE)
   return(df_new)
}

for (x in c("TopSummary_sublin", "TopSummary_parent_lin", "summary_sublin", "summary_parent_lin"))
{
  print(x)
  df1 = read.xlsx(paste0("init/",x,"_diseases_P0.25_Z2.xlsx"), sheetIndex = 1)
  df2 = read.xlsx(paste0("init/",x,"_canonical_pathways_P0.25_Z2.xlsx"), sheetIndex = 1)
  df = rbind(df1, df2)
  colnames(df)[1]="Pathway"
  ann_df = add_annotations(df)
  # add annotaitons of lineages
  if (grepl("^summary.*",x)==TRUE)
  {
    ann_df$HBC = sapply(ann_df$Clusters, FUN=grepl, pattern='HBC')
    ann_df$Trophoblasts = sapply(ann_df$Clusters, FUN=grepl, pattern='VCT|SCT|EVT')
    ann_df$Monocytes = sapply(ann_df$Clusters, FUN=grepl, pattern='monocyte|Monocytes')
    ann_df$mMac = sapply(ann_df$Clusters, FUN=grepl, pattern='mMac|PAMM')
    ann_df$DC = sapply(ann_df$Clusters, FUN=grepl, pattern='cDC|DC')
    ann_df$B = sapply(ann_df$Clusters, FUN=grepl, pattern='B_cells')
    ann_df$Endothelial = sapply(ann_df$Clusters, FUN=grepl, pattern='Endo')
    ann_df$fFB = sapply(ann_df$Clusters, FUN=grepl, pattern='fFB')
    ann_df$fPV = sapply(ann_df$Clusters, FUN=grepl, pattern='fPV')
    ann_df$dS = sapply(ann_df$Clusters, FUN=grepl, pattern='dS')
  }
  write.xlsx(x=ann_df, file=paste0(x,"_P0.25_Z2.xlsx"))
}

# Creating supp table S9 - All IPA results
#----------------------------------------------------------
setwd("/projects/Pregnancy_placenta/Results/22_Final_figures/SuppTables/")
df_dis = read.table(paste0(code_location,"IPA_res_parent_lin_diseases.txt"), sep = "\t", header = TRUE, quote = "", comment.char = "")
df_canonical = read.table(paste0(code_location,"IPA_res_parent_lin_canonical_pathways.txt"), sep = "\t", header = TRUE, quote = "", comment.char = "")
df_dis["Canonical Pathway or Disease"] = "Disease"
df_canonical["Canonical Pathway or Disease"] = "Canonical Pathway"


