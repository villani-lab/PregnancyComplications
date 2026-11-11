# Running DE analysis on all lineages across 13 comparisons

library(readxl)
library(ggplot2)
require(ggh4x)
require(forcats)
library(ComplexHeatmap)
library(reshape2)
library(circlize)

code_dir = "/projects/Pregnancy_placenta/Code_R/"
source(paste0(code_dir, "/Functions/DE_DA_functions.R"))

#                                              DESeq2
#----------------------------------------------------------------------------------------------------------------
# Running DESeq2 on every cluster in every lineage for every comparison
# Then compiling all these DEGs into a single dile
#----------------------------------------------------------------------------------------------------------------
pairs_test <- data.frame(controls_comp=c("Cont_labored", "Cont_unlabored"),
                         PTB_lbr=c("PTB","Cont_labored"),
                         DM1_lbr=c("DM1","Cont_labored"),
                         FGR_lbr=c("FGR","Cont_labored"),
                         PrePE_lbr=c("PrePE","Cont_labored"),
                         TermPE_lbr = c("TermPE","Cont_labored"),
                         DM1_unlbr=c("DM1","Cont_unlabored"),
                         FGR_unlbr=c("FGR","Cont_unlabored"),
                         PrePE_unlbr=c("PrePE","Cont_unlabored"),
                         TermPE_unlbr = c("TermPE","Cont_unlabored"),
                         PTB_lbr=c("PTB","Cont_unlabored"),
                         PrePE_ptb = c("PrePE", "Cont_PTB"),
                         PTB_contPTB=c("PTB","Cont_PTB") )
lin_dict = c("B_Plasma"="B_plasma.h5mu",
            "DCs"="DCs_6subsets.h5mu",
            "dS"="dS.h5mu",
            "Endo"="Endo.h5mu",
            "EVT"="EVT.h5mu",
            "fFB"="fFB.h5mu",
            "fPV"="fPV.h5mu",
            "HBC"="HBC_13subsets.h5mu",
            "Mon_mMac"="mMac_Mon_12subsets.h5mu",
            "NK"="NK.h5mu",
            "T"="T_cells.h5mu",
            "Trophoblasts"="Trophoblasts.h5mu")

setwd("/projects/Pregnancy_placenta/Results/23_Supp_tables/DEG_comparisons/")
samp_condition_file = "/projects/Pregnancy_placenta/Results/13_DE/sample_condition_list.csv"

for (lin in names(lin_dict))
{
  pb_path = "/projects/Pregnancy_placenta/Results/23_Supp_tables/DEG_comparisons/"
  dir.create(lin)
  setwd(lin)
  DE_pregnancy_wrapper(pb_counts_file= paste0(pb_path, lin, "_pb_counts.csv"),
                       meta_file= paste0(pb_path, lin, "_pb_meta.csv"),
                       pairs_test = pairs_test,
                       cell_counts_file= paste0(pb_path, lin, "_pb_cells_per_samp_clust.csv"),
                       samp_condition_file=samp_condition_file,
                       run_title=lin,
                       cluster_col='cluster',
                       plot_heatmap=TRUE)

  summarize_DE_results(path=".",
                       title=lin,
                       comparisons = sapply(pairs_test, function(x) {paste0(x[1], "_vs_", x[2])}),
                       lineages=c(lin))
  setwd("..")
}


#                                                    Summary
#-------------------------------------------------------------------------------------------------------------------
# Creating a supplementary table that contains the top 200 genes by abs(log(FC)) that have a nominal p-value of 0.01
#-------------------------------------------------------------------------------------------------------------------
df_summary = data.frame()
for (lin in names(lin_dict))
{
  print(lin)
  lin_summary = read.table(paste0("/projects/Pregnancy_placenta/Results/23_Supp_tables/DEG_comparisons/", lin, "/AllDEGs_",lin,".csv"), sep=",", header = TRUE)
  clusters = unique(lin_summary[,"cluster"])
  comparisons = unique(lin_summary[,"contrast"])
  for (clust in clusters)
  {
      for (comp in comparisons)
      {
        curr_df = lin_summary[(lin_summary["contrast"]==comp) & (lin_summary["cluster"]==clust), ]
        curr_df_filt = curr_df[curr_df["pvalue"] <= 0.05,]
        curr_df_filt = na.omit(curr_df_filt)
        curr_df_filt = curr_df_filt[order(curr_df_filt$abslog2FoldChange),]
        df_summary = rbind(df_summary, curr_df_filt[0:200,])
      }
  }
}





