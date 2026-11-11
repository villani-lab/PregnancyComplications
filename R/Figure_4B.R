# Figure 4B) PTB pathway enrichment heatmap
# Heatmap of IPA pathways enrichemnt results in PTB vs cont_PTB across cell types

code_dir="/projects/Pregnancy_placenta/Code_Publication/R/"
source(paste0(code_dir, "Functions/figure_functions.R"))
code_path = "/projects/Pregnancy_placenta/Code_R/Figures/"
setwd("/projects/Pregnancy_placenta/Results/22_Final_figures/Figure_4/")

# Merged - canonical pathways and diseases and biofunctions
df_canonical = read.table(paste0(code_path, "/AUX_files/IPA_res_sublin_canonical_pathways.txt"), sep="\t", header =1)
df_disease = read.table(paste0(code_path, "/AUX_files/IPA_res_sublin_diseases.txt"), sep="\t", header =1)
df = rbind(cbind(df_canonical, Type="Canonical_Pathways"), cbind(df_disease, Type="Diseases_and_biofunctions"))
PTB_df = df[(df$Comparison == "PTB_vs_Cont_PTB"),]

plot_heatmap_PW(df_comp_merged=PTB_df,
                        zscore_min_cp = 2,
                        padj_thresh_cp = 0.01,
                        zscore_min_db = 3,
                        padj_thresh_db = 0.0001,
                        save_pdf = TRUE,
                        output_prefix = paste0("Fig4B_PTB_PW_heatmap.pdf"),
                        order_vector=c("HBC_1","HBC_10","HBC_11","HBC_13","HBC_2","HBC_4","HBC_6","HBC_9","CD14_monocyte","CD16_monocyte","cDC","mMac_upper","fPV2","fPV3","mEndo_1","VCT_cycling_1","VCT_cycling_2","VCTfusion_1","VCT_1"),
                        orderColBy = "Cells",
                        pdf_width = 15,
                        pdf_height = 24,
                        cluster_names_col="Cluster",
                        filter_more="N",
                        plot_title="PTB vs PTB_cont",
                        pw_annot_file = paste0(code_path, "/AUX_files/PTB_merged_Pathway_annotation_v5.xlsx"))

