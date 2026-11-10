## Prep libraries----
source("/projects/Pregnancy_placenta/Code_R/Figures/figure_functions.R")
setwd("/projects/Pregnancy_placenta/Results/20_Figures/Pathways/PTB_vs_PTBcont_FC/")

###############################################################
#                          PTB
###############################################################
plot_heatmap_genes(PW_list="Organization of cytoskeleton",
                   condition= "PTB",
                   control= "Cont_PTB",
                   cell_subsets=c("fPV2","fPV3"),
                   width=4, height=15,
                   pdf_title="Organization of cytoskeleton PV")

plot_heatmap_genes(PW_list="Angiogenesis",
                   condition= "PTB",
                   control= "Cont_PTB",
                   cell_subsets=c("HBC_1","HBC_4_IL1B_MMP19","fPV2","fPV3"),
                   width=4, height=15)

plot_heatmap_genes(PW_list="Angiogenesis",
                   condition= "PTB",
                   control= "Cont_PTB",
                   cell_subsets=c("Endo_m1","HBC_1","HBC_4_IL1B_MMP19","mMac_upper"),
                   width=4, height=15,
                   pdf_title="Angiogenesis_endo")
plot_heatmap_genes(PW_list="Extracellular matrix organization",
                   condition= "PTB",
                   control= "Cont_PTB",
                   cell_subsets=c("Endo_m1"),
                   width=4, height=15,
                   pdf_title="ECM_endo")


plot_heatmap_genes(PW_list="Angiogenesis",
                   condition= "PTB",
                   control= "Cont_PTB",
                   cell_subsets=c("HBC_1","HBC_4_IL1B_MMP19","fPV2","fPV3","Endo_m1","mMac_upper","cDC"),
                   width=4, height=15, suffix="_full")

plot_heatmap_genes(PW_list="Vasculogenesis",
                   condition= "PTB",
                   control= "Cont_PTB",
                   cell_subsets=c("HBC_1","HBC_4_IL1B_MMP19","fPV2","fPV3","Endo_m1","cDC"),
                   width=4, height=15, suffix="_full")

plot_heatmap_genes(PW_list=c("Vasculogenesis","Angiogenesis"),
                   condition= "PTB",
                   control= "Cont_PTB",
                   cell_subsets=c("HBC_1","HBC_4_IL1B_MMP19","fPV2","fPV3","Endo_m1","cDC"),
                   width=4, height=15, suffix="_full")

plot_heatmap_genes(PW_list=c("Angiogenesis"),
                   condition= "PTB",
                   control= "Cont_PTB",
                   cell_subsets=c("HBC_1","HBC_4_IL1B_MMP19","cDC","mMac_upper","fPV2","fPV3","Endo_m1"),
                   width=4, height=15, suffix="_sigOnly", onlySig=TRUE)


plot_heatmap_genes(PW_list="PPAR Signaling",
                   condition= "PTB",
                   control= "Cont_PTB",
                   cell_subsets=c("HBC_1","HBC_4_IL1B_MMP19","HBC_6","HBC_11","HBC_2_CXCL2_CXCL3","HBC_9_CXCL10_IFIT44","HBC_10_cycling","mixed_HBC_13_S100A8"),
                   width=4, height=5)

plot_heatmap_genes(PW_list="FXR/RXR Activation",
                   condition= "PTB",
                   control= "Cont_PTB",
                   cell_subsets=c("HBC_1","HBC_4_IL1B_MMP19","HBC_6","HBC_11","HBC_2_CXCL2_CXCL3","HBC_9_CXCL10_IFIT44","HBC_10_cycling","mixed_HBC_13_S100A8"),
                   width=4, height=5)


plot_heatmap_genes(PW_list=c("IL-17 Signaling","IL-33 Signaling Pathway","IL-6 Signaling","IL-8 Signaling","Interleukin-1 family signaling","Interleukin-4 and Interleukin-13 signaling"),
                   condition= "PTB",
                   control= "Cont_PTB",
                   cell_subsets=c("HBC_1","HBC_4_IL1B_MMP19","HBC_6","HBC_11","HBC_2_CXCL2_CXCL3","HBC_9_CXCL10_IFIT44","HBC_10_cycling","mixed_HBC_13_S100A8"),
                   width=4, height=5,
                   pdf_title="IL_signaling",
                   cluster_rows = TRUE)

plot_heatmap_genes(PW_list=c("Pathogen Induced Cytokine Storm Signaling Pathway"),
                   condition= "PTB",
                   control= "Cont_PTB",
                   cell_subsets=c("HBC_1","HBC_4_IL1B_MMP19","HBC_6","HBC_11","HBC_2_CXCL2_CXCL3","HBC_9_CXCL10_IFIT44","HBC_10_cycling","mixed_HBC_13_S100A8","mMac_upper"),
                   width=4, height=5,
                   pdf_title="",
                   cluster_rows = TRUE)


plot_heatmap_genes(PW_list=c("Activation of lymphatic system cells"),
                   condition= "PTB",
                   control= "Cont_PTB",
                   cell_subsets=c("HBC_4_IL1B_MMP19","HBC_6","HBC_2_CXCL2_CXCL3"),
                   width=4, height=5,
                   pdf_title="",
                   cluster_rows = TRUE)

plot_heatmap_genes(PW_list=c("Wound Healing Signaling Pathway"),
                   condition= "PTB",
                   control= "Cont_PTB",
                   cell_subsets=c("HBC_1","HBC_4_IL1B_MMP19","fPV3","Endo_m1"),
                   width=4, height=5,
                   pdf_title="",
                   cluster_rows = TRUE)

plot_heatmap_genes(PW_list=c("Apoptosis"),
                   condition= "PTB",
                   control= "Cont_PTB",
                   cell_subsets=c("fPV3","fPV2"),
                   width=4, height=15,
                   pdf_title="",
                   cluster_rows = TRUE)

plot_heatmap_genes(PW_list=c("Migration of cells"),
                   condition= "PTB",
                   control= "Cont_PTB",
                   cell_subsets=c("fPV3","fPV2"),
                   width=4, height=15,
                   pdf_title="",
                   cluster_rows = TRUE)

plot_heatmap_genes(PW_list=c("Apoptosis"),
                   condition= "PTB",
                   control= "Cont_PTB",
                   cell_subsets=c("VCTfusion_1","VCT_1","VCTp_1","VCTp_2"),
                   width=4, height=15,
                   pdf_title="Apoptosis_trophoblasts",
                   cluster_rows = TRUE)

plot_heatmap_genes(PW_list=c("Invasion of cells"),
                   condition= "PTB",
                   control= "Cont_PTB",
                   cell_subsets=c("VCTfusion_1","VCT_1","VCTp_1","VCTp_2"),
                   width=4, height=15,
                   pdf_title="",
                   cluster_rows = TRUE)

plot_heatmap_genes(PW_list=c("Cellular response to heat stress"),
                   condition= "PTB",
                   control= "Cont_PTB",
                   cell_subsets=c("HBC_1","HBC_4_IL1B_MMP19","HBC_6","HBC_11","HBC_2_CXCL2_CXCL3","HBC_9_CXCL10_IFIT44","HBC_10_cycling","mixed_HBC_13_S100A8","mMac_upper"),
                   width=4, height=8,
                   pdf_title="",
                   cluster_rows = TRUE)

plot_heatmap_genes(PW_list=c("Cellular response to heat stress","HIF1α Signaling","HSP90 chaperone cycle for steroid hormone receptors in the presence of ligand","NOD1/2 Signaling Pathway"),
                   condition= "PTB",
                   control= "Cont_PTB",
                   cell_subsets=c("HBC_1","HBC_4_IL1B_MMP19","HBC_6","HBC_11","HBC_2_CXCL2_CXCL3","HBC_10_cycling","mixed_HBC_13_S100A8","mMac_upper","CD16_monocyte","CD14_monocyte","cDC"),
                   width=4, height=8,
                   pdf_title="Stress2",
                   cluster_rows = TRUE)

plot_heatmap_genes(PW_list=c("Cellular response to heat stress"),
                   condition= "PTB",
                   control= "Cont_PTB",
                   cell_subsets=c("HBC_1","HBC_4_IL1B_MMP19","HBC_6","mMac_upper","CD16_monocyte","CD14_monocyte","cDC","VCT_1","VCTfusion_1","VCTp_1"),
                   width=4, height=8,
                   pdf_title="Stress3",
                   cluster_rows = TRUE)

plot_heatmap_genes(PW_list=c("Activation of lymphatic system cells"),
                   condition= "PTB",
                   control= "Cont_PTB",
                   cell_subsets=c("HBC_4_IL1B_MMP19","HBC_6","HBC_2_CXCL2_CXCL3"),
                   width=4, height=8,
                   pdf_title="",
                   cluster_rows = TRUE)

plot_heatmap_genes(PW_list=c("Cell movement"),
                   condition= "PTB",
                   control= "Cont_PTB",
                   cell_subsets=c("HBC_4_IL1B_MMP19","HBC_6","HBC_1","mMac_upper","VCT_1","VCTfusion_1","VCTp_1"), #,"fPV2","fPV3"),
                   width=4, height=17,
                   pdf_title="",
                   cluster_rows = TRUE)
plot_heatmap_genes(PW_list=c("Cell movement"),
                   condition= "PTB",
                   control= "Cont_PTB",
                   cell_subsets=c("HBC_4_IL1B_MMP19","HBC_1","mMac_upper","VCT_1","VCTfusion_1","VCTp_1","Endo_m1","fPV2","fPV3"),
                   width=4, height=17,
                   pdf_title="Cell movement Sig",
                   cluster_rows = TRUE,
                    onlySig=TRUE)


plot_heatmap_genes(PW_list=c("Migration of cells"),
                   condition= "PTB",
                   control= "Cont_PTB",
                   cell_subsets=c("HBC_4_IL1B_MMP19","HBC_1","mMac_upper","VCT_1","VCTfusion_1","VCTp_1","Endo_m1","fPV2","fPV3"),
                   width=4, height=17,
                   pdf_title="",
                   cluster_rows = TRUE,
                   onlySig=TRUE)

plot_heatmap_genes(PW_list=c("Cell viability"),
                   condition= "PTB",
                   control= "Cont_PTB",
                   cell_subsets=c("HBC_4_IL1B_MMP19","HBC_1","mMac_upper","VCT_1","VCTfusion_1","VCTp_1","Endo_m1","fPV2","fPV3"),
                   width=4, height=17,
                   pdf_title="",
                   cluster_rows = TRUE,
                   onlySig=TRUE)



###############################################################
#                          PE
###############################################################
source("/projects/Pregnancy_placenta/Code_R/Figures/figure_functions.R")
setwd("/projects/Pregnancy_placenta/Results/20_Figures/Pathways/PE_gene_dotplot/Full_Gene_List_heatmaps/")


plot_heatmap_genes(PW_list=c("Adhesion of myeloid cells"),
                   condition= "TermPE",
                   control= "Cont_unlabored",
                   cell_subsets=c("HBC_3","mixed_HBC_13_S100A8","CD14_monocyte"),
                   width=4, height=5,
                   pdf_title="",
                   cluster_rows = TRUE,
                   onlySig=TRUE)

plot_heatmap_genes(PW_list=c("Adhesion of immune cells"),
                   condition= "TermPE",
                   control= "Cont_unlabored",
                   cell_subsets=c("HBC_3","mixed_HBC_13_S100A8","CD14_monocyte"),
                   width=4, height=5,
                   pdf_title="",
                   cluster_rows = TRUE,
                   onlySig=FALSE)
plot_heatmap_genes(PW_list=c("Binding of professional phagocytic cells"),
                   condition= "TermPE",
                   control= "Cont_unlabored",
                   cell_subsets=c("HBC_3","mixed_HBC_13_S100A8","HBC_2_CXCL2_CXCL3"),
                   width=4, height=5,
                   pdf_title="",
                   cluster_rows = TRUE,
                   onlySig=FALSE)

plot_heatmap_genes(PW_list=c("Migration of vascular endothelial cells"),
                   condition= "TermPE",
                   control= "Cont_unlabored",
                   cell_subsets=c("HBC_3","HBC_2_CXCL2_CXCL3"),
                   width=6, height=5,
                   pdf_title="",
                   cluster_rows = TRUE,
                   onlySig=TRUE)
plot_heatmap_genes(PW_list=c("Engulfment of phagocytes"),
                   condition= "TermPE",
                   control= "Cont_unlabored",
                   cell_subsets=c("HBC_3","HBC_2_CXCL2_CXCL3","HBC_4_IL1B_MMP19"),
                   width=5, height=5,
                   pdf_title="",
                   cluster_rows = TRUE,
                   onlySig=TRUE)
plot_heatmap_genes(PW_list=c("Innate immune response"),
                   condition= "TermPE",
                   control= "Cont_unlabored",
                   cell_subsets=c("HBC_3","HBC_2_CXCL2_CXCL3","HBC_4_IL1B_MMP19"),
                   width=5, height=5,
                   pdf_title="",
                   cluster_rows = TRUE,
                   onlySig=TRUE)
plot_heatmap_genes(PW_list=c("Inflammatory response"),
                   condition= "TermPE",
                   control= "Cont_unlabored",
                   cell_subsets=c("HBC_3","HBC_2_CXCL2_CXCL3","HBC_4_IL1B_MMP19"),
                   width=5, height=7,
                   pdf_title="",
                   cluster_rows = TRUE,
                   onlySig=TRUE)
plot_heatmap_genes(PW_list=c("Inflammatory response"),
                   condition= "PrePE",
                   control= "Cont_unlabored",
                   cell_subsets=c("HBC_3","HBC_2_CXCL2_CXCL3"),
                   width=5, height=7,
                   pdf_title="",
                   cluster_rows = TRUE,
                   onlySig=TRUE)

plot_heatmap_genes(PW_list=c("Interferon alpha/beta signaling"),
                   condition= "TermPE",
                   control= "Cont_unlabored",
                   cell_subsets=c("HBC_3","HBC_2_CXCL2_CXCL3","HBC_4_IL1B_MMP19"),
                   width=5, height=7,
                   pdf_title="",
                   cluster_rows = TRUE,
                   onlySig=TRUE)
plot_heatmap_genes(PW_list=c("ISG15 antiviral mechanism"),
                   condition= "TermPE",
                   control= "Cont_unlabored",
                   cell_subsets=c("HBC_3","HBC_2_CXCL2_CXCL3","HBC_4_IL1B_MMP19"),
                   width=5, height=7,
                   pdf_title="",
                   cluster_rows = TRUE,
                   onlySig=TRUE)
plot_heatmap_genes(PW_list=c("Leukocyte migration"),
                   condition= "TermPE",
                   control= "Cont_unlabored",
                   cell_subsets=c("HBC_3","HBC_2_CXCL2_CXCL3","mixed_HBC_13_S100A8"),
                   width=5, height=7,
                   pdf_title="",
                   cluster_rows = TRUE,
                   onlySig=TRUE)
plot_heatmap_genes(PW_list=c("Cell movement of phagocytes"),
                   condition= "TermPE",
                   control= "Cont_unlabored",
                   cell_subsets=c("HBC_3","HBC_2_CXCL2_CXCL3","mixed_HBC_13_S100A8"),
                   width=5, height=7,
                   pdf_title="",
                   cluster_rows = TRUE,
                   onlySig=TRUE)
plot_heatmap_genes(PW_list=c("Synthesis of reactive oxygen species"),
                   condition= "TermPE",
                   control= "Cont_unlabored",
                   cell_subsets=c("HBC_3","HBC_2_CXCL2_CXCL3","mixed_HBC_13_S100A8"),
                   width=5, height=7,
                   pdf_title="",
                   cluster_rows = TRUE,
                   onlySig=FALSE)
plot_heatmap_genes(PW_list=c("Activation of leukocytes"),
                   condition= "TermPE",
                   control= "Cont_unlabored",
                   cell_subsets=c("HBC_3","HBC_2_CXCL2_CXCL3","mixed_HBC_13_S100A8","HBC_4_IL1B_MMP19","CD14_monocyte"),
                   width=5, height=7,
                   pdf_title="",
                   cluster_rows = TRUE,
                   onlySig=FALSE)
plot_heatmap_genes(PW_list=c("Angiogenesis"),
                   condition= "TermPE",
                   control= "Cont_unlabored",
                   cell_subsets=c("HBC_3","HBC_2_CXCL2_CXCL3","mixed_HBC_13_S100A8"),
                   width=5, height=7,
                   pdf_title="",
                   cluster_rows = TRUE,
                   onlySig=FALSE)
plot_heatmap_genes(PW_list=c("Leukopoiesis"),
                   condition= "TermPE",
                   control= "Cont_unlabored",
                   cell_subsets=c("CD14_monocyte"),
                   width=5, height=7,
                   pdf_title="",
                   cluster_rows = TRUE,
                   onlySig=FALSE)
plot_heatmap_genes(PW_list=c("Leukopoiesis"),
                   condition= "TermPE",
                   control= "Cont_unlabored",
                   cell_subsets=c("CD14_monocyte"),
                   width=5, height=7,
                   pdf_title="",
                   cluster_rows = TRUE,
                   onlySig=FALSE)
plot_heatmap_genes(PW_list=c("Neutrophil degranulation"),
                   condition= "PrePE",
                   control= "Cont_unlabored",
                   cell_subsets=c("CD14_monocyte","SCT"),
                   width=5, height=7,
                   pdf_title="Imm_response_SCT",
                   cluster_rows = TRUE,
                   onlySig=FALSE)

setwd("/projects/Pregnancy_placenta/Results/20_Figures/Pathways/PE_gene_dotplot/Full_Gene_List_heatmaps/Trophoblasts/")
source("/projects/Pregnancy_placenta/Code_R/Figures/figure_functions.R")

plot_heatmap_genes(PW_list=c("Inflammatory response"),
                   condition= "TermPE",
                   control= "Cont_unlabored",
                   cell_subsets=c("SCT","VCT_1","VCT_3","VCTfusion_1"),
                   width=5, height=7,
                   pdf_title="",
                   cluster_rows = TRUE,
                   onlySig=FALSE)

plot_heatmap_genes(PW_list=c("Viral Infection"),
                   condition= "TermPE",
                   control= "Cont_unlabored",
                   cell_subsets=c("SCT","VCT_1","VCT_3","VCTfusion_1"),
                   width=5, height=7,
                   pdf_title="",
                   cluster_rows = TRUE,
                   onlySig=TRUE)
plot_heatmap_genes(PW_list=c("Leukopoiesis"),
                   condition= "TermPE",
                   control= "Cont_unlabored",
                   cell_subsets=c("SCT","VCT_1","VCT_3","VCTfusion_1"),
                   width=5, height=7,
                   pdf_title="",
                   cluster_rows = TRUE,
                   onlySig=TRUE)
plot_heatmap_genes(PW_list=c("Leukopoiesis"),
                   condition= "PrePE",
                   control= "Cont_unlabored",
                   cell_subsets=c("SCT"),
                   width=5, height=7,
                   pdf_title="",
                   cluster_rows = TRUE,
                   onlySig=TRUE)

plot_heatmap_genes(PW_list=c("Activation of leukocytes"),
                   condition= "PrePE",
                   control= "Cont_unlabored",
                   cell_subsets=c("SCT"),
                   width=5, height=7,
                   pdf_title="",
                   cluster_rows = TRUE,
                   onlySig=TRUE)
plot_heatmap_genes(PW_list=c("Leukopoiesis"),
                   condition= "PrePE",
                   control= "Cont_unlabored",
                   cell_subsets=c("SCT"),
                   width=5, height=7,
                   pdf_title="",
                   cluster_rows = TRUE,
                   onlySig=TRUE)
plot_heatmap_genes(PW_list=c("Leukocyte migration"),
                   condition= "PrePE",
                   control= "Cont_unlabored",
                   cell_subsets=c("SCT"),
                   width=5, height=7,
                   pdf_title="",
                   cluster_rows = TRUE,
                   onlySig=TRUE)
plot_heatmap_genes(PW_list=c("Warburg Effect Signaling Pathway"),
                   condition= "PrePE",
                   control= "Cont_unlabored",
                   cell_subsets=c("SCT"),
                   width=3, height=4,
                   pdf_title="",
                   cluster_rows = TRUE,
                   onlySig=TRUE)
plot_heatmap_genes(PW_list=c("HIF1α Signaling"),
                   condition= "PrePE",
                   control= "Cont_unlabored",
                   cell_subsets=c("SCT"),
                   width=3, height=4,
                   pdf_title="",
                   cluster_rows = TRUE,
                   onlySig=TRUE)



plot_heatmap_genes(PW_list=c("Fatty acid metabolism"),
                   condition= "TermPE",
                   control= "Cont_unlabored",
                   cell_subsets=c("VCT_1","VCTfusion_1","VCTfusion_2"),
                   width=5, height=7,
                   pdf_title="",
                   cluster_rows = TRUE,
                   onlySig=TRUE)
plot_heatmap_genes(PW_list=c("Fatty acid metabolism"),
                   condition= "PrePE",
                   control= "Cont_unlabored",
                   cell_subsets=c("VCT_1","VCTfusion_1","VCTfusion_2"),
                   width=5, height=7,
                   pdf_title="",
                   cluster_rows = TRUE,
                   onlySig=TRUE)

plot_heatmap_genes(PW_list=c("Angiogenesis"),
                   condition= "TermPE",
                   control= "Cont_unlabored",
                   cell_subsets=c("SCT","VCT_1","VCT_3","VCTfusion_1"),
                   width=5, height=7,
                   pdf_title="",
                   cluster_rows = TRUE,
                   onlySig=TRUE)
plot_heatmap_genes(PW_list=c("Angiogenesis"),
                   condition= "PrePE",
                   control= "Cont_unlabored",
                   cell_subsets=c("SCT","VCT_1","VCTfusion_1","VCTfusion_2","VCTccc"),
                   width=5, height=7,
                   pdf_title="",
                   cluster_rows = TRUE,
                   onlySig=TRUE)



setwd("/projects/Pregnancy_placenta/Results/20_Figures/Pathways/PE_gene_dotplot/Full_Gene_List_heatmaps/Endo//")
source("/projects/Pregnancy_placenta/Code_R/Figures/figure_functions.R")




plot_heatmap_genes(PW_list=c("Angiogenesis","Migration of vascular endothelial cells"),
                   condition= "PrePE",
                   control= "Cont_unlabored",
                   cell_subsets=c("Endo_f1","Endo_f2","Endo_f3"),
                   width=5, height=12,
                   pdf_title="Angio_related_Endo",
                   cluster_rows = TRUE,
                   onlySig=TRUE)

plot_heatmap_genes(PW_list=c('Growth of blood vessel'),#,"Migration of vascular endothelial cells"),
                   condition= "PrePE",
                   control= "Cont_unlabored",
                   cell_subsets=c("Endo_f1","Endo_f2","Endo_f3"),
                   width=5, height=7,
                   pdf_title="",#"Angio_related_Endo",
                   cluster_rows = TRUE,
                   onlySig=TRUE)

plot_heatmap_genes(PW_list=c('Leukopoiesis'),
                   condition= "PrePE",
                   control= "Cont_unlabored",
                   cell_subsets=c("Endo_f1","Endo_f2","Endo_f3"),
                   width=5, height=7,
                   pdf_title="",#"Angio_related_Endo",
                   cluster_rows = TRUE,
                   onlySig=TRUE)

plot_heatmap_genes(PW_list=c('Neutrophil degranulation'),
                   condition= "PrePE",
                   control= "Cont_unlabored",
                   cell_subsets=c("Endo_f1","Endo_f3"),
                   width=5, height=7,
                   pdf_title="",#"Angio_related_Endo",
                   cluster_rows = TRUE,
                   onlySig=TRUE)

plot_heatmap_genes(PW_list=c('Viral Infection'),
                   condition= "PrePE",
                   control= "Cont_unlabored",
                   cell_subsets=c("Endo_f1","Endo_f2","Endo_f3"),
                   width=5, height=15,
                   pdf_title="",#"Angio_related_Endo",
                   cluster_rows = TRUE,
                   onlySig=TRUE)
plot_heatmap_genes(PW_list=c('Inflammatory response'),
                   condition= "PrePE",
                   control= "Cont_unlabored",
                   cell_subsets=c("Endo_f1","Endo_f2","Endo_f3"),
                   width=5, height=10,
                   pdf_title="",#"Angio_related_Endo",
                   cluster_rows = TRUE,
                   onlySig=TRUE)
plot_heatmap_genes(PW_list=c('Permeability of vasculature'),
                   condition= "PrePE",
                   control= "Cont_unlabored",
                   cell_subsets=c("Endo_f1","Endo_f2","Endo_f3"),
                   width=5, height=7,
                   pdf_title="",#"Angio_related_Endo",
                   cluster_rows = TRUE,
                   onlySig=TRUE)



