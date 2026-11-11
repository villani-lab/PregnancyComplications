suppressPackageStartupMessages({
  library(openxlsx)
  library(readxl)
  library(tidyr)
  library(stringr)
  library(ComplexHeatmap)
  library(colorRamp2)
  library(dplyr)
  library(data.table)})



# This functions shows the expression of genes that are part of an IPa pathway that came up as significant.
# It will not work if the pathway is not significant for the specific condition and cell type that are given, so double check that if it fails.
# NOTE: this version reads the IPA export files from disk itself (<cell_type>-<condition>_vs_<control>-diseases.xls,
# falling back to the "-PW.xls" canonical-pathway file, in loc_of_IPA_gene_files). A different function with a
# similar purpose that works on an already-loaded IPA data frame is get_log2fc_of_pathway_genes_from_df()
# in funcs_gene_pathway.R; the two are intentionally kept separate and have different arguments.
get_log2fc_of_pathway_genes_from_IPA_files = function(pathway, #pathway of interest
                                       cell_type,
                                       condition,
                                       control,
                                       loc_of_IPA_gene_files,
                                       df_loc) # contains the info as to where the FC files are
{
  # Read pathway genes from "diseases and biofunctions" file
  file_name = paste0(cell_type,"-",condition,"_vs_",control,"-diseases.xls")
  df = as.data.frame(read_excel(col_names = TRUE,skip = 1,path = paste0(loc_of_IPA_gene_files,"/",file_name)))
  df_path_int = df[df[,'Diseases or Functions Annotation']== pathway,]
  if (dim(df_path_int)[1]==0)
  { # Read pathway genes from "Canonical pathways" file
    file_name = paste0(cell_type,"-",condition,"_vs_",control,"-PW.xls")
    df = as.data.frame(read_excel(col_names = TRUE,skip = 1,path = paste0(loc_of_IPA_gene_files,"/",file_name)))
    df_path_int = df[df[,'Ingenuity Canonical Pathways']== pathway,]
    if (dim(df_path_int)[1]==0)
    {
      print(paste0(pathway," is not present in either diseases or PW files"))
      return()
    }
  }

  # get the DEGs of the comparison of interest of the cell
  mask = (df_loc$Cell_type==cell_type & df_loc$name == paste0(condition,"_vs_",control))
  loc_int = df_loc[mask,]
  info_on_loc = str_split(loc_int$name_ipa,"_",4)[[1]]
  parent_or_sub = info_on_loc[[1]][1]
  if (parent_or_sub=="Parent")
      parent_or_sub = "Parentlin"
  csv_file = paste0("/projects/Pregnancy_placenta/Results/13_DE_PAL/output_csv_IPA/redo_", parent_or_sub,"/",loc_int$name_ipa,".csv")
  log2fc_gene_df = read.csv(csv_file)
  vals_int = log2fc_gene_df[,grepl(paste0(cell_type,".",condition,"_vs_",control),colnames(log2fc_gene_df))]
  rownames(vals_int) = log2fc_gene_df$GeneName

  # only keep the genes that are in the pathway
  vals_int = vals_int[which(rownames(vals_int) %in% str_split(df_path_int$Molecules,",")[[1]]),]
  vals_int = data.frame(vals_int)
  colnames(vals_int) = c("log2fc","padj")
  vals_int$ast = stars.pval(1*vals_int$padj)
  vals_int$gene_name = rownames(vals_int)
  rownames(vals_int) = NULL
  return(vals_int)
}


stars.pval <- function(x){
  stars <- c("***", "**", "*", "n.s.")
  cutoffs <- c(0, 0.001, 0.01, 0.05, 1)
  i <- findInterval(x, cutoffs, left.open = TRUE, rightmost.closed = TRUE)
  stars[i]
}



plot_heatmap_genes = function(PW_list, condition, control, cell_subsets, width=5, height=10, suffix="", cluster_rows = FALSE, pdf_title="", onlySig=FALSE)
{
  if (length(PW_list)==1)
    PW_list = c(PW_list)
  # Reading FC files
  parent_lin = read.csv("/projects/Pregnancy_placenta/Results/13_DE_PAL/output_csv_IPA/redo_Parentlin/Parent_Test_3.csv")
  sub_lin = read.csv("/projects/Pregnancy_placenta/Results/13_DE_PAL/output_csv_IPA/redo_sublin/sublin_Test_3.csv")
  all_lin_loc = rbind(data.frame(parent_lin),data.frame(sub_lin))

  loc_of_IPA_gene_files = "/projects/Pregnancy_placenta/Results/20_Figures/Pathways/IPA_genes/"
  merged_pathways = data.frame()
  for (PW in PW_list)
  {
    print(PW)
    for(cell_type in cell_subsets)
    {
      print(cell_type)
      curr_df = get_log2fc_of_pathway_genes_from_IPA_files(pathway = PW,
                                  cell_type = cell_type,
                                  condition = condition,
                                  control = control,
                                  loc_of_IPA_gene_files = loc_of_IPA_gene_files,
                                  df_loc = all_lin_loc)

      curr_df$condition = cell_type
      merged_pathways = rbind(merged_pathways, curr_df)
    }
  }

  # Some genes may repeat in different pathways - lets keep one row per gene
  merged_pathways = unique(merged_pathways)
  # If filtration is needed - filter
  if (onlySig)
    merged_pathways = merged_pathways[merged_pathways["ast"]!="n.s.",]

  # make heatmap
  log2fc_df = as.data.frame(merged_pathways%>%
                             pivot_wider(id_cols = 'gene_name',
                                         names_from = 'condition',
                                         values_from = "log2fc",values_fill = NA))
  ast_df = as.data.frame(merged_pathways %>%
                             pivot_wider(id_cols = 'gene_name',
                                         names_from = 'condition',
                                         values_from = "ast",values_fill = NA))

  log2fc_mat = as.matrix(log2fc_df[,-1])
  rownames(log2fc_mat) = log2fc_df$gene_name
  log2fc_mat[is.na(log2fc_mat)]=0

  ast_mat = as.matrix(ast_df[,-1])
  rownames(ast_mat) = ast_df$gene_name
  ast_mat[is.na(ast_mat)]="n.s."

  minv <- min(0, min(log2fc_mat, na.rm = TRUE))
  maxv <- max(0, max(log2fc_mat, na.rm = TRUE))
  if (minv < 0 && maxv > 0) {
    # both negative and positive: blue -> white -> red, with white at 0
    col_fun <- colorRamp2(c(minv, 0, maxv), c("blue", "white", "red"))
  } else if (minv >= 0) {
    # only positive: white -> red
    col_fun <- colorRamp2(c(0, maxv), c("white", "red"))
  } else {
    # only negative: blue -> white
    col_fun <- colorRamp2(c(minv, 0), c("blue", "white"))
  }

  print("Plotting")
  if (pdf_title=="")
  {
    pdf_title = paste(PW_list, collapse = "-")
    pdf_title = gsub("/","-", pdf_title)
  }
  pdf(paste0(pdf_title,"_",condition,"_vs_",control,suffix,".pdf"), width=width, height=height)
  print(Heatmap(log2fc_mat,
                col= col_fun,
                name = "log2FC",
                show_row_names = TRUE,
                row_names_side = "left",
                row_names_gp = gpar(fontsize = 8),
                column_title = paste0(pdf_title," - ",condition," vs ",control),
                cluster_columns = FALSE,
                cluster_rows = cluster_rows,
                cell_fun = function(j, i, x, y, w, h, fill) {
                 if (!is.na(ast_mat[i, j]) && ast_mat[i, j] != "n.s.") {
                   grid::grid.text(ast_mat[i, j], x, y - grid::unit(0.6, "mm"), just = c("center", "center"))
                 }}))
  dev.off()
}





# Plots a heatmap of PWs (y axis) for various comparisons (<subset>-<condition>_<control>, x axis)
# Color denotes z-scores, p-values appear as asterkis in the heatmap
# If orderColBy=="Cells" than the columns will be sorted by Lineage. Columns are currently clustered so the cell subsets are not sorted, though this can be changed in the Heatmap function by not allowing column clustering.
# withNAs - should pathwyas with "NA" as an annotaiton be removed from the data or not
# zscore_min_cp, padj_thresh_cp - z-score and p-val thresholds for canonical pathways (as appreas in df_comp_merged["Type"]
# zscore_min_db, padj_thresh_db - z-score and p-val thresholds for diseases and biofunctions (as appreas in df_comp_merged["Type"]
plot_heatmap_PW <- function(df_comp_merged,
                                    padj_or_pval_ast = 'padj',
                                    comp_int = NULL,
                                    zscore_min_cp = 0,
                                    padj_thresh_cp = NULL,
                                    zscore_min_db = 0,
                                    padj_thresh_db = NULL,
                                    save_pdf = TRUE,
                                    output_prefix = "",
                                    output_dir = ".",
                                    orderColBy = "Cells", # How to order the columns in - by "Cells" will order it by lineage, put something else if you want to use orcer_vector to order the columns
                                    order_vector=c(),  # vector of column names to order the heatmap by
                                    pdf_width = 9,
                                    pdf_height = 6,
                                    cluster_names_col = "Clust_Comparison",
                                    split_cols = c(PrePE="Grey", TermPE="black"),
                                    filter_more="Y",
                                    plot_title="",
                                    pw_annot_file=NULL,
                                    row_display_mode = "pathway_name", # or pathway_annotation_group
                                    font_family = "Arial",
                                    split_by_cond_OR_ctrl = "Cond",
                                    withNAs = FALSE)
  {
  df_comp_merged = na.omit(df_comp_merged)  # omiting rows with z-score==NA
  df_comp_merged$pval_ast = stars.pval(10^(-1*df_comp_merged$log10_pval))
  df_comp_merged$qval_ast = stars.pval(10^(-1*df_comp_merged$log10_padj))

  # how many column splits do we have within lienage...
  cond_number = length(unique(df_comp_merged[[cluster_names_col]])) # we need to know the number of conditions we have to change ploting
  split_group_col = "Condition"
  if (cond_number>1) {
    if(split_by_cond_OR_ctrl=="Cond"){
      df_comp_merged["Cluster_cond"] = paste(df_comp_merged[,"Cluster"], df_comp_merged[,"Condition"], sep="_")
      cluster_names_col = "Cluster_cond"
      split_group_col = "Condition"
    }
    if(split_by_cond_OR_ctrl=="Ctrl"){
      df_comp_merged["Cluster_cond"] = paste(df_comp_merged[,"Cluster"], df_comp_merged[,"Control"], sep="_")
      cluster_names_col = "Cluster_cond"
      split_group_col = "Control"
    }
  }


  # Filtering
  if(filter_more =="N") # this filtering will *keep* any cell in the row if there is at least one cell subset that is significant for this pathway
  {
      filtered_df <- df_comp_merged  %>%
        group_by(Canonical_Pathways) %>%
        filter(any(((log10_padj >= -log10(padj_thresh_cp)) & (abs(Zscore) > zscore_min_cp) & (Type == "Canonical_Pathways")) |
                   ((log10_padj >= -log10(padj_thresh_db)) & (abs(Zscore) > zscore_min_db) & (Type == "Diseases_and_biofunctions"))))
  }
  if(filter_more =="Y") # this filtering will filter out any cell in the matrix that is lower than the thresholds.
  {
    filtered_df<- df_comp_merged %>% dplyr::filter(abs(Zscore) > zscore_min)
    filtered_df <- filtered_df %>% dplyr::filter(log10_padj >= -log10(padj_thresh))
  }

  dup_keys <- filtered_df %>%
    dplyr::count(Canonical_Pathways, .data[[cluster_names_col]], name = "n") %>%
    dplyr::filter(n > 1)
  if (nrow(dup_keys) > 0)
  {
    message("Collapsing ", nrow(dup_keys), " duplicated pathway/column combinations by keeping the row with the largest absolute Z-score.")
    filtered_df <- filtered_df %>%
      dplyr::group_by(Canonical_Pathways, .data[[cluster_names_col]]) %>%
      dplyr::arrange(dplyr::desc(abs(Zscore)), dplyr::desc(log10_padj), .by_group = TRUE) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup()
  }

  # Prepare Z-score matrix for ploting
  zscore_df = as.data.frame(filtered_df %>%
                            pivot_wider(id_cols = 'Canonical_Pathways',
                                        names_from = cluster_names_col,
                                        values_from = "Zscore",values_fill = NA))

  # Prepare P-value matrix for ploting
  if(padj_or_pval_ast=="padj")
  {
    ast_df = as.data.frame(filtered_df %>%
                           pivot_wider(id_cols = 'Canonical_Pathways',
                                       names_from = cluster_names_col,
                                       values_from = "qval_ast",values_fill = NA))
  }
  else
  {
    ast_df = as.data.frame(filtered_df %>%
                             pivot_wider(id_cols = 'Canonical_Pathways',
                                         names_from = cluster_names_col,
                                         values_from = "pval_ast",values_fill = NA))
  }

  zscore_mat = as.matrix(zscore_df[,-1])
  rownames(zscore_mat) = zscore_df$Canonical_Pathways

  ast_mat = as.matrix(ast_df[,-1])
  rownames(ast_mat) = ast_df$Canonical_Pathways

  if (orderColBy == "Cells") # Column order by cell type or given vector <order_vector>
  {
    # conversion cell -> lineage
    if (cond_number>1)
    {
      conv_cell_lin = data.frame(unique(filtered_df[,c("Cluster","Lineage",split_group_col,cluster_names_col)]))
      if (length(order_vector)>0)
      {
        print(" ** sorting ** ")
        missing <- setdiff(order_vector, colnames(zscore_mat))
        if (length(missing) > 0) stop("order_vector has unknown columns: ", paste(missing, collapse=", "))
        rownames(conv_cell_lin) = paste(conv_cell_lin$Cluster, conv_cell_lin[[split_group_col]], sep="_")
        conv_cell_lin = conv_cell_lin[order_vector,]

        zscore_mat = zscore_mat[,order_vector]
        ast_mat = ast_mat[,order_vector]
      }
      else
      {
        conv_cell_lin = conv_cell_lin[order(conv_cell_lin[[split_group_col]], conv_cell_lin$Lineage, conv_cell_lin$Cluster),]
        col_order <- conv_cell_lin[,cluster_names_col]
        zscore_mat = zscore_mat[,col_order]
        ast_mat = ast_mat[,col_order]
      }
      conv_cell_lin$Split_lin = paste0(conv_cell_lin[[split_group_col]], "_", conv_cell_lin$Lineage)
      column_split_cols = conv_cell_lin$Split_lin
      # --- prevent column_split from reordering slices alphabetically ---
      column_split_cols <- factor(column_split_cols, levels = unique(column_split_cols))
    }
    else
    {
        conv_cell_lin = unique(filtered_df[,c("Cluster","Lineage")])
        if (length(order_vector)!=0)
        {
          print("sorting by order_vector")
          conv_cell_lin <- conv_cell_lin[order(match(conv_cell_lin$Cluster, order_vector)), ]
        }
        else
          conv_cell_lin = conv_cell_lin[order(conv_cell_lin$Lineage, conv_cell_lin$Cluster),]
        col_order <- conv_cell_lin$Cluster
        zscore_mat = zscore_mat[,col_order]
        ast_mat = ast_mat[,col_order]
        column_split_cols = conv_cell_lin$Lineage
        column_split_cols <- factor(column_split_cols, levels = unique(column_split_cols))
    }
  }
  else # order by control type
  {
      control_celltype = sapply(colnames(zscore_mat), function(x) paste({strsplit(x, "_vs_|-")[[1]][c(3,1)]},collapse="-"))
      col_order = colnames(zscore_mat)[order(control_celltype)]
      zscore_mat = zscore_mat[,col_order]
      ast_mat = ast_mat[,col_order]
      column_split_cols = sapply(colnames(zscore_mat), function(x) {strsplit(x, "_vs_")[[1]][2]})
  }

  # Generate column annotations
  if ((orderColBy=="Cells") | (length(order_vector)>0)) # We'll annotated lineages
  {
    if (cond_number>1)
    {
      lin_colors = c(Myeloid = "#590902", HBC = "#BF3604", Trophoblasts = "#60B364", Stromal = "#022873", Endothelial = "#6B0AC9") #fPV = "#056CF2",fFB = "#0597F2", T = "#F27405", NK = "#F2AE2E", B = "#F2668B",    EVT = "#638c3e",
      #lin_colors = c(Myeloid = "#BF3604", HBC = "#BF3604", Trophoblasts = "#60B364", Stromal = "#022873", Endothelial = "#6B0AC9") #fPV = "#056CF2",fFB = "#0597F2", T = "#F27405", NK = "#F2AE2E", B = "#F2668B",    EVT = "#638c3e",
      cond_colors = split_cols
      lin_colors = lin_colors[unique(conv_cell_lin$Lineage)]
      split_values = unique(conv_cell_lin[[split_group_col]])
      cond_colors = cond_colors[split_values]
      missing_split_values = split_values[is.na(cond_colors) | is.na(names(cond_colors))]
      if (length(missing_split_values) > 0)
      {
        extra_cols = grDevices::hcl.colors(length(missing_split_values), palette = "Set 2")
        names(extra_cols) = missing_split_values
        cond_colors = c(cond_colors[!(is.na(cond_colors) | is.na(names(cond_colors)))], extra_cols)
      }
      cond_colors = cond_colors[split_values]
      annot_df = data.frame(split_value = conv_cell_lin[[split_group_col]],
                            Lineage = conv_cell_lin$Lineage,
                            check.names = FALSE)
      colnames(annot_df)[1] = split_group_col
      annot_colors = list(Lineage = lin_colors)
      annot_colors[[split_group_col]] = cond_colors
      annot_legend_param = list(Lineage = list(nrow = 1,
                                               title_position = "topcenter",
                                               direction = "horizontal",
                                               title_gp = grid::gpar(fontfamily = font_family, fontface = "plain"),
                                               labels_gp = grid::gpar(fontfamily = font_family, fontface = "plain")))
      annot_legend_param[[split_group_col]] = list(title_gp = grid::gpar(fontfamily = font_family, fontface = "plain"),
                                                   labels_gp = grid::gpar(fontfamily = font_family, fontface = "plain"))
      Col_annotation <- ComplexHeatmap::HeatmapAnnotation(df = annot_df,
                                                          col = annot_colors,
                                                          annotation_legend_param = annot_legend_param,
                                                          annotation_name_gp = grid::gpar(fontfamily = font_family, fontface = "plain"),
                                                          show_annotation_name = FALSE,
                                                          border = TRUE)
    }
    else
    {
      #lin_colors = c(Myeloid = "#590902", HBC = "#BF3604", Trophoblasts = "#60B364", Stromal = "#022873", Endothelial = "#6B0AC9") #fPV = "#056CF2",fFB = "#0597F2", T = "#F27405", NK = "#F2AE2E", B = "#F2668B",    EVT = "#638c3e",
      lin_colors = c(Myeloid = "#BF3604", HBC = "#BF3604", Trophoblasts = "#60B364", Stromal = "#022873", Endothelial = "#6B0AC9") #fPV = "#056CF2",fFB = "#0597F2", T = "#F27405", NK = "#F2AE2E", B = "#F2668B",    EVT = "#638c3e",
      lin_colors = lin_colors[unique(conv_cell_lin$Lineage)]
      Col_annotation <- ComplexHeatmap::HeatmapAnnotation(Lineage = conv_cell_lin$Lineage,
                                                         col = list(Lineage = lin_colors),
                                                         annotation_legend_param = list(Lineage = list(nrow = 1,
                                                                                                       title_position = "topcenter",
                                                                                                       direction = "horizontal",
                                                                                                       title_gp = grid::gpar(fontfamily = font_family, fontface = "plain"),
                                                                                                       labels_gp = grid::gpar(fontfamily = font_family, fontface = "plain"))),
                                                         annotation_name_gp = grid::gpar(fontfamily = font_family, fontface = "plain"),
                                                         show_annotation_name = FALSE,
                                                         border = TRUE)

    }
  }
  else
  {
    mat_column_annotation_split = stringr::str_split_fixed(colnames(zscore_mat), "_vs_", 2)[, 2]
    cond2_blocks <- c("Cont_labored" = "#0073C2FF",
                    "Cont_unlabored" =  "#EFC000FF",
                    "Cont_PTB" = "#CD534CFF")
    Col_annotation <- ComplexHeatmap::HeatmapAnnotation(Ctrl = mat_column_annotation_split,
                                                       col = list(Ctrl = cond2_blocks),
                                                       annotation_legend_param = list(
                                                         Ctrl = list(nrow = 1,
                                                         title_position = "topcenter",
                                                         direction = "horizontal",
                                                         title_gp = grid::gpar(fontfamily = font_family, fontface = "plain"),
                                                         labels_gp = grid::gpar(fontfamily = font_family, fontface = "plain"))),
                                                       annotation_name_gp = grid::gpar(fontfamily = font_family, fontface = "plain"),
                                                       show_annotation_name = FALSE,
                                                       border = TRUE)
  }


  if (!is.null(pw_annot_file)) # row annotations and row reordering - grouping by pathways annotations
  {
    #pw_annot_col = c(, "#FF1493", "#3F007D", "#7FFF00", "#7F0000", "#FFD700", "black", "#009E73", "#56B4E9", "#F0E442", "#CC79A7")
    #names(pw_annot_col) = c("Hormones","Inflammatory","Stress","Immune response","Pathogen/danger sensing","Activation","Cell death","Metabolism","Repair","Viability","Inhibition","ECM","Contraction","Collagen-related","?","Cell viability", "Invasiveness", "Angiogenesis","Migration", "Cytoskeleton organization","Proliferation", "Neuro-related", "Adhesion", "Fibrosis","Engulfment", "Neoplasia", "Infection", "Damage")
    pw_annot_col = c("#b15928","#fb9a99","#e31a1c","#f16913","#feb24c","#a6cee3","#1f78b4","#b2df8a","#33a02c","#cab2d6","#6a3d9a","#ffff99","#7fcdbb","#bc80bd","#d9d9d9","#737373", "#00BCD4", "#005F73","#000000")
    names(pw_annot_col) = c("Cell death and survival","Migration and Invasion","Cytoskeleton","Tissue Repair/ECM","Cell membrane","Metabolic/ Hormone","Angiogenesis","Inflammation/ Immune Activation","Stress","Collagen","Autophagy","Mitochondria","Metabolism","Hormones","Neuro-related","Engulfment", "Adhesion","Angiogenesis/Vascular","Not_Annotated")

    pw_annot = data.frame(read_excel(pw_annot_file, .name_repair = "minimal"))
    rownames(pw_annot) = pw_annot$Pathway
    if (withNAs==FALSE)
      pw_annot = pw_annot[pw_annot$Annotation!="NA",]
      pw_annot_col =append( pw_annot_col, c("black"))
      names(pw_annot_col)[length(pw_annot_col)] = "NA"
    rel_pathways = pw_annot$Pathway
    rel_pathways = gsub("−", "-", rel_pathways) # something funky with the hyphens

    # keep only pathways that exist in the matrices
    keep <- rel_pathways %in% rownames(zscore_mat)
    if (any(!keep))
      {message("Dropping ", sum(!keep), " annotated pathways not found in filtered data, including: \n", paste(rel_pathways[!keep], "\n"))}
    pw_annot <- pw_annot[keep, , drop = FALSE]
    rel_pathways <- rel_pathways[keep]

    zscore_mat = zscore_mat[rel_pathways,]
    ast_mat = ast_mat[rel_pathways, ]
    row_annotation <- ComplexHeatmap::rowAnnotation(pathway_annotation=pw_annot$Annotation,
                                                    col = list(pathway_annotation=pw_annot_col),
                                                    annotation_name_gp = grid::gpar(fontfamily = font_family, fontface = "plain"),
                                                    annotation_legend_param = list(pathway_annotation = list(title_gp = grid::gpar(fontfamily = font_family, fontface = "plain"),
                                                                                                              labels_gp = grid::gpar(fontfamily = font_family, fontface = "plain"))))
    cluster_rows=FALSE
    if (row_display_mode == "pathway_annotation_group")
    {
      row_split_rows = factor(pw_annot$Annotation, levels = unique(pw_annot$Annotation))
    }
    else
    {
      row_split_rows = NULL
    }
  }
  else
  {
    row_annotation=NULL
    cluster_rows=TRUE
    row_split_rows = NULL
  }

  zscore_col <- circlize::colorRamp2(c(min(zscore_mat, na.rm = TRUE), 0, max(zscore_mat, na.rm = TRUE)),c("#1207A3", "#FFFEFE", "#BB0103"))
  zscore_mat[is.na(zscore_mat)] = 0  # so I can cluster rows

  ht <- ComplexHeatmap::Heatmap(zscore_mat,
                                cluster_rows = cluster_rows,
                                cluster_columns = FALSE,
                                name = "Zscore",
                                border_gp = grid::gpar(col = "black", lwd = 1),
                                rect_gp = grid::gpar(col = "black", lwd = 0.5),
                                top_annotation = Col_annotation,
                                left_annotation = row_annotation,
                                row_title_side = "right",
                                width = ncol(zscore_mat)*unit(5 ,"mm"),
                                height = nrow(zscore_mat)*unit(5, "mm"),
                                column_split = column_split_cols,
                                row_split = row_split_rows,
                                row_title_rot = 0,
                                row_names_side = "left",
                                row_names_gp = grid::gpar(fontfamily = font_family, fontface = "plain"),
                                column_names_gp = grid::gpar(fontfamily = font_family, fontface = "plain"),
                                row_title_gp = grid::gpar(fontfamily = font_family, fontface = "bold"),
                                column_title_gp = grid::gpar(fontfamily = font_family, fontface = "plain"),
                                col =zscore_col,
                                show_row_dend = FALSE,
                                column_title = paste0(plot_title, " Canonical Z= ", zscore_min_cp, " Canonical adjP=",padj_thresh_cp,  " Diseases Z= ", zscore_min_db, " Diseases adjP=",padj_thresh_db),
                                heatmap_legend_param = list(title_gp = grid::gpar(fontfamily = font_family, fontface = "plain"),
                                                            labels_gp = grid::gpar(fontfamily = font_family, fontface = "plain")),
                               cell_fun = function(j, i, x, y, w, h, fill) {
                                 if (!is.na(ast_mat[i, j]) && ast_mat[i, j] != "n.s.") {
                                   grid::grid.text(ast_mat[i, j], x, y, gp = grid::gpar(fontfamily = font_family, fontface = "plain"))
                                 }
                               })

  if (save_pdf)
  {
    grDevices::cairo_pdf(file = file.path(output_dir, paste0(output_prefix, ".pdf")),
                         width = pdf_width,
                         height = pdf_height,
                         family = font_family)
    ComplexHeatmap::draw(ht, annotation_legend_side = "top",
                         column_title = "Conditions",
                         column_title_side = "bottom",
                         column_title_gp = grid::gpar(fontfamily = font_family, fontface = "bold"))
    dev.off()
    message("Saved heatmap to PDF.")
  } else
  {
    ComplexHeatmap::draw(ht, annotation_legend_side = "top",
                         column_title = "Conditions",
                         column_title_side = "bottom",
                         column_title_gp = grid::gpar(fontfamily = font_family, fontface = "bold"))
  }
}






# Plotting 3 heatmaps for DEGs summary - raw number of DEGs, binary thereshold and a filtered heatmap without any rows that include only 0's
DEGS_heatmp_wrapper = function(DEGs_count_df, lin_subset_degs, output_str, clust_order)
{
  # Converting DEGs number to a table
  DEGs_count = dcast(DEGs_count_df, Lineage + Subset ~ Comparison, value.var = "DEG_no")
  DEGs_count[is.na(DEGs_count)] = 0
  # sort by lineage
  DEGs_count$Lineage = factor(DEGs_count$Lineage, levels = c("Myeloid","NK","T","B_Plasma","dStromal","Endothelial","fFB","fPV","Trophoblasts","EVT"), ordered=TRUE)
  DEGs_count = DEGs_count[order(DEGs_count$Lineage),]
  rownames(DEGs_count) = DEGs_count$Subset
  DEGs_count = DEGs_count[,-1]
  DEGs_count = DEGs_count[,-1]

  library(viridisLite)
  library(circlize)
  breaks <- unique(quantile(as.matrix(DEGs_count), probs = c(0, 0.7, 0.75, 0.8, 0.95, 0.99, 1),na.rm = TRUE))
  max_val <- max(as.matrix(DEGs_count), na.rm = TRUE)
  DEGs_count = DEGs_count[clust_order,]

  plot_heatmap(df=DEGs_count,
               lin_subset_df=lin_subset_degs,
               output_prefix=paste0("DEGs_test3_heatmap_",output_str),
               mat_name = "DEGs_no",
               #mat_values = breaks,
               #mat_colors = rev(viridisLite::magma(length(breaks))),
               mat_values = c(0, max_val),
               mat_colors = c("white", "black"),
               cluster_rows=F,
               pdf_width=7,
               pdf_height=15)
  write.table(DEGs_count, paste0("DEGs_count_pval0.15_abslog2FC0_ploting_raw.58",output_str,".csv"), row.names = T, col.names = T, quote = F, sep = ",")

}





# Plots a summary heatmap of DA across all lineages and conditions
plot_heatmap = function(df, lin_subset_df, output_prefix, mat_name, pdf_width=5, pdf_height=5, binary=TRUE, mat_values= NULL, mat_colors = NULL, cluster_rows=F)
{
  condition_df = data.frame(Comparison = colnames(df), Condition = sapply(colnames(df), function(x) strsplit(x, "-")[[1]][1]))
  colnames(df) = factor(colnames(df), levels = condition_df$Comparison, ordered = T)

  if (binary)
  {
    ht <- Heatmap(df,
          name=mat_name,
          cluster_columns = F,
          cluster_rows = cluster_rows,
          col=colorRamp2(mat_values, mat_colors),
          height = grid::unit(4, "in"),)
  }
  else
  {
    ht <- Heatmap(df,
      name=mat_name,
      cluster_columns = F,
      cluster_rows = cluster_rows,
      height = grid::unit(4, "in"),)
  }
  pdf(file = file.path(paste0(output_prefix, ".pdf")), width = pdf_width, height = pdf_height)
  ComplexHeatmap::draw(ht, annotation_legend_side = "top",
                           column_title = "Conditions",
                           column_title_side = "bottom",
                           column_title_gp = grid::gpar(fontface = "bold"))
  dev.off()
}



#################################################### DE #############################################################



# Runs DESeq2 per comparison pair per cluster
#
# Runs QC on data of a given lineage, in a loop subsets the data to the required complication and runs DESeq2 on each comparisons given in pairs_test, for every cluster separately
#
# pb_counts_file - the file containing pseudo-bulk counts per cluster per sample (created with a python code DifferentialExpression.py)
# meta_file - the meta information file for the counts (created with a python code DifferentialExpression.py)
# samp_condition_file - a csv file with the mapping of sample to condition
# cell_counts_file - the file containing the cell counts per sample per cluster to allow filtering clusters/ samples with too few cells (created with a python code DifferentialExpression.py)
# pairs_test - list of lists defining the pair of conditions to test
# run_title - lineage name, will appear in plots and file names
# cluster_col - the column names in seurat's meta-data attribute that contains the leiden clusters
# allCells - a boolean value indicating whether to test all cells together across conditions (as opposed to splititng by cluster)
# min_cells - the minimum amount of cells per cluster we require to run DE for that cluster
#
# When I refer to pb-samples in the documentation of this function I mean the unique ID in the PB data - that would be <sample>:<cluster>
#
DE_pregnancy_wrapper <- function(pb_counts_file, meta_file, cell_counts_file, samp_condition_file, pairs_test, run_title, cluster_col, allCells=FALSE,  plot_heatmap=FALSE, min_cells = 20)
{
    count_mtx <- as.matrix(read.csv(pb_counts_file, row.names = 1,check.names=FALSE))
    meta_data <- read.csv(meta_file, row.names = 1)
    cells_per_samp <- read.csv(cell_counts_file, row.names = 1)
    samp_condition = read.table(samp_condition_file, sep=",", header=TRUE)

    # updating the meta data rownames
    rownames(meta_data) <- paste0(meta_data$sample,":", meta_data$condition,":", meta_data$cluster)

    # Making sure that the meta file and the count file are ordered similarly
    .check_meta_and_count_match <- function(count_mtx, meta_data)
    {
       if (any(rownames(meta_data) != colnames(count_mtx)))
        {
          print("Error: meta file and count file are not ordered similarly")
          return(FALSE)
        }
        else
          return(TRUE)
    }

    if (!.check_meta_and_count_match(count_mtx, meta_data))
      return()

    clusters <- sort(unique(meta_data$cluster))

    for(pair in pairs_test)
    {
        p1 <- pair[1]
        p2 <- pair[2]
        curr_tit = paste0(run_title,"_",p1,"_vs_",p2)
        OUT <- openxlsx::createWorkbook()

        if (!allCells)
        {
            for (clust in clusters) # clust="DCs_1"
            {
                print(clust)
                model_res_file=paste0("tmp_DEseq_res/",curr_tit,"_clust_",clust,".rds")
                if(!file.exists(model_res_file))
                {
                    # Making sure we have enough cell per sample per cluster
                    df <- cells_per_samp[cells_per_samp[cluster_col] == clust,]
                    samp_past_filter <- df[df$count >= min_cells, "sample"]

                    if (length(samp_past_filter) > 6) # we need at least 6 samples that are valid for comparison
                    {
                        # focusing on the bp-samples that passed the filter + specific cluster
                        clust_meta <- meta_data[(meta_data$cluster == clust) & (meta_data$sample %in% samp_past_filter),]
                        # Focsuing on the current 2 conditions
                        current_samples = samp_condition[samp_condition$Condition_group %in% c(p1, p2),"Sample"]
                        clust_meta <- clust_meta[clust_meta$sample %in% current_samples,]
                        if (dim(clust_meta)[1] > 6) # we need at least 6 samples that are valid for comparison
                        {
                            # Creating a column in the meta data df that represents the comparison we are currently testing
                            clust_meta$cond_groups = p1
                            p2_samples = samp_condition[samp_condition$Condition_group == p2,"Sample"]
                            clust_meta[clust_meta$sample %in%  p2_samples, "cond_groups"] = p2

                            # Filtering the count matrix too
                            clust_count <- count_mtx[,rownames(clust_meta)]

                            if (!.check_meta_and_count_match(clust_count, clust_meta))
                                return()

                            # Do we have at least 3 samples per condition in the pair of conditions we are currently comparing?
                            enough_samples_per_cond <- (dim(clust_meta)[1] > 4) & (all(table(clust_meta$cond_groups) > 2) & (length(unique(clust_meta$cond_groups))==2))

                            if (enough_samples_per_cond)
                            {
                                print(paste0("  ===> Comparing ",p1," vs ",p2," - Cluster ", clust))
                                Running_DESeq_per_comp(clust_count, clust_meta, OUT, model_res_file, p1, p2, sheet_name=clust, contrast_column="cond_groups")
                                openxlsx::saveWorkbook(OUT, paste0("DEGs_",curr_tit,".xlsx"), overwrite=TRUE)
                            }
                            else
                              {print(paste0("  ===> Not enough valid samples per condition for ",p1," vs ",p2," - Cluster ", clust))}
                        }
                    else
                       {print(paste0("  ===> Not enough valid samples per condition for ",p1," vs ",p2," - Cluster ", clust))}
                }
                  else
                    {print(paste0("  ===> Not enough valid samples per condition for ",p1," vs ",p2," - Cluster ", clust))}
                }
                else
                    {print(paste0("  ===> Already done ",p1," vs ",p2," - Cluster ", clust))}
            }
        }
        else
        {
            model_res_file=paste0("tmp_DEseq_res/",curr_tit,"_all.rds")
            print(paste0("  ===> Comparing ",p1," vs ",p2," - all cells"))
            Running_DESeq_per_comp(count_mtx, meta_data, OUT, clust, model_res_file, p1, p2, sheet_name="all", contrast_column="cond_groups")
            openxlsx::saveWorkbook(OUT, paste0("DEGs_",curr_tit,".xlsx"), overwrite=TRUE)
        }
    }

}




# Running DESEq2 for a given count data and meta data
#
# OUT - workbook object to write the results to, created by openxlsx::createWorkbook()
# clust - names of cuurent cluster or "all" in case all cells of a lineage are tested together
# model_res_file - a name of a file in which to write the model results per cluster.
# p1, p2 - the conditions that are currently being contrasted
Running_DESeq_per_comp <- function(clust_count, clust_meta, OUT, model_res_file, p1, p2, sheet_name, contrast_column="condition")
{
    clust_meta[,contrast_column] <- factor(clust_meta[,contrast_column], levels = c(p2, p1))
    if (contrast_column == "condition")
    {
      dds <- DESeqDataSetFromMatrix(countData = clust_count, colData = clust_meta, design = ~ condition)
    }
    else if (contrast_column == "cond_groups")
    {
      dds <- DESeqDataSetFromMatrix(countData = clust_count, colData = clust_meta, design = ~ cond_groups)
    }
    else
    {
      print("Error: <contrast_column> can be 'condition' or 'cond_gorups' only" )
      return()
    }
    #keep <- rowSums(counts(dds)) >= 5 # Pre-filtering genes to reduce number of tests
    #dds <- dds[keep,]
    dds <- DESeq(dds)

    print("DESEQ2 default result:")
    # p2: baseline, p1: level to compare:
    res <- results(dds, contrast=c(contrast_column, p1, p2))
    # more strict filtering is applied on the mean of normalized counts within the results function
    res$abslog2FoldChange=abs(res$log2FoldChange)
    res$norm_sumcounts=rowSums(counts(dds,normalized=T))
    print("calculating variance stabiliziing data:")
    print("calculating normalized variance per gene:")

    dir.create("tmp_DEseq_res", showWarnings = FALSE)
    saveRDS(res, model_res_file) # res = readRDS("res_clust1.rds")

    res <- data.frame(res)
    res=res[order(res$padj),]
    openxlsx::addWorksheet(OUT, sheet_name)
    openxlsx::writeData(OUT, sheet = as.character(sheet_name), x = res, rowNames =TRUE)
}







summarize_DE_results <- function(path, title, comparisons , lineages)
{
  all_DEGs = data.frame()

  print("Aggregating all DE results...")
  for (lin in lineages)
  {
    print(lin)
    for (contrast in comparisons)
    {
      file_name = list.files(path=path, pattern=paste0("^DEGs_", lin, ".*_", contrast, ".xlsx"))
      if (length(file_name)>1) stop(paste0("More than one file found in ", path))

      subsets = excel_sheets(path = paste0(path, file_name))
      for (subset in subsets)
      {
        print(subset)
        df = as.data.frame(read_excel(paste0(path, file_name), sheet = subset))
        colnames(df)[1] = "Gene"
        df["contrast"] = contrast
        df["cluster"] = subset
        df["Lineage"] = lin
        all_DEGs = rbind(all_DEGs, df[,c("Gene","Lineage","cluster","contrast","baseMean","log2FoldChange","pvalue","padj","abslog2FoldChange")])
      }
    }
  }
  write.csv(all_DEGs, paste0("AllDEGs_",title,".csv"), row.names = FALSE)

  print("Plot bar plots...")
  all_DEGs["FDR_0.05"] = all_DEGs$padj < 0.05
  all_DEGs[is.na(all_DEGs["FDR_0.05"]), "FDR_0.05"] = FALSE
  all_DEGs["FDR_0.1"] = all_DEGs$padj < 0.1
  all_DEGs[is.na(all_DEGs["FDR_0.1"]), "FDR_0.1"] = FALSE

  n_DEGs = all_DEGs %>% dplyr::group_by(cluster, contrast) %>%
                              dplyr::reframe("0.05"=sum(FDR_0.05),
                                               "0.1"=sum(FDR_0.1),
                                               "Lineage"=unique(Lineage)) %>% arrange(Lineage) %>%
                              tidyr::pivot_longer(names_to="FDR",values_to="n_DE",c("0.05","0.1"))
  n_DEGs$FDR=as.factor(n_DEGs$FDR)

  write.csv(n_DEGs, paste0("nDEGs_",title,".csv"), row.names = FALSE)

  n_DEGs

  pdf(paste0("DEGs_barplot_",title,".pdf"), width=10, height=5)
  for(contrast in comparisons)
  {
    curr_df = n_DEGs[n_DEGs$contrast==contrast,]
    print(ggplot(data=curr_df, aes(x=n_DE, y=fct_reorder(cluster, n_DE), fill=FDR)) +
        geom_bar(stat="identity", position="dodge")+
        #xlim(0,1100)+
        facet_nested(~Lineage,scale="free", space="free_x") +
        geom_vline(xintercept = 0)+coord_flip()+
        theme_classic()+
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.line.x = element_blank())+ scale_fill_brewer(palette = "Set1")+
        ggtitle(contrast))
  }
  dev.off()

}
