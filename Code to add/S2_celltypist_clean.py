import scanpy as sc 
import celltypist as ct 
import pandas as pd 
import muon as mu 
import time 
import numpy as np 
from celltypist import models


# ----------
# ------ dir org 
projdir="/projects/home/rpragya/pregnancy_plac/label_comparisons/"
inputdir=projdir+"data/"
savedir=projdir+"results/"
placdir="/projects/Pregnancy_placenta/Results/21_Final_objects/GEO_objects/"
preddir=savedir+"predictions/"


# ----------
# ------ 1. predict normand labels using vento-tormo et al 2018 dataset 
# load internal data 
mudata = mu.read(placdir+"AllCells.h5mu")
adata = mudata['GEX'].copy()

# internal data is normalised to 1e6, celltypist needs 1e4, renormalise
adata.layers['counts'] = adata.raw.X.copy() # save raw counts
tmp = sc.AnnData(adata.layers['counts'].copy(), obs=adata.obs, var=adata.var)
sc.pp.normalize_total(tmp, target_sum=1e4)
sc.pp.log1p(tmp)
# modify var_names to match model
tmp.var['featurekey'] = tmp.var_names.copy()
tmp.var['featurekey'] = tmp.var['featurekey'].map(lambda x: x.replace("G:", ""))
tmp.var_names = tmp.var['featurekey'].copy()

# run classifier 
t_start = time.time()
predictions = ct.annotate(tmp, model='Human_Placenta_Decidua.pkl', majority_voting=True)
t_end = time.time()
print(f"Time elapsed: {t_end - t_start}/60 minutes") # 591 seconds

# save predictions for plotting
predictions.predicted_labels.to_csv(preddir+"normand_pred_using_vt/predicted_labels.csv")
predictions.decision_matrix.to_csv(preddir+"normand_pred_using_vt/decision_matrix.csv")
predictions.probability_matrix.to_csv(preddir+"normand_pred_using_vt/probability_matrix.csv")
tmp.obs[['Cell_subset', 'Lineage']].to_csv(preddir+"normand_pred_using_vt/reference.csv")


# ----------
# ------ 2. predict normand labels using shook et al 2024 dataset 
# -- load shook data 
shook = sc.read(inputdir+"batorsky_2024/raw_data_w_annots.h5ad")

# normalize to 1e4
shook.layers['counts'] = shook.X.copy() # save copy of raw counts 
sc.pp.normalize_total(shook, target_sum=1e4)
sc.pp.log1p(shook)

# train model on annotated cell subsets 
t_start = time.time()
model = ct.train(shook, labels='ident',
                 top_genes=300, 
                 n_jobs=10, max_iter=100)
t_end = time.time()
print(f"Time elapsed: {(t_end - t_start)/60} minutes") # 0.6 mins 
# save model 
model.write(savedir+'models/shook_HBC_allfeatures.pkl')

# filter internal data to only hbcs, pamms and mmacs for comparison with shook et al
celltypes_list = ['mMac', 'HBC', 'Monocytes']
adata_subset = tmp[tmp.obs.Lineage.isin(celltypes_list)].copy()

# CellTypist prediction with over-clustering and majority-voting.
t_start = time.time()
pred_shook = ct.annotate(adata_subset, model = savedir+'models/shook_HBC_allfeatures.pkl', majority_voting = True)
t_end = time.time()
print(f"Time elapsed: {(t_end - t_start)/60} minutes") # 14 minutes 

# save pred_shook  
pred_shook.predicted_labels.to_csv(preddir+"normand_pred_using_shook/predicted_labels.csv")
pred_shook.decision_matrix.to_csv(preddir+"normand_pred_using_shook/decision_matrix.csv")
pred_shook.probability_matrix.to_csv(preddir+"normand_pred_using_shook/probability_matrix.csv")
adata_subset.obs[['Cell_subset', 'Lineage']].to_csv(preddir+"normand_pred_using_shook/reference.csv")

# visualise predictions using celltypist dotplot wrapper 
# dp = ct.dotplot(pred_shook, use_as_reference = 'Cell_subset', use_as_prediction = 'majority_voting', return_fig=True, show=False)
# dp.savefig(savedir+"figures/normand_pred_using_shook.pdf", bbox_inches='tight') # majority voting
# dp = ct.dotplot(pred_shook, use_as_reference = 'Cell_subset', use_as_prediction = 'predicted_labels', return_fig=True, show=False)
# dp.savefig(savedir+"figures/normand_pred_using_shook_all.pdf", bbox_inches='tight') # predicted_labels


# ----------
# ------ 3. train normand model and predict vt labels using it 
immune = ['mMac','HBC','CD14_Mon','CD16_Mon','cDC','pDC','NK','T_naive','CD4_eff',
          'CD8_eff','Treg','T_unconventional','T_cycling','B_cells','Plasma']
adata_imm = tmp[tmp.obs['Parent_lineage'].isin(immune)].copy()

# downsample data to 3000 cells per subset (cell subsets with <3000 cells will keep all cells)
train = ct.samples.downsample_adata(adata_imm, mode = 'each', n_cells = 3000, by = 'Cell_subset', random_state=0, return_index = False)
train.obs.to_csv(inputdir+"normand_training_set.csv")

t_start = time.time()
model = ct.train(train, labels='Cell_subset', use_SGD = True, 
                mini_batch = True, balance_cell_type=True,
                n_jobs=10, max_iter=100)
t_end = time.time()
print(f"Time elapsed: {(t_end - t_start)/60} minutes") # 13 mins 
model.write(savedir+'models/normand_immune_downsampled_w_SGD.pkl')

# test model by predicting cells in internal data not used for training
train_bc = train.obs_names
test = adata_imm[~(adata_imm.obs_names.isin(train_bc))].copy()
# predict 
t_start = time.time()
pred_test = ct.annotate(test, model = savedir+'models/normand_immune_downsampled_w_SGD.pkl', majority_voting = True)
# 29670 features used for prediction
t_end = time.time()
print(f"Time elapsed: {t_end - t_start} seconds") # 32 seconds

# visualise results 
dp = ct.dotplot(pred_test, use_as_reference = 'Cell_subset', use_as_prediction = 'predicted_labels', return_fig=True, show=False)
dp.savefig(savedir+"figures/normand_pred_using_normand_test.pdf", bbox_inches='tight') # predicted_labels


# predict vt labels using internal data 
# load vt data 
vt_raw = sc.read(inputdir+"ventotormo_2018/raw_data.h5ad")
# raw data contains data from pbmcs, filter that
vt_filt = vt_raw[~(vt_raw.obs.location=='Blood')].copy()
# normalise 
vt_filt.layers['counts'] = vt_filt.X.copy()
scales_counts = sc.pp.normalize_total(vt_filt, target_sum=None, inplace=False)
# log1p transform
vt_filt.layers["log1p_norm"] = sc.pp.log1p(scales_counts["X"], copy=True)
vt_filt.X = vt_filt.layers["log1p_norm"].copy() 
# make var compatible with model 
split = vt_filt.var_names.str.split('_', n=1, expand=True)
vt_filt.var['gene_symbol'] = split.get_level_values(0)
vt_filt.var['ensembl_id'] = split.get_level_values(1)
vt_filt.var['old_var_names'] = vt_filt.var_names.copy()
vt_filt.var_names = vt_filt.var['gene_symbol'].copy()
vt_filt.var_names_make_unique()
vt_filt.var.index.name = 'featurekey'

# predict 
t_start = time.time()
predictions = ct.annotate(vt_filt, model = savedir+'models/normand_immune_downsampled_w_SGD.pkl')
# 19362 features used for prediction
t_end = time.time()
print(f"Time elapsed: {t_end - t_start} seconds") # 21 seconds

# visualise results 
dp = ct.dotplot(predictions, use_as_reference = 'annotation', use_as_prediction = 'predicted_labels', return_fig=True, show=False)
dp.savefig(savedir+"figures/vt_pred_using_normand_downsampled_w_SGD.pdf", bbox_inches='tight') # predicted_labels


# downsample keeping the proportions constant 