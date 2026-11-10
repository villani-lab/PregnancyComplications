"""Figure S11 correction, classification, cross-K, plotting and final exports.

Intermediate arrays remain in memory. Original cell order and the final
within-lineage permutation implementation are preserved. This module never
imports Cellpose or launches a worker.
"""
from pathlib import Path
import hashlib
import json
import re
from io import BytesIO

import anndata as ad
import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial.distance import cdist
from numba import njit

import figs11_segmentation_helpers as segmentation

MARKERS = segmentation.MARKERS
GATED_MARKERS = ('CD3', 'CD8A', 'Granzyme_B', 'CD163', 'CD274')
PRIMARY_SUBTYPE = "CD8A+/Granzyme B+ T"
SCHEMES = [dict(name="global_both", title="Both subtype labels shuffled", shuffle_t=True, shuffle_m=True)]


# 1. Load samples and thresholds
# ----------------------------------------------------------------------------

def dense(matrix):
    return matrix.toarray() if sparse.issparse(matrix) else np.asarray(matrix)

def json_value(value):
    """Convert NumPy provenance to JSON, representing unavailable numbers as null."""
    if isinstance(value, dict):
        return {str(key): json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_value(item) for item in value]
    if isinstance(value, np.ndarray):
        return json_value(value.tolist())
    if isinstance(value, np.generic):
        return json_value(value.item())
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value

def read_thresholds(path, sample_id):
    path = segmentation.require_orcd_path(path)
    content = path.read_bytes()
    document = json.loads(content)
    if document.get('schema_version') != 1:
        raise ValueError('Expected Figure S11 thresholds schema_version=1.')
    sample = document['samples'][sample_id]
    if set(sample['markers']) != set(MARKERS[1:]):
        raise ValueError('Thresholds must contain the six reviewed proteins.')
    for marker, setting in sample['markers'].items():
        q, cutoff = setting['quantile'], setting['threshold_raw_units']
        if not np.isfinite(q) or not 0 <= q <= 1 or not np.isfinite(cutoff) or cutoff < 0:
            raise ValueError(f'Invalid pixel threshold for {marker}.')
        if marker in GATED_MARKERS:
            if not np.isfinite(setting['z_cutoff']) or setting['operator'] not in ('>', '>='):
                raise ValueError(f'Invalid cell gate for {marker}.')
        elif setting['z_cutoff'] is not None or setting['operator'] is not None:
            raise ValueError('CD68 is retained for display and must not define a lineage.')
    return sample, dict(path=str(path), sha256=hashlib.sha256(content).hexdigest(),
                        schema_version=document['schema_version'])

def load_cohort(thresholds_path):
    """Read all sample identifiers from the supplied JSON, keeping file order."""
    document = json.loads(Path(thresholds_path).read_text())
    if document.get('schema_version') != 1 or not document.get('samples'):
        raise ValueError('Expected the supplied Figure S11 threshold JSON.')
    samples = []
    for sample_id in document['samples']:
        match = re.fullmatch(r'TMA(\d+)_ROI(\d+)_Row(\d+)', sample_id)
        if match is None:
            raise ValueError(f'Invalid sample identifier: {sample_id}')
        tma, roi, row = map(int, match.groups())
        segmentation.sample_identity(tma, roi, row)
        read_thresholds(thresholds_path, sample_id)
        samples.append(dict(sample_id=sample_id, TMA=tma, ROI=roi, ROW=row))
    return samples

def load_sample(tma, roi, row, *, segmentation_run, thresholds_path,
                manuscript_segmentation=None, image_relative=None):
    """Load raw segmentation geometry and flags, never a thresholded V1/V2 table.

    The segmentation run record identifies a new export or the retained
    manuscript inventory. Both preserve their saved cell IDs, order and crop.
    """
    segmentation.require_compute()
    sample_id = segmentation.sample_identity(tma, roi, row)
    folder = segmentation.result_directory(segmentation_run, sample_id)
    table_path = folder / 'all_detections_with_eligibility.csv'
    geometry_path = folder / 'all_expanded_geometry_with_eligibility.parquet'
    record_path = folder / 'provenance.json'
    table = pd.read_csv(table_path, dtype={'analysis_cell_id': str, 'source_cell_id': str},
                        keep_default_na=False)
    table.index = pd.Index(table.analysis_cell_id, name='cell_id')
    if table.index.has_duplicates or (table.index.str.len() == 0).any():
        raise ValueError('Segmentation cell identifiers must be unique and nonempty.')
    if not pd.api.types.is_bool_dtype(table.analysis_eligible):
        raise TypeError('analysis_eligible must contain boolean flags.')
    geometry = gpd.read_parquet(geometry_path)
    geometry.index = pd.Index(geometry.analysis_cell_id.astype(str), name='cell_id')
    if geometry.index.has_duplicates:
        raise ValueError('Duplicate geometry identifiers.')
    geometry = geometry.reindex(table.index)
    eligible = table.analysis_eligible.to_numpy(bool)
    if not eligible.any() or geometry.loc[eligible].geometry.isna().any():
        raise ValueError('Every eligible cell needs its saved geometry.')
    if not np.array_equal(geometry.loc[eligible].nucleus_label, table.loc[eligible].nucleus_label):
        raise ValueError('Nuclear labels and ordered cell geometry disagree.')
    if 'analysis_eligible' in geometry and not np.array_equal(geometry.analysis_eligible, table.analysis_eligible):
        raise ValueError('Geometry and inventory eligibility flags disagree.')
    record = json.loads(record_path.read_text())
    if manuscript_segmentation is None:
        manuscript_segmentation = record.get('schema') != 'FigS11_segmentation_v1'
    image = segmentation.selected_image(tma, roi, image_relative=image_relative)
    # Historical path strings are not used to select, match or open files.
    scale = float(record['scale_um_per_pixel'])
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError('Saved positive cell calibration is required.')
    channels, shape = segmentation.inspect_image_channels(image)
    crop = record.get('saved_crop')
    bounds = segmentation.checked_bounds(None if crop is None else crop['bounds_x0_x1_y0_y1'], shape)
    required_crop = segmentation.recorded_crop(sample_id)
    if required_crop is not None and bounds != required_crop:
        raise ValueError('TMA10 ROI4 requires its recorded manuscript crop.')
    if 'eligible' in record and int(record['eligible']) != int(eligible.sum()):
        raise ValueError('Saved eligible count disagrees with the cell inventory.')
    source = None
    source_path = None
    source_hash = None
    if manuscript_segmentation:
        if record.get('schema') == 'FigS11_segmentation_v1':
            raise ValueError('Use manuscript_segmentation=False for a new segmentation export.')
        if sample_id != 'TMA8_ROI7_Row5':
            source_path = (segmentation.PROJECTS / 'Preg_proteomics_review/Data/3_segmentation_dapi/3_segmentation'
                           / f'TMA{tma}_ROI{roi}' / 'data/tables/cell_table_dapi.h5ad')
    else:
        if record.get('schema') != 'FigS11_segmentation_v1' or record.get('sample_id') != sample_id:
            raise ValueError('Expected a matching Figure S11 segmentation export.')
        source_path = folder / 'raw_cells.h5ad'
    if source_path is not None:
        source_bytes = source_path.read_bytes()
        source_hash = hashlib.sha256(source_bytes).hexdigest()
        source = ad.read_h5ad(BytesIO(source_bytes))
        del source_bytes
        if 'raw' not in source.layers or source.var_names.has_duplicates:
            raise ValueError('The original segmentation requires its raw layer and unique channels.')
    thresholds, threshold_record = read_thresholds(thresholds_path, sample_id)
    inputs = dict(sample_id=sample_id, TMA=int(tma), ROI=int(roi), ROW=int(row),
                  registered_image=str(image), segmentation_run=segmentation_run,
                  segmentation_inventory=str(table_path), geometry=str(geometry_path),
                  segmentation_inventory_sha256=hashlib.sha256(table_path.read_bytes()).hexdigest(),
                  geometry_sha256=hashlib.sha256(geometry_path.read_bytes()).hexdigest(),
                  segmentation_provenance=str(record_path),
                  segmentation_provenance_sha256=hashlib.sha256(record_path.read_bytes()).hexdigest(),
                  original_source_table=None if source_path is None else str(source_path),
                  original_source_table_sha256=source_hash,
                  image_shape_yx=list(shape), analysis_window_pixels=(np.asarray(bounds, float)-.5).tolist(),
                  scale_um_per_pixel=scale, manuscript_segmentation=bool(manuscript_segmentation),
                  thresholds=threshold_record, sample_thresholds=thresholds)
    return dict(table=table, geometry=geometry, source=source, inputs=inputs,
                thresholds=thresholds, channels=channels)


# 2. Correct pixel backgrounds and measure cell intensities
# ----------------------------------------------------------------------------

def measure_corrected_cells(sample):
    """Measure saved cutoffs once; retain all detections and original raw channels."""
    table, geometry, inputs = sample['table'], sample['geometry'], sample['inputs']
    valid = geometry.geometry.notna() & geometry.geometry.is_valid & ~geometry.geometry.is_empty
    regions, counts = segmentation.native_pixel_regions(geometry.loc[valid], inputs['image_shape_yx'])
    native = np.full((len(table), len(MARKERS)), np.nan, dtype=np.float64)
    corrected = np.full_like(native, np.nan)
    fractions = np.full_like(native, np.nan)
    specs = segmentation.marker_acquisitions(inputs['sample_id'])
    acquisitions = {}
    for column, marker in enumerate(MARKERS):
        if marker == 'DAPI':
            plane, acquisition = segmentation.load_segmentation_dapi(inputs['registered_image'])
            expected_cycle = 0 if inputs['sample_id'] == 'TMA7_ROI7_Row5' else 1
            if acquisition['cycle'] != expected_cycle:
                raise ValueError('DAPI acquisition differs from the manuscript record.')
            cutoff = 0.
        else:
            images, metadata = segmentation.load_raw_markers(inputs['registered_image'],
                                                              {marker: specs[marker]}, require_cycles=True)
            plane = images[marker]
            acquisition = metadata.reset_index().to_dict('records')[0]
            cutoff = float(sample['thresholds']['markers'][marker]['threshold_raw_units'])
        raw_values, corrected_values, positive_fraction = [], [], []
        for index, ys, xs, mask in regions:
            pixels = plane[ys, xs][mask].astype(np.float64)
            raw_values.append(pixels.sum() / counts[index])
            corrected_values.append(np.maximum(pixels-cutoff, 0).sum() / counts[index])
            positive_fraction.append(np.count_nonzero(pixels > cutoff) / counts[index])
        native[valid.to_numpy(), column] = raw_values
        corrected[valid.to_numpy(), column] = corrected_values
        fractions[valid.to_numpy(), column] = positive_fraction
        acquisitions[marker] = json_value(acquisition)
        del plane
        if marker != 'DAPI':
            del images
    obs = table.copy()
    obs['cell_id'] = obs.index.astype(str)
    obs['area'] = geometry.geometry.area.to_numpy()
    obs['area_um2'] = obs.area * inputs['scale_um_per_pixel']**2
    for name in ('TMA', 'ROI', 'ROW'):
        obs[name] = inputs[name]
    obs['sample_name'] = f"TMA{inputs['TMA']}_ROI{inputs['ROI']}"
    obs['measurement_pixel_count'] = 0
    obs.loc[valid, 'measurement_pixel_count'] = counts
    obs['cell_polygon_wkb_hex'] = geometry.geometry.to_wkb(hex=True).fillna('').to_numpy()
    result = ad.AnnData(X=corrected, obs=obs, var=pd.DataFrame(index=pd.Index(MARKERS, name='marker')))
    result.layers['uncorrected_native'] = native
    result.layers['fraction_pixels_above_background'] = fractions
    # Preserve original segmentation measurements without assigning old values to new detections.
    source = sample['source']
    if source is None:
        full_raw, raw_var, available = native.copy(), result.var.copy(), valid.to_numpy()
        result.layers['raw'] = native.copy()
    else:
        original = dense(source.layers['raw'])
        source_ids = (source.obs.cell_id.astype(str) if 'cell_id' in source.obs
                      else pd.Index(source.obs_names.astype(str)))
        lookup = pd.Series(np.arange(source.n_obs), index=source_ids)
        if lookup.index.has_duplicates:
            raise ValueError('Original raw cell identifiers are duplicated.')
        available = table.source_cell_id.astype(str).ne('').to_numpy()
        selected_ids = table.loc[available, 'source_cell_id'].astype(str)
        if selected_ids.duplicated().any() or set(selected_ids) != set(source_ids):
            raise ValueError('Every original source cell must occur once in the detection inventory.')
        positions = lookup.loc[selected_ids].to_numpy(int)
        full_raw = np.full((len(table), source.n_vars), np.nan,
                           dtype=np.result_type(original.dtype, np.float32))
        full_raw[available] = original[positions]
        raw_var = source.var.copy()
        marker_positions = source.var_names.get_indexer(MARKERS)
        # Newly measured exports store every image channel under a unique index.
        result.layers['raw'] = (full_raw[:, marker_positions].copy()
                                if (marker_positions >= 0).all() else native.copy())
    result.obs['raw_layer_available'] = available
    result.raw = ad.AnnData(X=full_raw, obs=pd.DataFrame(index=result.obs_names.copy()), var=raw_var)
    centers = geometry.geometry.centroid
    result.obsm['spatial'] = np.column_stack([centers.x.to_numpy(), centers.y.to_numpy()])
    result.obsm['spatial_um'] = result.obsm['spatial'] * inputs['scale_um_per_pixel']
    result.uns['sample_id'] = inputs['sample_id']
    result.uns['analysis_window_pixels'] = np.asarray(inputs['analysis_window_pixels'])
    result.uns['scale_um_per_pixel'] = inputs['scale_um_per_pixel']
    result.uns['inputs_json'] = json.dumps(json_value(inputs), sort_keys=True)
    result.uns['acquisitions_json'] = json.dumps(acquisitions, sort_keys=True)
    result.uns['thresholds_json'] = json.dumps(sample['thresholds'], sort_keys=True)
    result.uns['measurement_semantics'] = 'mean(max(native pixel - saved cutoff, 0)); denominator includes corrected zeros'
    eligible = result.obs.analysis_eligible.to_numpy(bool)
    if not np.isfinite(result.X[eligible]).all():
        raise ValueError('An eligible cell has missing measurements.')
    return result


# 3. Apply cell gates and assign parent populations
# ----------------------------------------------------------------------------

def assign_margin_lineage(t_gate, macrophage_gate, t_margin, macrophage_margin):
    """Compare margins only for cells passing both gates; exact ties stay unresolved."""
    t_gate, macrophage_gate = np.asarray(t_gate, bool), np.asarray(macrophage_gate, bool)
    t_margin, macrophage_margin = np.asarray(t_margin, float), np.asarray(macrophage_margin, float)
    if len({v.shape for v in [t_gate, macrophage_gate, t_margin, macrophage_margin]}) != 1:
        raise ValueError("Gates and scores must align.")
    if not np.isfinite(t_margin).all() or not np.isfinite(macrophage_margin).all():
        raise ValueError("Nonfinite margins require review.")
    both = t_gate & macrophage_gate
    return np.select(
        [both & (t_margin > macrophage_margin), both & (t_margin < macrophage_margin),
         both, t_gate, macrophage_gate],
        ["T_cell", "mMac_candidate", "unresolved_tie", "T_cell", "mMac_candidate"],
        default="neither",
    )

def classify_cells(adata, *, cofactor=5., ddof=1):
    """Apply the final saved gates in memory, fitting only eligible cells."""
    if cofactor != 5. or ddof != 1:
        raise ValueError('The manuscript uses arcsinh cofactor 5 and ddof=1.')
    settings = json.loads(adata.uns['thresholds_json'])['markers']
    eligible = adata.obs.analysis_eligible.to_numpy(bool)
    if eligible.sum() < 2:
        raise ValueError('At least two eligible cells are needed for sample SD scaling.')
    # Preserve the reviewed six-protein array selection and reduction order.
    indices = adata.var_names.get_indexer(MARKERS[1:])
    values = dense(adata.X)[eligible][:, indices].astype(float)
    protein_transformed = np.arcsinh(values / cofactor)
    means, sds = protein_transformed.mean(axis=0), protein_transformed.std(axis=0, ddof=ddof)
    if not np.isfinite(sds).all() or (sds <= 0).any():
        raise ValueError('A constant/nonfinite protein requires review before applying z-score gates.')
    transformed = np.full(adata.shape, np.nan, dtype=float)
    z = np.full_like(transformed, np.nan)
    transformed[np.ix_(eligible, indices)] = protein_transformed
    z[np.ix_(eligible, indices)] = (protein_transformed - means) / sds
    calls = {}
    for marker in GATED_MARKERS:
        index = adata.var_names.get_loc(marker)
        setting = settings[marker]
        call = z[:, index] > setting['z_cutoff'] if setting['operator'] == '>' else z[:, index] >= setting['z_cutoff']
        calls[marker] = call & eligible
        adata.obs['positive_' + marker] = calls[marker]
    lineage = np.full(adata.n_obs, 'excluded', dtype=object)
    t_margin = z[:, adata.var_names.get_loc('CD3')] - settings['CD3']['z_cutoff']
    m_margin = z[:, adata.var_names.get_loc('CD163')] - settings['CD163']['z_cutoff']
    lineage[eligible] = assign_margin_lineage(calls['CD3'][eligible], calls['CD163'][eligible],
                                             t_margin[eligible], m_margin[eligible])
    t, m = lineage == 'T_cell', lineage == 'mMac_candidate'
    adata.obs['resolved_lineage'] = lineage
    states = dict(T_cell=t, mMac_candidate=m, CD8_T_cell=t & calls['CD8A'],
                  GZB_T_cell=t & calls['Granzyme_B'],
                  CD8_GZB_T_cell=t & calls['CD8A'] & calls['Granzyme_B'],
                  CD274_mMac_candidate=m & calls['CD274'], unresolved_tie=lineage == 'unresolved_tie')
    for name, value in states.items():
        adata.obs['resolved_' + name] = value
    adata.obs['CD3_CD163_double_positive_before_assignment'] = calls['CD3'] & calls['CD163']
    adata.layers['arcsinh_c5'] = transformed
    adata.layers['zscore_c5_ddof1'] = z
    adata.uns['classification_json'] = json.dumps(json_value(dict(
        cofactor=cofactor, ddof=ddof, fit_population='all eligible cells within sample',
        marker_order=list(MARKERS[1:]), arcsinh_mean=means, arcsinh_sd=sds,
        lineage_policy='larger_z_margin_exact_ties_unresolved')))
    return adata

def population_counts(adata):
    populations = dict(Detected=adata.n_obs, Eligible=int(adata.obs.analysis_eligible.sum()),
                       T_cells=int(adata.obs.resolved_T_cell.sum()),
                       CD8A_GZMB_T_cells=int(adata.obs.resolved_CD8_GZB_T_cell.sum()),
                       Macrophages=int(adata.obs.resolved_mMac_candidate.sum()),
                       CD274_macrophages=int(adata.obs.resolved_CD274_mMac_candidate.sum()),
                       CD3_CD163_before_assignment=int(adata.obs.CD3_CD163_double_positive_before_assignment.sum()))
    return pd.DataFrame({'Population': list(populations), 'Cells': list(populations.values())})


# 4. Run the Figure S11 cross-K analysis
# ----------------------------------------------------------------------------

def analyze_crossk(adata, *, seed, n_permutations, radii_um):
    if seed != 20260930 or n_permutations != 9999 or not np.array_equal(radii_um, np.arange(10., 151., 5.)):
        raise ValueError('Preserve manuscript seed 20260930, 9,999 permutations and radii 10:5:150 µm.')
    analyzed = adata[adata.obs.analysis_eligible.to_numpy(bool)].copy()
    global_stats, radius_stats, details = run_scheme(
        analyzed, SCHEMES[0], window_pixels=adata.uns['analysis_window_pixels'],
        seed=seed, n_permutations=n_permutations, radii_um=radii_um, return_details=True)
    test = upper_global_test(details['K'], details['expected_K'])
    global_stats['max_T_radius_um'] = details['radii_um'][int(np.argmax(test['standardized_observed']))]
    global_stats.insert(0, 'sample_id', adata.uns['sample_id'])
    radius_stats.insert(0, 'sample_id', adata.uns['sample_id'])
    for column in ('TMA', 'ROI', 'ROW'):
        global_stats[column] = int(adata.obs[column].iloc[0])
    global_stats['eligible_cells'] = analyzed.n_obs
    global_stats['CD3_CD163_before_assignment'] = int(adata.obs.CD3_CD163_double_positive_before_assignment.sum())
    global_stats['post_hoc_gate_revision'] = bool(json.loads(adata.uns['thresholds_json'])['post_hoc_gate_revision'])
    # Threshold settings remain in JSON, not in the Table S20 statistics export.
    global_stats = global_stats.drop(columns=['seed', 'simulations'])
    details['ordered_cell_ids'] = analyzed.obs_names.to_numpy(str)
    return dict(core_stats=global_stats, radius_stats=radius_stats, details=details)

def prepare_inputs(adata, *, window_pixels, radii_um=None):
    radii = np.arange(10., 151., 5.) if radii_um is None else np.asarray(radii_um, float)
    if radii.ndim != 1 or len(radii) == 0 or not np.isfinite(radii).all() or np.any(np.diff(radii) <= 0) or radii[0] <= 0:
        raise ValueError("Use increasing positive radii.")
    pixel_sizes = np.sqrt(adata.obs.area_um2.to_numpy() / adata.obs.area.to_numpy())
    if not np.isfinite(pixel_sizes).all() or (pixel_sizes <= 0).any():
        raise ValueError('Positive finite pixel-to-micron calibration is required.')
    pixel_size = float(np.median(pixel_sizes))
    np.testing.assert_allclose(pixel_sizes, pixel_size, rtol=1e-5)
    xy = np.asarray(adata.obsm["spatial"], dtype=float) * pixel_size
    window = np.asarray(window_pixels, dtype=float) * pixel_size
    if (not np.isfinite(xy).all() or (xy[:, 0] < window[0]).any() or (xy[:, 0] > window[1]).any()
            or (xy[:, 1] < window[2]).any() or (xy[:, 1] > window[3]).any()):
        raise ValueError("Coordinates fall outside the recorded analysis window.")
    t_mask = adata.obs.resolved_T_cell.to_numpy(bool).copy()
    m_mask = adata.obs.resolved_mMac_candidate.to_numpy(bool).copy()
    if (t_mask & m_mask).any():
        raise ValueError("Distinct parent lineages are required; no self-pairs.")
    t_xy, m_xy = xy[t_mask], xy[m_mask]
    codes = (adata.obs.resolved_CD8_T_cell.to_numpy(np.uint8)[t_mask]
             + 2 * adata.obs.resolved_GZB_T_cell.to_numpy(np.uint8)[t_mask])
    mac = adata.obs.resolved_CD274_mMac_candidate.to_numpy(bool)[m_mask]
    n_t_primary = int((codes == 3).sum())
    if n_t_primary == 0 or not mac.any():
        raise ValueError("Every tested population must contain cells.")
    if len(np.unique(xy[t_mask | m_mask], axis=0)) != (t_mask | m_mask).sum():
        raise ValueError("Duplicate centroids require inspection before distance analysis.")
    distances = cdist(t_xy, m_xy)
    ti, mi = np.where(distances <= radii[-1])
    bins = np.searchsorted(radii, distances[ti, mi], side="left").astype(np.int32)
    width_um, height_um = window[1]-window[0], window[3]-window[2]
    area = width_um * height_um
    dx, dy = np.abs(t_xy[ti]-m_xy[mi]).T
    weights = area / ((width_um-dx) * (height_um-dy))
    return dict(radii=radii, xy=xy, t_xy=t_xy, m_xy=m_xy, codes=codes, mac=mac,
                n_t_primary=n_t_primary, n_mac=int(mac.sum()), window=window, area=area,
                pixel_size=pixel_size, distances=distances, ti=ti.astype(np.int32),
                mi=mi.astype(np.int32), bins=bins, weights=weights,
                t_names=np.asarray(adata.obs_names[t_mask], dtype=str),
                m_names=np.asarray(adata.obs_names[m_mask], dtype=str))

def permute_labels(values, groups, *, shuffle, n_permutations, rng):
    result = np.tile(values, (n_permutations+1, 1))
    if shuffle:
        for group in np.unique(groups):
            indices = np.flatnonzero(groups == group)
            if len(indices) < 2 or np.unique(values[indices]).size < 2:
                continue
            for i in range(1, n_permutations+1):
                result[i, indices] = rng.permutation(values[indices])
    # Audit every code count within every exchangeability group, including joint T labels.
    for group in np.unique(groups):
        indices = np.flatnonzero(groups == group)
        for code in np.unique(values):
            np.testing.assert_array_equal((result[:, indices] == code).sum(axis=1),
                                           np.full(n_permutations+1, (values[indices] == code).sum()))
    return np.ascontiguousarray(result)

@njit(cache=True)
def evaluate_permutations(t_codes, m_positive, ti, mi, edge_bins, weights, n_radii):
    """Evaluate only the CD8A+/Granzyme B+ T versus CD274+ mMac pair."""
    n_trials = t_codes.shape[0]
    weighted = np.zeros((n_trials, n_radii))
    pairs = np.zeros((n_trials, n_radii), dtype=np.int32)
    for trial in range(n_trials):
        for edge in range(len(ti)):
            if m_positive[trial, mi[edge]] and t_codes[trial, ti[edge]] == 3:
                weighted[trial, edge_bins[edge]] += weights[edge]
                pairs[trial, edge_bins[edge]] += 1
        for radius in range(1, n_radii):
            weighted[trial, radius] += weighted[trial, radius-1]
            pairs[trial, radius] += pairs[trial, radius-1]
    return weighted, pairs

def upper_global_test(curves, expected):
    """One-sided global max-T test; pooling is symmetric over all curves.

    The analytic null mean depends only on fixed group counts and coordinates.
    The pooled SD includes the observed curve, so all B+1 curves are treated
    exchangeably. The Monte Carlo p value counts ties conservatively.
    """
    sd = curves.std(axis=0, ddof=1)
    active = sd > 1e-12
    deviations = np.zeros_like(curves)
    deviations[:, active] = (curves[:, active]-expected[active]) / sd[active]
    max_t = deviations.max(axis=1)
    p = float(np.count_nonzero(max_t >= max_t[0]) / len(max_t))
    return dict(p=p, observed_stat=float(max_t[0]),
                standardized_observed=deviations[0], test_radii=active)

def analytic_expectations(data):
    """Exact primary-pair expectation for the single within-lineage shuffle."""
    probability = np.full(len(data["ti"]),
                          np.mean(data["codes"] == 3) * np.mean(data["mac"]))
    raw = np.bincount(data["bins"], weights=probability, minlength=len(data["radii"])).cumsum()
    weighted = np.bincount(data["bins"], weights=probability*data["weights"],
                           minlength=len(data["radii"])).cumsum()
    expected_k = weighted * data["area"] / (data["n_t_primary"]*data["n_mac"])
    return expected_k, raw

def run_scheme(adata, scheme, *, window_pixels, n_permutations=9999, seed=20260930, radii_um=None, return_details=False):
    """Run the requested whole-ROI within-lineage shuffle for the primary pair only."""
    if scheme != SCHEMES[0]:
        raise ValueError("Only the global_both shuffle is used in Figure S11.")
    if n_permutations < 2:
        raise ValueError("At least two permutations are required.")
    data = prepare_inputs(adata, window_pixels=window_pixels, radii_um=radii_um)
    t_groups = np.zeros(len(data["t_xy"]), dtype=int)
    m_groups = np.zeros(len(data["m_xy"]), dtype=int)
    rng = np.random.default_rng(seed)
    t_codes = permute_labels(data["codes"], t_groups, shuffle=True, n_permutations=n_permutations, rng=rng)
    m_positive = permute_labels(data["mac"], m_groups, shuffle=True, n_permutations=n_permutations, rng=rng)
    weighted, pairs = evaluate_permutations(t_codes, m_positive, data["ti"], data["mi"], data["bins"],
                                            data["weights"], len(data["radii"]))
    curves = weighted * data["area"] / (data["n_t_primary"]*data["n_mac"])
    expected, expected_pairs = analytic_expectations(data)
    # Check the implementation against exact expected counts under this shuffle.
    se = pairs[1:].std(axis=0, ddof=1) / np.sqrt(n_permutations)
    assert np.all(np.abs(pairs[1:].mean(axis=0)-expected_pairs) <= 6*se+1e-7)
    test = upper_global_test(curves, expected)
    global_table = pd.DataFrame([dict(
        scheme=scheme["name"], scheme_title=scheme["title"], T_subtype=PRIMARY_SUBTYPE,
        n_T_parent=len(data["codes"]), n_mMac_parent=len(data["mac"]), n_T_subtype=data["n_t_primary"],
        n_CD274_mMac=data["n_mac"], p_global=test["p"], max_T=test["observed_stat"],
        simulations=n_permutations, seed=seed,
    )])
    curve_rows = []
    for k, radius in enumerate(data["radii"]):
        curve_rows.append(dict(scheme=scheme["name"], T_subtype=PRIMARY_SUBTYPE, radius_um=radius,
                               observed_K=curves[0,k], expected_K=expected[k],
                               K_enrichment=curves[0,k]/expected[k] if expected[k]>0 else np.nan,
                               null_pointwise_lo=np.quantile(curves[1:,k], .025),
                               null_pointwise_hi=np.quantile(curves[1:,k], .975),
                               observed_pairs=int(pairs[0,k]), expected_pairs=expected_pairs[k]))
    tables = (global_table, pd.DataFrame(curve_rows))
    if return_details:
        details = dict(K=curves, expected_K=expected, radii_um=data['radii'],
                       analysis_window_um=data['window'].copy(), analysis_area_um2=data['area'],
                       seed=seed, n_permutations=n_permutations,
                       bit_generator=type(rng.bit_generator).__name__)
        return (*tables, details)
    return tables


# 5. Plot the results
# ----------------------------------------------------------------------------

def plot_gate_distributions(adata):
    settings = json.loads(adata.uns['thresholds_json'])['markers']
    eligible = adata.obs.analysis_eligible.to_numpy(bool)
    fig, axes = plt.subplots(1, len(GATED_MARKERS), figsize=(15, 3), layout='constrained')
    for axis, marker in zip(axes, GATED_MARKERS):
        values = adata.layers['zscore_c5_ddof1'][eligible, adata.var_names.get_loc(marker)]
        setting = settings[marker]
        axis.hist(values, bins=60, color='0.65', linewidth=0)
        axis.axvline(setting['z_cutoff'], color='black', ls='--')
        axis.set(title=f"{marker}\n{setting['operator']} {setting['z_cutoff']:g}", xlabel='z-score', ylabel='Cells')
        axis.spines[['top', 'right']].set_visible(False)
    fig.suptitle(adata.uns['sample_id'])
    return fig

def plot_cell_populations(adata):
    fig, axes = plt.subplots(1, 2, figsize=(10, 5), layout='constrained')
    xy = adata.obsm['spatial_um']
    eligible = adata.obs.analysis_eligible.to_numpy(bool)
    for axis in axes:
        axis.scatter(xy[eligible, 0], xy[eligible, 1], s=1, c='0.85', rasterized=True)
    groups = [('resolved_T_cell', 'T cells', '#3575ab'), ('resolved_mMac_candidate', 'mMac', '#d47924')]
    targets = [('resolved_CD8_GZB_T_cell', 'CD8A+/GZMB+ T', '#3575ab'),
               ('resolved_CD274_mMac_candidate', 'CD274+ mMac', '#d47924')]
    for axis, states in zip(axes, (groups, targets)):
        for key, label, color in states:
            mask = adata.obs[key].to_numpy(bool)
            axis.scatter(xy[mask, 0], xy[mask, 1], s=5, c=color, label=f'{label} (n={mask.sum():,})', rasterized=True)
        window = np.asarray(adata.uns['analysis_window_pixels']) * adata.uns['scale_um_per_pixel']
        axis.set(xlim=window[:2], ylim=(window[3], window[2]), xlabel='x (µm)', ylabel='y (µm)', aspect='equal')
        axis.legend(frameon=False, fontsize=8)
    fig.suptitle(adata.uns['sample_id'])
    return fig

def plot_crossk(results):
    radius = results['radius_stats']
    core = results['core_stats'].iloc[0]
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), layout='constrained')
    axes[0].fill_between(radius.radius_um, radius.null_pointwise_lo, radius.null_pointwise_hi,
                         color='0.85', label='95% pointwise null range')
    axes[0].plot(radius.radius_um, radius.expected_K, color='0.4', ls='--', label='Analytic expectation')
    axes[0].plot(radius.radius_um, radius.observed_K, color='#3575ab', label='Observed')
    axes[0].set(xlabel='Radius (µm)', ylabel='Cross-K (µm²)', title=f'Global p = {core.p_global:.4g}')
    axes[0].legend(frameon=False, fontsize=8)
    axes[1].plot(radius.radius_um, radius.K_enrichment, color='#3575ab')
    axes[1].axhline(1, color='0.4', ls='--')
    axes[1].set(xlabel='Radius (µm)', ylabel='Observed / expected K', title='Cumulative enrichment')
    for axis in axes:
        axis.spines[['top', 'right']].set_visible(False)
    fig.suptitle(core.sample_id)
    return fig


# 6. Save individual and combined outputs
# ----------------------------------------------------------------------------

def save_results(adata, results, *, run_id, figures=None):
    """Save one final cell object and complete cross-K source data; never overwrite."""
    segmentation.require_compute()
    folder = segmentation.result_directory(run_id, adata.uns['sample_id'])
    folder.mkdir(parents=True, exist_ok=False)
    details = results['details']
    record = dict(sample_id=adata.uns['sample_id'], inputs=json.loads(adata.uns['inputs_json']),
                  classification=json.loads(adata.uns['classification_json']),
                  acquisitions=json.loads(adata.uns['acquisitions_json']),
                  seed=details['seed'], n_permutations=details['n_permutations'],
                  radii_um=details['radii_um'], bit_generator=details['bit_generator'],
                  analysis_window_um=details['analysis_window_um'],
                  cell_order='saved segmentation inventory order, filtered by analysis_eligible; no sorting',
                  permutation_order='all joint CD8A/GZMB T draws, then all CD274 mMac draws, one generator',
                  execution=segmentation.execution_record(), code=segmentation.helper_fingerprints(),
                  historical_equivalence_verified=False)
    record = json_value(record)
    adata.uns['run_provenance_json'] = json.dumps(record, sort_keys=True, allow_nan=False)
    adata.write_h5ad(folder / 'final_cells.h5ad', compression='gzip')
    results['core_stats'].to_csv(folder / 'core_stats.csv', index=False)
    results['radius_stats'].to_csv(folder / 'radius_stats.csv', index=False)
    np.savez_compressed(folder / 'crossk_curves.npz', observed_K=details['K'][0],
                        permutation_K=details['K'][1:], expected_K=details['expected_K'],
                        radii_um=details['radii_um'], ordered_cell_ids=details['ordered_cell_ids'])
    adata.obs.loc[details['ordered_cell_ids'],
                   ['resolved_lineage', 'resolved_CD8_T_cell', 'resolved_GZB_T_cell',
                    'resolved_CD274_mMac_candidate']].to_csv(folder / 'ordered_analysis_cells.csv', index_label='cell_id')
    segmentation.write_json(folder / 'run_settings.json', record)
    segmentation.write_environment_record(folder)
    for name, figure in (figures or {}).items():
        if not name or Path(name).name != name:
            raise ValueError('Use a simple figure basename.')
        figure.savefig(folder / f'{name}.pdf', bbox_inches='tight')
        figure.savefig(folder / f'{name}.png', dpi=300, bbox_inches='tight')
    return folder

def save_cohort_results(core_tables, radius_tables, *, output_folder, sample_ids):
    """Combine per-sample results without pooling cells or running extra tests."""
    output_folder = segmentation.require_orcd_path(output_folder)
    core = pd.concat(core_tables, ignore_index=True)
    radii = pd.concat(radius_tables, ignore_index=True)
    if core.sample_id.tolist() != sample_ids:
        raise ValueError('Combined results do not match the requested sample order.')
    if radii.sample_id.drop_duplicates().tolist() != sample_ids:
        raise ValueError('Radius results are missing a sample or changed sample order.')
    core_path = output_folder / 'all_samples_core_stats.csv'
    radii_path = output_folder / 'all_samples_radius_stats.csv'
    if core_path.exists() or radii_path.exists():
        raise FileExistsError('Combined outputs already exist; use a new run ID.')
    core.to_csv(core_path, index=False)
    radii.to_csv(radii_path, index=False)
    return core
