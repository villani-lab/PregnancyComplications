"""Figure S11 image access, DAPI segmentation and geometry helpers.

"""

from collections import Counter
from datetime import datetime, timezone
from importlib import metadata
from math import ceil, floor
from pathlib import Path
import hashlib
import inspect
import json
import os
import re
import subprocess
import sys
import xml.etree.ElementTree as ET

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
from shapely.geometry import LinearRing, MultiPolygon, Polygon, box
from shapely.errors import GEOSException
from skimage.draw import polygon as polygon_pixels
from skimage.measure import find_contours, regionprops, regionprops_table
from skimage.segmentation import expand_labels
import tifffile

PROJECTS = Path('/orcd/pool/007/lauffen_shared/Projects')
REFERENCE = PROJECTS / 'Pregnancy_complications_Miltenyi'
CODE_ROOT = REFERENCE / 'Code_PAL'
RESULTS = PROJECTS / 'Preg_proteomics_review/Results/0_SupFig8_colocalization'
MARKERS = ('DAPI', 'CD3', 'CD8A', 'Granzyme_B', 'CD163', 'CD68', 'CD274')
FALLBACK_PIXEL_SIZE_UM = 0.325


# 1. Read registered images and select marker channels
# ----------------------------------------------------------------------------

def load_segmentation_dapi(image_path):
    """Read the FIRST DAPI acquisition, matching KEEP_FIRST_DAPI_ONLY.

    The first stored acquisition is selected independently of protein channels.
    This preserves the cycle-0 DAPI acquisition for TMA7 ROI7.
    """
    channel_table, shape = inspect_image_channels(image_path)
    manifest = pd.read_csv(Path(image_path).with_suffix('.channels.csv'))
    if (manifest.channel_0_based.tolist() != channel_table.channel_index.tolist()
            or manifest.marker.tolist() != channel_table.image_channel.tolist()):
        raise ValueError('DAPI channel manifest disagrees with the registered image.')
    dapi = manifest.loc[manifest.marker.astype(str).str.contains('dapi', case=False)]
    if dapi.empty:
        raise ValueError('No DAPI channel in the registered image.')
    row = dapi.iloc[0]
    index = int(row.channel_0_based)
    with tifffile.TiffFile(image_path) as tif:
        plane = tif.series[0].levels[0].pages[index].asarray()
    _validate_plane(plane)
    if tuple(plane.shape) != tuple(shape):
        raise ValueError('DAPI plane dimensions disagree with its image.')
    plane.setflags(write=False)
    return plane, dict(channel_index=index, cycle=int(row.cycle), name=str(row.marker),
                       source_file=str(row.source_file),
                       rule='first DAPI in stored order; KEEP_FIRST_DAPI_ONLY=true')

def resolve_registered_image(project_root, sample_name, image_path=None):
    """Resolve one registered full-resolution image; never choose among runs."""
    if image_path is not None:
        path = Path(image_path).expanduser().resolve()
        if not path.is_file():
            raise FileNotFoundError(path)
        return path
    match = re.fullmatch(r"(TMA[0-9]+)_ROI([0-9]+)", sample_name)
    if not match:
        raise ValueError(f"Unexpected sample name: {sample_name}")
    tma, roi = match.groups()
    root = Path(project_root) / "Data/2_ashlar_registered"
    paths = sorted({p.resolve() for folder in root.glob(f"Run*_{tma}")
                    for name in (f"ROI_{roi}", f"ROI{roi}")
                    for p in (folder / name).glob("*.ome.tif*") if p.is_file()})
    if len(paths) != 1:
        raise ValueError(f"Expected one registered image for {sample_name}; found {paths}. "
                         "Set threshold_image_path explicitly to the image used for segmentation.")
    return paths[0]

def inspect_image_channels(image_path, remap_log=None):
    """Read OME channel names and the actual export remapping, without pixels.

    Supported: a full-resolution OME series with one YX page per channel and
    singleton time/Z axes. Unsupported layouts stop rather than guessing axes.
    """
    with tifffile.TiffFile(image_path) as tif:
        series = tif.series[0]
        dims = dict(zip(series.axes, series.shape))
        if (not series.axes.endswith("YX") or not {"Y", "X"}.issubset(dims)
            or any(n != 1 for axis, n in dims.items() if axis not in "CYX")):
            raise ValueError(f"Expected 2D channels, received {series.axes}, {series.shape}.")
        count = dims.get("C", 1)
        if len(series.pages) != count:
            raise ValueError("This TIFF is not stored as one YX page per channel.")
        if not tif.ome_metadata:
            raise ValueError("OME channel metadata are required; refusing positional guesses.")
        pixels = ET.fromstring(tif.ome_metadata).find("{*}Image/{*}Pixels")
        channels = [] if pixels is None else pixels.findall("{*}Channel")
        names = [c.get("Name", "").strip() for c in channels]
        if len(names) != count or not all(names):
            raise ValueError("OME metadata must name every channel.")
        shape = (dims["Y"], dims["X"])
    markers = [None] * count
    if remap_log is not None and Path(remap_log).is_file():
        log = pd.read_csv(remap_log, keep_default_na=False)
        if not {"original_name", "canonical_name"}.issubset(log.columns):
            raise ValueError("Remap log needs original_name and canonical_name columns.")
        original = log.original_name.astype(str).tolist()
        canonical = log.canonical_name.astype(str).tolist()
        if not all(original) or not all(canonical):
            raise ValueError("Empty channel names in remap log.")
        # Repeat the exporter suffix rule, but never use its deduplicated row
        # position as an index into the original (multi-DAPI) image.
        seen = Counter()
        used = set()
        for raw, name in zip(original, canonical):
            unique_name = name if seen[name] == 0 else f"{name}_{seen[name]}"
            seen[name] += 1
            if "dapi" in name.casefold():
                continue
            candidates = [i for i, value in enumerate(names) if value == raw and i not in used]
            if names.count(raw) != original.count(raw) or not candidates:
                raise ValueError(f"Cannot map {raw!r} from remap log to OME channels exactly.")
            index = candidates[0]  # Same-named acquisitions retain their recorded order.
            markers[index] = unique_name
            used.add(index)
    else:
        counts = Counter(names)
        markers = [name if counts[name] == 1 else None for name in names]
    return pd.DataFrame({"channel_index": range(count), "image_channel": names,
                         "table_marker": markers}), shape

def checked_bounds(bounds, shape):
    if bounds is None:
        return [0, shape[1], 0, shape[0]]
    if len(bounds) != 4 or any(isinstance(v, bool) or int(v) != v for v in bounds):
        raise ValueError("Bounds must be integer (x0, x1, y0, y1) pixel coordinates.")
    x0, x1, y0, y1 = map(int, bounds)
    if not (0 <= x0 < x1 <= shape[1] and 0 <= y0 < y1 <= shape[0]):
        raise ValueError(f"Bounds {bounds} are outside image shape {shape}.")
    return [x0, x1, y0, y1]

CHANNEL_NAME_ALIASES = {
    ('cd45', 18): ('CD45', 'CD45_CST_AF647'),
    ('cd45_2b11', 24): ('2B11_PD726', 'CD45_2B11', 'CD45_2B11_PD7_26'),
}

def _marker_key(name):
    """Allow spelling differences in separators/case, without dropping clone names."""
    return re.sub(r"[\s_-]+", "_", str(name).strip().casefold())

def load_raw_markers(image_path, markers, channels_csv=None, *, require_cycles=False):
    """Read selected level-0 YX planes using verified OME/manifest channel names.

    ``markers`` maps output labels to a cycle (or None), or to a specification
    such as {"names": ["CD4"], "cycle": 10}. Distinct labels preserve repeated
    acquisitions, e.g. CD4_cycle10 and CD4_cycle21. Explicit names can include
    known export aliases; multiple matches are rejected, never picked by order.
    Confirmed CD45 aliases are also recognized here for older notebook settings,
    but only for their recorded cycles and a requested name in the same group.

    A manifest must describe every channel in stored order and agree with OME
    names. ``require_cycles=True`` requires this manifest for verified cycle
    selection. Otherwise, unique OME names are sufficient, with unknown cycle.
    Panel row order is never used as TIFF channel order.
    """
    image_path = Path(image_path).expanduser().resolve()
    manifest_path = (Path(channels_csv).expanduser() if channels_csv is not None
                     else image_path.with_suffix(".channels.csv"))
    if channels_csv is not None and not manifest_path.is_file():
        raise FileNotFoundError(manifest_path)
    if require_cycles and not manifest_path.is_file():
        raise FileNotFoundError(
            f"Cycle-specific selection requires the matching channel manifest: {manifest_path}. "
            "Set CHANNELS_CSV explicitly if it is stored elsewhere."
        )
    if not markers:
        raise ValueError("Provide at least one marker.")
    with tifffile.TiffFile(image_path) as tif:
        series = tif.series[0].levels[0]
        dims = dict(zip(series.axes, series.shape))
        if (not series.axes.endswith("YX") or not {"Y", "X"}.issubset(dims)
                or any(n != 1 for axis, n in dims.items() if axis not in "CYX")):
            raise ValueError(f"Expected 2D channels, got {series.axes} {series.shape}.")
        count = dims.get("C", 1)
        if len(series.pages) != count or not tif.ome_metadata:
            raise ValueError("Expected OME metadata and one full-resolution YX page per channel.")
        pixels = ET.fromstring(tif.ome_metadata).find("{*}Image/{*}Pixels")
        names = [] if pixels is None else [
            ch.get("Name", "").strip() for ch in pixels.findall("{*}Channel")
        ]
        if len(names) != count:
            raise ValueError("OME channel count does not match the TIFF.")
        if manifest_path.is_file():
            channels = pd.read_csv(manifest_path, keep_default_na=False)
            required = {"channel_0_based", "marker", "cycle"}
            if not required.issubset(channels):
                raise ValueError(f"Channel manifest requires {sorted(required)}.")
            if channels.channel_0_based.tolist() != list(range(count)):
                raise ValueError("Channel manifest must list every TIFF channel in stored order.")
            channels["marker"] = channels.marker.astype(str).str.strip()
            cycles = pd.to_numeric(channels.cycle, errors="raise")
            if not np.isfinite(cycles).all() or (cycles % 1 != 0).any():
                raise ValueError("Manifest cycles must be finite integers.")
            channels["cycle"] = cycles.astype(int)
            if not channels.marker.str.len().gt(0).all():
                raise ValueError("Channel manifest contains an empty marker name.")
            if any(name and name != channels.marker.iloc[i] for i, name in enumerate(names)):
                raise ValueError("Channel manifest and OME names disagree.")
            mapping_source = str(manifest_path.resolve())
        else:
            if not all(names):
                raise ValueError("Unnamed OME channels: supply the matching .channels.csv; "
                                 "do not infer channel indices from the antibody panel.")
            channels = pd.DataFrame({"channel_0_based": range(count), "marker": names,
                                     "cycle": [None] * count})
            mapping_source = "OME channel names (cycle unavailable)"

        raw_images, rows = {}, []
        for marker, specification in markers.items():
            if isinstance(specification, dict):
                if set(specification) - {"names", "cycle"}:
                    raise ValueError("Marker specifications accept only names and cycle.")
                expected_cycle = specification.get("cycle")
                aliases = specification.get("names", [marker])
                aliases = [aliases] if isinstance(aliases, str) else list(aliases)
            else:
                expected_cycle, aliases = specification, [marker]
            if not aliases or any(not isinstance(name, str) or not name.strip() for name in aliases):
                raise ValueError(f"Provide nonempty channel names for {marker}.")
            if expected_cycle is not None and (
                    isinstance(expected_cycle, bool)
                    or not isinstance(expected_cycle, (int, np.integer)) or expected_cycle < 1):
                raise ValueError(f"Expected cycle for {marker} must be a positive integer.")
            if require_cycles and expected_cycle is None:
                raise ValueError(f"Provide an expected cycle for {marker}.")
            alias_keys = {_marker_key(a) for a in aliases}
            known_aliases = CHANNEL_NAME_ALIASES.get((_marker_key(marker), expected_cycle), ())
            known_keys = {_marker_key(a) for a in known_aliases}
            if manifest_path.is_file() and alias_keys & known_keys:
                alias_keys.update(known_keys)
            match = channels[channels.marker.map(_marker_key).isin(alias_keys)]
            if expected_cycle is not None and manifest_path.is_file():
                match = match[match.cycle == expected_cycle]
            if len(match) != 1:
                available = channels
                if expected_cycle is not None and manifest_path.is_file():
                    available = channels.loc[channels.cycle == expected_cycle]
                listed = available[['channel_0_based', 'marker', 'cycle']].to_dict('records')
                other_matches = channels.loc[channels.marker.map(_marker_key).isin(alias_keys),
                                             ['channel_0_based', 'marker', 'cycle']].to_dict('records')
                recorded_cycles = sorted(channels.cycle.dropna().unique().tolist())
                raise ValueError(
                    f"No unique channel for {marker}, cycle {expected_cycle}: {len(match)} matches. "
                    f"Requested names: {aliases}; accepted names: {sorted(alias_keys)}. "
                    f"Available channels for this cycle: {listed}. "
                    f"Requested-name matches across all cycles: {other_matches}. "
                    f"Cycles recorded in this image: {recorded_cycles}. "
                    f"Image: {image_path}. Channel mapping: {mapping_source}."
                )
            row = match.iloc[0]
            index = int(row.channel_0_based)
            if index in [r["channel_0_based"] for r in rows]:
                raise ValueError("Requested markers map to the same channel.")
            plane = series.pages[index].asarray()
            _validate_plane(plane)
            if plane.shape != (dims["Y"], dims["X"]):
                raise ValueError("Channel page does not match the OME image dimensions.")
            # Prevent accidental in-place modification in later notebook cells.
            plane.setflags(write=False)
            raw_images[marker] = plane
            rows.append({"marker": marker, "image_channel": row.marker, "cycle": row.cycle,
                         "channel_0_based": index, "shape_yx": plane.shape,
                         "dtype": str(plane.dtype), "raw_min": float(plane.min()),
                         "raw_max": float(plane.max()), "image_path": str(image_path),
                         "mapping_source": mapping_source, "pyramid_level": 0})
    return raw_images, pd.DataFrame(rows).set_index("marker")

def _validate_plane(raw):
    if (not isinstance(raw, np.ndarray) or raw.ndim != 2 or raw.size == 0
            or not np.issubdtype(raw.dtype, np.number) or np.iscomplexobj(raw)
            or not np.isfinite(raw).all() or np.any(raw < 0)):
        raise ValueError("Expected a nonempty, finite, nonnegative 2D intensity array.")


# 2. Convert segmentation masks to measurement geometry
# ----------------------------------------------------------------------------

def native_pixel_regions(polygons, shape):
    """Rasterize each polygon once in its true image coordinates.

    Unlike the historical chunk implementation, the raster origin is the actual
    clipped view origin. Pixels corrected to zero remain in the denominator.
    This function preserves the source geometry and raw measurements.
    """
    regions = []
    counts = np.zeros(len(polygons), dtype=np.int64)
    height, width = shape
    for index, cell in enumerate(polygons.geometry):
        bx, by, ex, ey = cell.bounds
        x0, y0 = (max(0, floor(bx)), max(0, floor(by)))
        x1, y1 = (min(width, ceil(ex) + 1), min(height, ceil(ey) + 1))
        if x1 <= x0 or y1 <= y0:
            raise ValueError(f'Cell {index} is outside the current image.')
        mask = np.zeros((y1 - y0, x1 - x0), dtype=bool)
        for part in cell.geoms if isinstance(cell, MultiPolygon) else [cell]:
            vertices = np.asarray(part.exterior.coords)
            rr, cc = polygon_pixels(vertices[:, 1] - y0, vertices[:, 0] - x0, mask.shape)
            mask[rr, cc] = True
            for ring in part.interiors:
                vertices = np.asarray(ring.coords)
                rr, cc = polygon_pixels(vertices[:, 1] - y0, vertices[:, 0] - x0, mask.shape)
                mask[rr, cc] = False
        counts[index] = int(mask.sum())
        if not counts[index]:
            raise ValueError(f'Cell {index} contains no image pixels.')
        regions.append((index, slice(y0, y1), slice(x0, x1), mask))
    return (regions, counts)

def _mask_polygons(mask_path):
    """Reconstruct the complete expanded-mask geometry before matching final cells."""
    try:
        labels = tifffile.memmap(mask_path, mode='r')
    except ValueError:
        labels = tifffile.imread(mask_path)
    if (labels.ndim != 2 or not np.issubdtype(labels.dtype, np.integer)
            or labels.min() < 0 or labels.max() > np.iinfo(np.int32).max):
        raise ValueError('Expected a nonnegative 2D integer segmentation mask.')
    cells = []
    for region in regionprops(np.asarray(labels, dtype=np.int32)):
        rings = [LinearRing(c[:, ::-1]) for c in find_contours(np.pad(region.image, 1), .5)
                 if len(c) >= 4]
        parts = [Polygon(r) for r in rings if r.is_ccw]
        if not parts:
            raise ValueError(f'No cell boundary for mask label {region.label}.')
        offset = np.array([region.bbox[1]-1, region.bbox[0]-1])
        parts = [Polygon(np.asarray(p.exterior.coords)+offset) for p in parts]
        cells.append(parts[0] if len(parts) == 1 else MultiPolygon(parts))
    if not cells:
        raise ValueError('The segmentation mask contains no cells.')
    mean_radius = float(np.mean(np.sqrt(np.array([p.area for p in cells])/np.pi)))
    radius, tolerance = .1*mean_radius, min(.03*mean_radius, 1.)
    smoothed = []
    for cell in cells:
        try:
            cell = cell.buffer(-radius).buffer(2*radius).buffer(-radius).simplify(tolerance)
        except GEOSException:
            # Preserve the original contour-cleanup behavior.
            continue
        cell = shapely.make_valid(cell)
        if isinstance(cell, Polygon):
            cell = Polygon(cell.exterior) if cell.interiors else cell
        else:
            parts = [p for p in getattr(cell, 'geoms', ()) if isinstance(p, Polygon)]
            cell = max(parts, key=lambda p: p.area) if parts else Polygon()
        if not cell.is_empty:
            smoothed.append(cell)
    return gpd.GeoDataFrame(geometry=smoothed)



# 3. Resolve sample paths and record runs
# ----------------------------------------------------------------------------

def require_orcd_path(path):
    """Reject runtime paths outside the analysis project root."""
    path = Path(path)
    if not path.is_absolute() or '..' in path.parts or not path.is_relative_to(PROJECTS):
        raise ValueError(f'Use an absolute path inside {PROJECTS}: {path}')
    return path


def require_compute():
    """Data processing belongs in an allocated ORCD compute session."""
    if not os.environ.get('SLURM_JOB_ID') or not PROJECTS.is_dir():
        raise RuntimeError('Run this workflow in an allocated ORCD compute session.')


def sample_identity(tma, roi, row):
    if any(isinstance(v, bool) or int(v) != v or v < 1 for v in (tma, roi, row)):
        raise ValueError('TMA, ROI and ROW must be positive integers.')
    return f'TMA{int(tma)}_ROI{int(roi)}_Row{int(row)}'


def result_directory(run_id, sample_id):
    if not re.fullmatch(r'[A-Za-z0-9_-]+', run_id):
        raise ValueError('Use a simple, new run identifier.')
    if not re.fullmatch(r'TMA\d+_ROI\d+_Row\d+', sample_id):
        raise ValueError('Invalid sample identifier.')
    return RESULTS / run_id / sample_id


def selected_image(tma, roi, *, image_relative=None):
    """Use the established sample resolver; overrides select one acquisition."""
    sample = f'TMA{int(tma)}_ROI{int(roi)}'
    override = None
    if image_relative is not None:
        relative = Path(image_relative)
        if relative.is_absolute() or '..' in relative.parts:
            raise ValueError('An acquisition override must be project-relative.')
        override = require_orcd_path(REFERENCE / relative)
        run_parts = [part for part in relative.parts if re.fullmatch(r'Run\w+_TMA\d+', part)]
        if (not any(part.endswith(f'_TMA{int(tma)}') for part in run_parts)
                or not any(part in (f'ROI_{int(roi)}', f'ROI{int(roi)}') for part in relative.parts)):
            raise ValueError('The image override must belong to the selected TMA/ROI.')
    return resolve_registered_image(REFERENCE, sample, image_path=override)


def marker_acquisitions(sample_id):
    """Recorded S11 cycles, including the Pt_11 exceptions."""
    specs = dict(CD3=dict(names=['CD3'], cycle=13),
                 CD8A=dict(names=['CD8a'], cycle=24),
                 Granzyme_B=dict(names=['Granzyme_B'], cycle=17),
                 CD163=dict(names=['CD163'], cycle=12),
                 CD68=dict(names=['CD68'], cycle=25),
                 CD274=dict(names=['CD274'], cycle=19))
    if sample_id == 'TMA5_ROI20_Row1':
        specs['CD8A']['cycle'], specs['CD68']['cycle'] = 25, 26
    return specs


def recorded_crop(sample_id):
    return [0, 4250, 200, 2550] if sample_id == 'TMA10_ROI4_Row6' else None


def write_json(path, record):
    require_orcd_path(path).write_text(json.dumps(record, indent=2, allow_nan=False) + '\n')


def execution_record():
    """Record the environment actually used, without claiming a historical match."""
    versions = {}
    for name in ('numpy', 'scipy', 'numba', 'pandas', 'anndata', 'geopandas',
                 'shapely', 'scikit-image', 'tifffile', 'imagecodecs',
                 'matplotlib', 'pyarrow', 'cellpose', 'torch'):
        try:
            versions[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            versions[name] = None
    commit = subprocess.run(['git', '-C', str(CODE_ROOT), 'rev-parse', 'HEAD'],
                            text=True, capture_output=True, check=False)
    return dict(python=sys.version, dependencies=versions,
                git_commit=commit.stdout.strip() if commit.returncode == 0 else None,
                slurm_job_id=os.environ.get('SLURM_JOB_ID'),
                generated_utc=datetime.now(timezone.utc).isoformat())


def helper_fingerprints():
    files = tuple(Path(__file__).with_name(name) for name in
                  ('run_figs11.py', 'figs11_segmentation_helpers.py', 'figs11_analysis_helpers.py'))
    return {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in files}


def write_environment_record(folder):
    """Save installed package versions; this is not an OS/CUDA build lock."""
    packages = {dist.metadata['Name']: dist.version for dist in metadata.distributions()
                if dist.metadata.get('Name')}
    lines = ['# Actual Python: ' + sys.version.replace('\n', ' '),
             '# Installed Python package versions; OS/CUDA builds are not locked.']
    lines.extend(f'{name}=={version}' for name, version in sorted(packages.items()))
    require_orcd_path(folder / 'environment.lock.txt').write_text('\n'.join(lines) + '\n')



# 4. Segment DAPI and expand nuclear labels
# ----------------------------------------------------------------------------

def get_pixel_size_um(path, fallback=None):
    """PhysicalSizeX from OME-XML converted to microns; else `fallback`."""
    fallback = FALLBACK_PIXEL_SIZE_UM if fallback is None else fallback
    try:
        with tifffile.TiffFile(path) as tif:
            if tif.ome_metadata:
                ns = {"ome": "http://www.openmicroscopy.org/Schemas/OME/2016-06"}
                px = ET.fromstring(tif.ome_metadata).find(".//ome:Pixels", ns)
                if px is not None and px.attrib.get("PhysicalSizeX"):
                    val = float(px.attrib["PhysicalSizeX"])
                    unit = px.attrib.get("PhysicalSizeXUnit", "µm")
                    val *= {"nm": 1e-3, "mm": 1e3, "m": 1e6}.get(unit, 1.0)
                    print(f"Pixel size from OME metadata: {val:.4f} µm/px")
                    return val
    except Exception as e:  # noqa: BLE001
        print(f"Could not read pixel size from metadata ({e}).")
    print(f"Using fallback pixel size: {fallback} µm/px")
    return fallback

def segment_dapi(dapi, *, expected_cellpose_version, fresh_tma8=False):
    """Preserve the original call structure; the TMA8 rerun has its own profile."""
    from cellpose import models
    import torch

    version = metadata.version('cellpose')
    if version != expected_cellpose_version:
        raise RuntimeError(f'Cellpose {version} is installed; requested {expected_cellpose_version}.')
    if fresh_tma8 and version != '4.0.6':
        raise ValueError('The documented TMA8 ROI7 rerun used Cellpose 4.0.6.')
    if not torch.cuda.is_available():
        raise RuntimeError('The requested segmentation GPU is unavailable.')
    lo, hi = np.percentile(dapi, [1, 99.5])
    if hi <= lo:
        raise ValueError('DAPI normalization range is empty.')
    normalized = np.clip((dapi.astype(np.float32) - lo) / (hi - lo), 0, 1)
    kwargs = dict(flow_threshold=.4, cellprob_threshold=0., min_size=1, max_size_fraction=.4)
    if fresh_tma8:
        np.random.seed(0)
        torch.manual_seed(0)
        model = models.CellposeModel(gpu=True, pretrained_model='cpsam')
        kwargs.update(normalize=True, batch_size=8)
    else:
        model = models.CellposeModel(gpu=True)
    supported = inspect.signature(model.eval).parameters
    effective = {key: value for key, value in kwargs.items() if key in supported}
    # Do not silently drop scientifically relevant settings when an API differs.
    if set(effective) != set(kwargs):
        raise TypeError(f'Cellpose does not accept the recorded parameters: {set(kwargs)-set(effective)}')
    nuclei = np.asarray(model.eval(normalized, **effective)[0], dtype=np.int32)
    if nuclei.shape != dapi.shape or nuclei.max() == 0:
        raise ValueError('Segmentation produced no nuclei or an unexpected image shape.')
    record = dict(cellpose_version=version, model_kwargs=effective,
                  model_selection='cpsam' if fresh_tma8 else 'original CellposeModel default',
                  resolved_model=str(getattr(model, 'pretrained_model', 'unavailable')),
                  normalization_percentiles=[1, 99.5], normalization_limits=[float(lo), float(hi)],
                  seed=0 if fresh_tma8 else None)
    return nuclei, record


def segmentation_table(nuclei, expanded, geometry, scale, sample_id, bounds):
    """Retain all detections; express final inclusion as flags, not deletion."""
    table = pd.DataFrame(regionprops_table(nuclei, properties=('label', 'area', 'centroid')))
    table = table.rename(columns={'label': 'nucleus_label', 'area': 'nuclear_area_px',
                                  'centroid-0': 'nuclear_centroid_y', 'centroid-1': 'nuclear_centroid_x'})
    suffix = 'redo_nucleus' if sample_id == 'TMA8_ROI7_Row5' else 'mask'
    table['analysis_cell_id'] = [f'{sample_id}_{suffix}_{int(i):06d}' for i in table.nucleus_label]
    table.index = pd.Index(table.analysis_cell_id, name='cell_id')
    table['source_cell_id'] = table.analysis_cell_id
    table['nuclear_area_um2'] = table.nuclear_area_px * scale**2
    for prefix, mask in [('nucleus', nuclei), ('expansion', expanded)]:
        edge_ids = np.unique(np.concatenate([mask[0], mask[-1], mask[:, 0], mask[:, -1]]))
        table[f'{prefix}_touches_image_edge'] = table.nucleus_label.isin(edge_ids[edge_ids > 0])
    table['polygon_available'] = table.nucleus_label.isin(geometry.nucleus_label)
    paired = geometry.set_index('nucleus_label').geometry.reindex(table.nucleus_label)
    paired.index = table.index
    table['polygon_valid'] = paired.is_valid & ~paired.is_empty & paired.notna()
    rectangle = box(bounds[0], bounds[2], bounds[1], bounds[3])
    inside, crop_edge = [], []
    for polygon in paired:
        usable = polygon is not None and not polygon.is_empty
        inside.append(usable and rectangle.covers(polygon))
        crop_edge.append(usable and polygon.intersects(rectangle.boundary))
    table['within_saved_analysis_window'] = inside
    if recorded_crop(sample_id) is None:
        table['within_saved_analysis_window'] = True
    table['expansion_touches_saved_crop_edge'] = crop_edge if recorded_crop(sample_id) is not None else False
    table['below_minimum_nuclear_area'] = table.nuclear_area_um2 < 3.
    table['exclusion_reason'] = np.select(
        [~table.within_saved_analysis_window, ~table.polygon_available | ~table.polygon_valid,
         table.nucleus_touches_image_edge, table.expansion_touches_image_edge,
         table.expansion_touches_saved_crop_edge, table.below_minimum_nuclear_area],
        ['outside_saved_analysis_window', 'unusable_expanded_geometry', 'nucleus_at_image_edge',
         'expansion_only_at_image_edge', 'expansion_at_saved_crop_edge', 'tiny_nuclear_fragment'], default='')
    table['analysis_eligible'] = table.exclusion_reason == ''
    paired_geometry = gpd.GeoDataFrame({'analysis_cell_id': table.analysis_cell_id,
                                      'nucleus_label': table.nucleus_label}, geometry=paired)
    return table, paired_geometry



# 5. Save a sample's segmentation and raw measurements
# ----------------------------------------------------------------------------

def run_segmentation(tma, roi, row, *, run_id, pixel_size_um, expected_cellpose_version,
                     image_relative=None):
    """Segment one registered sample and write only reusable segmentation outputs."""
    import anndata as ad

    require_compute()
    sample = sample_identity(tma, roi, row)
    if not np.isfinite(pixel_size_um) or pixel_size_um <= 0:
        raise ValueError('Supply the saved positive pixel calibration in micrometres.')
    image = selected_image(tma, roi, image_relative=image_relative)
    dapi, selection = load_segmentation_dapi(image)
    expected_dapi_cycle = 0 if sample == 'TMA7_ROI7_Row5' else 1
    if selection['cycle'] != expected_dapi_cycle:
        raise ValueError('The first DAPI cycle differs from the recorded manuscript acquisition.')
    bounds = checked_bounds(recorded_crop(sample), dapi.shape)
    folder = result_directory(run_id, sample)
    if folder.exists():
        raise FileExistsError(f'Choose a new run ID; existing output is preserved: {folder}')
    nuclei, model_record = segment_dapi(dapi, expected_cellpose_version=expected_cellpose_version,
                                       fresh_tma8=sample == 'TMA8_ROI7_Row5')
    expansion_pixels = int(round(5. / pixel_size_um))
    expanded = expand_labels(nuclei, distance=expansion_pixels).astype(np.int32)
    folder.mkdir(parents=True, exist_ok=False)
    tifffile.imwrite(folder / 'nuclear_mask.tif', nuclei, compression='zlib')
    tifffile.imwrite(folder / 'expanded_mask.tif', expanded, compression='zlib')
    geometry = _mask_polygons(folder / 'expanded_mask.tif')
    points = geometry.geometry.representative_point()
    labels = expanded[np.clip(points.y.to_numpy().astype(int), 0, dapi.shape[0]-1),
                      np.clip(points.x.to_numpy().astype(int), 0, dapi.shape[1]-1)]
    if (labels == 0).any() or len(np.unique(labels)) != len(labels):
        raise ValueError('Polygons do not map one-to-one to detected nuclei.')
    geometry['nucleus_label'] = labels
    table, geometry = segmentation_table(nuclei, expanded, geometry, pixel_size_um, sample, bounds)
    valid = table.polygon_available & table.polygon_valid
    regions, counts = native_pixel_regions(geometry.loc[valid], dapi.shape)
    channels, _ = inspect_image_channels(image)
    raw = np.full((len(table), len(channels)), np.nan, dtype=np.float64)
    with tifffile.TiffFile(image) as tif:
        for channel in channels.channel_index:
            plane = tif.series[0].levels[0].pages[int(channel)].asarray()
            means = np.zeros(len(counts), dtype=np.float64)
            for index, ys, xs, mask in regions:
                means[index] = plane[ys, xs][mask].astype(np.float64).sum() / counts[index]
            raw[valid.to_numpy(), int(channel)] = means
    channels.index = pd.Index([f'image_channel_{int(i):03d}' for i in channels.channel_index], name='channel')
    channels['table_marker'] = channels.table_marker.fillna('')
    source = ad.AnnData(X=raw, obs=table.copy(), var=channels)
    source.layers['raw'] = raw.copy()
    source.uns['sample_id'] = sample
    source.write_h5ad(folder / 'raw_cells.h5ad', compression='gzip')
    table.to_csv(folder / 'all_detections_with_eligibility.csv', index=False)
    geometry.to_parquet(folder / 'all_expanded_geometry_with_eligibility.parquet')
    record = dict(schema='FigS11_segmentation_v1', sample_id=sample, registered_image=str(image),
                  source_table=str(folder / 'raw_cells.h5ad'), scale_um_per_pixel=float(pixel_size_um),
                  image_shape_yx=list(dapi.shape), analysis_bounds_pixels=bounds,
                  saved_crop=None if recorded_crop(sample) is None else
                  dict(bounds_x0_x1_y0_y1=bounds, sample_id=sample, source_image_shape_yx=list(dapi.shape)),
                  nuclear_mask=str(folder / 'nuclear_mask.tif'), expanded_mask=str(folder / 'expanded_mask.tif'),
                  detected=len(table), eligible=int(table.analysis_eligible.sum()),
                  dapi_selection=selection, expansion_um=5., expansion_pixels=expansion_pixels,
                  model=model_record, execution=execution_record(), code=helper_fingerprints(),
                  identity_note='Fresh detections; original manuscript cell identities are not transferred.',
                  historical_equivalence_verified=False)
    write_json(folder / 'provenance.json', record)
    write_environment_record(folder)
    return folder
