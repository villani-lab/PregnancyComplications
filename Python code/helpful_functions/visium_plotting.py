
from pathlib import Path
from typing import Optional, Sequence, Tuple, Union

import matplotlib as mpl
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import squidpy as sq
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar


_PUBLICATION_FONT_FAMILY = "Arial"
_PUBLICATION_FONT_PATHS = {
    "regular": None,
    "bold": None,
    "italic": None,
    "bold_italic": None,
}


def configure_publication_fonts(
    regular_font_path="/home/plopez/software/fonts/Arial.ttf",
    bold_font_path="/home/plopez/software/fonts/Arial-Bold.ttf",
    italic_font_path="/home/plopez/software/fonts/Arial-Italic.ttf",
    bold_italic_font_path="/home/plopez/software/fonts/Arial-BoldItalic.ttf",
    font_family="Arial",
):
    """Register publication fonts and update Matplotlib defaults."""
    global _PUBLICATION_FONT_FAMILY, _PUBLICATION_FONT_PATHS

    font_paths = {
        "regular": regular_font_path,
        "bold": bold_font_path,
        "italic": italic_font_path,
        "bold_italic": bold_italic_font_path,
    }

    for key, font_path in font_paths.items():
        if font_path and Path(font_path).exists():
            fm.fontManager.addfont(str(font_path))
            font_paths[key] = str(font_path)
        else:
            font_paths[key] = None

    _PUBLICATION_FONT_FAMILY = font_family
    _PUBLICATION_FONT_PATHS = font_paths

    mpl.rcParams.update(
        {
            "font.family": font_family,
            "font.sans-serif": [font_family],
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


configure_publication_fonts()


def get_publication_fontproperties(style="normal", weight="normal", size=None):
    if style == "italic" and weight == "bold":
        font_path = _PUBLICATION_FONT_PATHS.get("bold_italic")
    elif style == "italic":
        font_path = _PUBLICATION_FONT_PATHS.get("italic")
    elif weight == "bold":
        font_path = _PUBLICATION_FONT_PATHS.get("bold")
    else:
        font_path = _PUBLICATION_FONT_PATHS.get("regular")

    if font_path is not None:
        return fm.FontProperties(fname=font_path, size=size)

    fallback_path = _PUBLICATION_FONT_PATHS.get("regular") or _PUBLICATION_FONT_PATHS.get("bold")
    if fallback_path is not None:
        return fm.FontProperties(fname=fallback_path, size=size)

    raise FileNotFoundError(
        "Publication font file not found for the requested style. "
        "Run configure_publication_fonts(...) with valid Arial font file paths."
    )


def get_spatial_library(adata):
    return next(iter(adata.uns["spatial"]))


def get_microns_per_pixel(adata):
    spatial_lib = get_spatial_library(adata)
    return adata.uns["spatial"][spatial_lib]["scalefactors"]["microns_per_pixel"]


def get_scale_bar_px(adata, scale_bar_um=100):
    return scale_bar_um / get_microns_per_pixel(adata)


def add_scale_bar(
    ax,
    scale_bar_px,
    *,
    label="",
    loc="lower right",
    pad=0.2,
    bbox_to_anchor=(1.0, -0.03),
    bbox_transform=None,
    size_vertical=0.7,
    fontproperties=None,
    zorder=None,
):
    if fontproperties is None:
        fontproperties = get_publication_fontproperties(weight="bold", size=12)

    scalebar = AnchoredSizeBar(
        ax.transData,
        scale_bar_px,
        label,
        loc=loc,
        pad=pad,
        bbox_to_anchor=bbox_to_anchor,
        bbox_transform=bbox_transform,
        color="black",
        frameon=False,
        size_vertical=size_vertical,
        fontproperties=fontproperties,
    )
    if zorder is not None:
        scalebar.set_zorder(zorder)
    ax.add_artist(scalebar)
    return scalebar


def add_rectangles(ax, rects, crop_coord, linewidth=2.5, linestyle="--", zorder=20):
    crop_x0, crop_y0, _, _ = crop_coord
    for (x0, y0, x1, y1), color in rects:
        ax.add_patch(
            Rectangle(
                (x0 - crop_x0, y0 - crop_y0),
                x1 - x0,
                y1 - y0,
                fill=False,
                edgecolor=color,
                linewidth=linewidth,
                linestyle=linestyle,
                zorder=zorder,
                clip_on=False,
            )
        )


def clean_spatial_axis(ax, title, *, fontsize=20, pad=4, fontweight="bold"):
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title(title, fontsize=fontsize, fontweight=fontweight, pad=pad)


def save_figure(fig, filename, output_dir, *, dpi=600, bbox_inches="tight", pad_inches=0.01, format=None):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    save_path = output_dir / filename
    fig.savefig(
        save_path,
        dpi=dpi,
        bbox_inches=bbox_inches,
        pad_inches=pad_inches,
        format=format,
    )
    return save_path


def spatial_scatter_with_coords_top_on_top(
    adata,
    color,
    img=True,
    library_key="library_id",
    cmap="YlOrRd",
    shape='hex',
    size=1,
    n_ticks=5,
    panel_width=6,
    panel_height=5,
    wspace=0.5,
    invert_y=True,
    alpha=1.0,           # dot transparency
    img_key=None,        # "hires"/"lowres" (auto if None)
    img_alpha=1.0,       # background image transparency
    show=True,
    **kwargs
):
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from scipy import sparse

    def _get_expr_vector(adata, color):
        if color in adata.obs.columns:
            return adata.obs[color].to_numpy(dtype=float)
        elif color in adata.var_names:
            X = adata[:, color].X
            if sparse.issparse(X):
                X = X.toarray()
            return np.asarray(X, dtype=float).ravel()
        raise KeyError(f"'{color}' not found in .obs or .var_names.")

    expr = _get_expr_vector(adata, color)
    lib_vals = adata.obs[library_key]
    lib_order = (lib_vals.cat.categories
                 if pd.api.types.is_categorical_dtype(lib_vals)
                 else pd.Index(pd.unique(lib_vals)))

    fig, axes = plt.subplots(
        1, len(lib_order),
        figsize=(panel_width * len(lib_order), panel_height),
        squeeze=False
    )
    axes = axes.flatten()

    marker = 'h' if str(shape).lower() == 'hex' else 'o'

    for ax, lib in zip(axes, lib_order):
        # ----- Optional tissue image -----
        if img:
            try:
                spatial_meta = adata.uns.get("spatial", {}).get(str(lib), {})
                images = spatial_meta.get("images", {})
                use_img_key = (img_key if (img_key and img_key in images)
                               else ("hires" if "hires" in images
                                     else (next(iter(images)) if images else None)))
                if use_img_key is not None:
                    im = images[use_img_key]
                    h, w = im.shape[:2]
                    sf = spatial_meta.get("scalefactors", {}).get(f"tissue_{use_img_key}_scalef", 1.0)
                    ax.imshow(
                        im, origin="upper", interpolation="nearest",
                        extent=(0, w / sf, h / sf, 0),
                        alpha=img_alpha, zorder=0
                    )
            except Exception:
                pass

        # ----- Spots: low→high so highs draw last -----
        mask = (lib_vals.values == lib)
        coords = adata.obsm['spatial'][mask, :]
        expr_lib = expr[mask]
        order = np.argsort(expr_lib)

        # Per-library bounds: [0, max]
        vmin_use = 0.0
        vmax_use = float(np.nanmax(expr_lib)) if expr_lib.size else 1.0

        sc = ax.scatter(
            coords[order, 0], coords[order, 1],
            c=expr_lib[order], cmap=cmap, s=size,
            edgecolor="none", alpha=alpha, marker=marker, zorder=1,
            vmin=vmin_use, vmax=vmax_use
        )

        # Limits + padding
        xmin, xmax = float(coords[:, 0].min()), float(coords[:, 0].max())
        ymin, ymax = float(coords[:, 1].min()), float(coords[:, 1].max())
        pad_frac = 0.06
        xr = xmax - xmin if xmax > xmin else 1.0
        yr = ymax - ymin if ymax > ymin else 1.0
        ax.set_xlim(xmin - xr * pad_frac, xmax + xr * pad_frac)
        ax.set_ylim(ymin - yr * pad_frac, ymax + yr * pad_frac)
        if invert_y:
            ax.invert_yaxis()
        ax.set_aspect('equal')

        # ticks/labels
        ax.set_xticks(np.linspace(xmin, xmax, n_ticks))
        ax.set_yticks(np.linspace(ymin, ymax, n_ticks))
        ax.tick_params(left=True, bottom=True, labelleft=True, labelbottom=True)
        ax.ticklabel_format(style='plain', useOffset=False)
        ax.set_xlabel("spatial1")
        ax.set_ylabel("spatial2")
        ax.set_title(f"{lib} — {color}")

        fig.colorbar(sc, ax=ax)

    fig.subplots_adjust(wspace=wspace)
    if show:
        plt.show()


def gene_vec(adata, gene):
    x = adata[:, gene].X
    if hasattr(x, "toarray"):   # sparse
        x = x.toarray()
    return np.asarray(x).ravel()

def plot_tma_crop_or_rect(
    adata,
    library_id: str,
    orientation: str,
    region: Optional[str] = None,         # e.g., "VC", "Decidua"
    color: Optional[Union[str, Sequence[str]]] = None,
    img: bool = True,
    groups=None,
    shape: str = "hex",
    size: float = 1.0,
    outline: bool = False,
    cmap: str = "YlOrRd",
    vmin:float = 0,
    vmax:float = 1,
    title: Optional[str] = None,
    mode: str = "crop",                   # "crop" or "rect"
    rect_color: str = "red",
    rect_lw: float = 2.0,
    zoom_rect: Optional[Tuple[int, int, int, int]] = None,   # (xmin, ymin, xmax, ymax)
    zoom_rect_is_global: bool = True,
    figsize=(4, 4),
    ax=None,                              # allow plotting into existing axis
    reset_categorical_colors: bool = True # NEW: fixes "all gray" categorical plots
):
    import pandas as pd
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    from matplotlib.patches import Rectangle
    import squidpy as sq

    # constructed tma table
    tma2_140_df = pd.DataFrame({
        'library_id': ["TMA2-140"]*4,
        'Patient': ['BMP-0122', 'BMP-0122', 'BMP-0064', 'BMP-0064'],
        'Region': ['Decidua', 'VC', 'Decidua', 'VC'],
        'Condition': ['PTB']*4,
        'xmin': [280, 280, 1000, 900],
        'xmax': [780, 780, 1500, 1400],
        'ymin': [900, 1530, 875, 1530],
        'ymax': [1400, 2030, 1375, 2030],
        'orientation': ['top_left','bottom_left','top_right','bottom_right']
    })

    tma2_145_df = pd.DataFrame({
        'library_id': ["TMA2-145"]*4,
        'Patient': ['BMP-0038', 'BMP-0038', 'BMP-0151', 'BMP-0151'],
        'Region': ['VC', 'Decidua', 'VC', 'Decidua'],
        'Condition': ['Cont_Unlbr']*4,
        'xmin': [350, 250, 1020, 950],
        'xmax': [850, 825, 1520, 1450],
        'ymin': [900, 1550, 975, 1550],
        'ymax': [1400, 2050, 1475, 2050],
        'orientation': ['top_left','bottom_left','top_right','bottom_right']
    })

    tma4_df = pd.DataFrame({
        'library_id': ["TMA4"]*4,
        'Patient': ['BMP-0023', 'BMP-0023', 'BMP-0077', 'BMP-0077'],
        'Region': ['VC', 'Decidua', 'VC', 'Decidua'],
        'Condition': ['Cont_Lbr']*4,
        'xmin': [425, 425, 1110, 1100],
        'xmax': [925, 925, 1610, 1610],
        'ymin': [955, 1650, 975, 1610],
        'ymax': [1455, 2150, 1475, 2110],
        'orientation': ['top_left','bottom_left','top_right','bottom_right']
    })

    tma9_df = pd.DataFrame({
        'library_id': ["TMA9"]*4,
        'Patient': ['BMP-0121', 'BMP-0121', 'BMP-0119', 'BMP-0119'],
        'Region': ['Decidua', 'VC', 'Decidua', 'VC'],
        'Condition': ['TermPE_FGR', 'TermPE_FGR', 'PrePE_FGR', 'PrePE_FGR'],
        'xmin': [520, 510, 1150, 1150],
        'xmax': [1020, 1020, 1650, 1650],
        'ymin': [900, 1500, 950, 1525],
        'ymax': [1400, 2025, 1450, 2025],
        'orientation': ['top_left','bottom_left','top_right','bottom_right']
    })

    coord_df = pd.concat([tma2_140_df, tma2_145_df, tma4_df, tma9_df], ignore_index=True)

    # --- lookup ---
    sub = coord_df[(coord_df["library_id"] == library_id) & (coord_df["orientation"] == orientation)]
    if region is not None and "Region" in coord_df.columns:
        sub = sub[sub["Region"] == region]

    if sub.shape[0] == 0:
        raise ValueError(f"No coordinates found for library_id={library_id}, orientation={orientation}, region={region}")
    if sub.shape[0] > 1:
        print("Warning: multiple rows matched; using the first.")

    row = sub.iloc[0]

    # IMPORTANT: squidpy crop_coord expects (x0, y0, x1, y1)
    crop = (int(row["xmin"]), int(row["ymin"]), int(row["xmax"]), int(row["ymax"]))

    if title is None:
        bits = [f"{library_id} | {orientation}"]
        for k in ["Region", "Condition", "Patient"]:
            if k in row.index and pd.notnull(row[k]):
                bits.append(str(row[k]))
        title = " | ".join(bits)

    # --- axis handling ---
    created_fig = False
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
        created_fig = True
    else:
        fig = ax.figure

    crop_coord = crop if mode == "crop" else None

    # -----------------------------
    # Categorical palette handling (FIXED for Squidpy)
    # -----------------------------
    palette = None
    color_arg = None if color is None else (list(color) if isinstance(color, (list, tuple)) else [color])

    if color_arg is not None and len(color_arg) == 1:
        key = str(color_arg[0])
        if key in adata.obs:
            # if it's "object" but looks categorical, convert to categorical
            if adata.obs[key].dtype.name == "object":
                adata.obs[key] = adata.obs[key].astype("category")

            if adata.obs[key].dtype.name == "category":
                # wipe old/bad palette if requested
                if reset_categorical_colors and f"{key}_colors" in adata.uns:
                    del adata.uns[f"{key}_colors"]

                # ensure we have enough colors
                from scanpy.plotting.palettes import default_102
                ncat = adata.obs[key].cat.categories.size

                # store Scanpy-style list in .uns
                colors = default_102[:ncat]
                adata.uns[f"{key}_colors"] = colors

                # Squidpy expects string or ListedColormap (NOT a python list)
                palette = mpl.colors.ListedColormap(colors)

    # Safety: if some other code path produced a list palette, convert it
    if isinstance(palette, list):
        palette = mpl.colors.ListedColormap(palette)

    # --- plot ---
    sq.pl.spatial_scatter(
        adata,
        color=color_arg,
        img=img,
        library_key="library_id",
        library_id=library_id,
        shape=shape,
        size=size,
        use_raw=False,
        outline=outline,
        vmin = vmin,
        vmax = vmax,
        groups=groups,
        cmap=cmap,
        palette=palette,
        title=title,
        crop_coord=crop_coord,
        ax=ax
    )

    # draw crop rectangle on full view
    if mode == "rect":
        x0, y0, x1, y1 = crop
        ax.add_patch(Rectangle((x0, y0), x1-x0, y1-y0,
                               linewidth=rect_lw, edgecolor=rect_color,
                               facecolor="none", zorder=100))

    # optional zoom rectangle
    if zoom_rect is not None:
        zx0, zy0, zx1, zy1 = zoom_rect

        if zoom_rect_is_global and mode == "crop":
            cx0, cy0, _, _ = crop
            zx0, zx1 = zx0 - cx0, zx1 - cx0
            zy0, zy1 = zy0 - cy0, zy1 - cy0

        ax.add_patch(Rectangle((zx0, zy0), zx1-zx0, zy1-zy0,
                               linewidth=rect_lw, edgecolor=rect_color,
                               facecolor="none", zorder=101))

    if created_fig:
        plt.tight_layout()

    return fig, ax, crop


# Computes the mean z-score per cell of a set of genes
# Scale each gene in the gene set across all cells in the data, then computes an average of z-scores of all the genes. Saves the score in the obs under <score_title> columns
def signature_score_per_cell(data, gene_set, score_title) :
    # Get rid of genes that aren't in data
    orig_len = len(gene_set)
    gene_set = [gene for gene in gene_set if gene in data.var_names]
    print(str(len(gene_set)) + "/" + str(orig_len)  + " of the gene set genes are measured" )

    # Limit the data to just those genes
    dat = data[:,gene_set].X

    dat = dat.toarray()
    mean = dat.mean(axis=0)
    var = dat.var(axis=0)
    std = np.sqrt(var)
    std[std==0]=0.0001

    with np.errstate(divide="ignore", invalid="ignore"):
        dat = (dat - mean) / std
    dat[dat < -5] = -5
    dat[dat > 5] = 5

    scores = dat.mean(axis = 1)
    data.obs[score_title] = scores


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from statsmodels.nonparametric.smoothers_lowess import lowess

def _get_gene_expr(adata, gene, layer=None, use_raw=True):
    """
    Returns a 1D numpy array of expression for `gene`.
    - If layer is provided, uses adata.layers[layer]
    - Else if use_raw and adata.raw exists, uses adata.raw
    - Else uses adata.X
    """
    if layer is not None:
        X = adata[:, gene].layers[layer]
    else:
        if use_raw and (adata.raw is not None):
            X = adata.raw[:, gene].X
        else:
            X = adata[:, gene].X

    # convert sparse -> dense if needed
    try:
        import scipy.sparse as sp
        if sp.issparse(X):
            X = X.toarray()
    except Exception:
        pass

    X = np.asarray(X).reshape(-1)
    return X

def plot_obs_vs_gene_overlay(
    adata_map,
    obs_col="HBC",
    gene="SPP1",
    layer=None,          # e.g., "norm_counts" if you want
    use_raw=True,
    frac=0.25,           # LOWESS smoothing span
    scatter_alpha=0.25,  # transparency for dots
    s=18,                # dot size
):
    """
    adata_map: dict like {
        "Preterm Birth": PTB_only,
        "Control (Labored)": Ctrl_lbr,
        "Control (Unlabored)": Ctrl_unlbr
    }
    """
    rows = []
    for label, ad in adata_map.items():
        x = ad.obs[obs_col].to_numpy()
        y = _get_gene_expr(ad, gene, layer=layer, use_raw=use_raw)
        rows.append(pd.DataFrame({"x": x, "y": y, "Condition": label}))

    df = pd.concat(rows, ignore_index=True)

    # Choose distinct colors (can modulate)
    color_map = {
        "Preterm Birth": "#D55E00",
        "Control (Labored)": "#0072B2",
        "Control (Unlabored)": "#009E73",
    }

    fig, ax = plt.subplots(figsize=(7.5, 5.5))

    # scatter + smooth per condition
    for cond, g in df.groupby("Condition", sort=False):
        c = color_map.get(cond, None)

        ax.scatter(
            g["x"], g["y"],
            s=s,
            alpha=scatter_alpha,
            label=cond,
            rasterized=True,
            color=c
        )

        # LOWESS smooth line
        # (sort x to make a nice curve)
        g2 = g.sort_values("x")
        sm = lowess(g2["y"], g2["x"], frac=frac, it=0, return_sorted=True)
        ax.plot(sm[:, 0], sm[:, 1], linewidth=2.5, color=c)

    ax.set_xlabel(obs_col)
    ax.set_ylabel(f"{gene} expression")
    ax.set_title(f"{obs_col} vs {gene} (overlay by condition)")
    ax.legend(frameon=False)
    plt.tight_layout()
    return fig, ax, df

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf


def prep_bin_covariates_filtered(
        adata,
        signature_col: str,
        condition_col: str = "Condition",
        patient_col: str = "Patient",
        villus_col: str = "which_VC",
        villus_exclude_value: str = "Not_in_VC",
        villi_type_col: str = "Type_of_Villi",
        villi_type_keep: str = "Terminal",
        hbc_col: str = "ab_parent__HBC",
        total_abund_col: str = "total_cell_abundance",
        conditions_keep=None,
        condition_order=None,
        eps: float = 1e-6,
):
    obs = adata.obs.copy()

    needed = [
        signature_col, condition_col, patient_col, villus_col,
        villi_type_col, hbc_col, total_abund_col
    ]
    missing = [c for c in needed if c not in obs.columns]
    if missing:
        raise ValueError(f"Missing columns in adata.obs: {missing}")

    # Ensure strings
    obs[villus_col] = obs[villus_col].astype(str)
    obs[villi_type_col] = obs[villi_type_col].astype(str)
    obs[condition_col] = obs[condition_col].astype(str)
    obs[patient_col] = obs[patient_col].astype(str)

    # 1) Filter: must be in a villus + Terminal villi only
    obs = obs[
        (obs[villus_col] != str(villus_exclude_value)) &
        (obs[villi_type_col] == str(villi_type_keep))
        ].copy()

    # 2) Filter: keep only selected Conditions
    if conditions_keep is not None:
        obs = obs[obs[condition_col].isin(list(conditions_keep))].copy()

    # Build covariates
    obs["HBC_frac"] = obs[hbc_col] / (obs[total_abund_col] + eps)
    obs["HBC_frac"] = obs["HBC_frac"].clip(lower=0)
    obs["log_total_cells"] = np.log1p(obs[total_abund_col])

    # Set categorical order for Condition (controls the reference)
    if condition_order is None:
        # If not provided, use sorted unique in data
        condition_order = sorted(obs[condition_col].unique().tolist())

    obs[condition_col] = pd.Categorical(obs[condition_col], categories=condition_order, ordered=True)

    # Drop missing model inputs
    model_cols = [signature_col, condition_col, patient_col, villus_col, "HBC_frac", "log_total_cells"]
    obs = obs.dropna(subset=model_cols).copy()

    if obs.shape[0] == 0:
        raise ValueError("No bins left after filtering. Check labels/filters.")

    return obs

import pandas as pd
import statsmodels.formula.api as smf

def _ensure_patient_villus(obs: pd.DataFrame, patient_col="Patient", villus_col="which_VC") -> pd.DataFrame:
    obs = obs.copy()
    if "patient_villus" not in obs.columns:
        if patient_col not in obs.columns or villus_col not in obs.columns:
            raise ValueError(f"Need columns '{patient_col}' and '{villus_col}' to build patient_villus.")
        obs["patient_villus"] = obs[patient_col].astype(str) + ":" + obs[villus_col].astype(str)
    return obs


def fit_signature_mixedlm_FULL(
        obs: pd.DataFrame,
        signature_col: str,
        condition_col: str = "Condition",
        patient_col: str = "Patient",
        villus_col: str = "which_VC",   # kept for compatibility
        reml: bool = True,
):
    obs = _ensure_patient_villus(obs, patient_col=patient_col, villus_col=villus_col)

    formula = f"{signature_col} ~ C({condition_col}) * HBC_frac"
    vc = {"villus_in_patient": "0 + C(patient_villus)"}

    m = smf.mixedlm(
        formula=formula,
        data=obs,
        groups=obs[patient_col],
        vc_formula=vc,
        re_formula="1",
    )
    r = m.fit(reml=reml, method="lbfgs", maxiter=300, disp=False)
    return r


def fit_signature_mixedlm_REDUCED(
        obs: pd.DataFrame,
        signature_col: str,
        condition_col: str = "Condition",
        patient_col: str = "Patient",
        villus_col: str = "which_VC",   # kept for compatibility
        reml: bool = True,
):
    obs = _ensure_patient_villus(obs, patient_col=patient_col, villus_col=villus_col)

    formula = f"{signature_col} ~ HBC_frac + C({condition_col})"
    vc = {"villus_in_patient": "0 + C(patient_villus)"}

    m = smf.mixedlm(
        formula=formula,
        data=obs,
        groups=obs[patient_col],
        vc_formula=vc,
        re_formula="1",
    )
    r = m.fit(reml=reml, method="lbfgs", maxiter=300, disp=False)
    return r

import numpy as np
import pandas as pd
import scipy.sparse as sp

def pathway_filter_table(
    adata,
    gene_set,
    min_frac=0.05,
    min_mean=0.05
):
    rows = []

    for gene in gene_set:

        if gene not in adata.var_names:
            rows.append({
                "gene": gene,
                "mean_expr": np.nan,
                "frac_nonzero": np.nan,
                "keep": False,
                "reason": "not_in_var_names"
            })
            continue

        X = adata[:, gene].X
        if sp.issparse(X):
            X = X.toarray().ravel()
        else:
            X = np.asarray(X).ravel()

        mean_expr = float(X.mean())
        frac_nonzero = float((X > 0).mean())

        reasons = []
        if frac_nonzero < min_frac:
            reasons.append(f"low_frac<{min_frac}")
        if mean_expr < min_mean:
            reasons.append(f"low_mean<{min_mean}")

        keep = len(reasons) == 0

        rows.append({
            "gene": gene,
            "mean_expr": mean_expr,
            "frac_nonzero": frac_nonzero,
            "keep": keep,
            "reason": " | ".join(reasons) if not keep else "kept"
        })

    df = pd.DataFrame(rows)
    return df

# Code for plotting marker genes
import os
import matplotlib.pyplot as plt
import squidpy as sq

def plot_spatial_lineage_markers(
    adata,
    objs_plot,
    celltype_markers_dict,
    crop_coord,
    include_cell2location=True,
    marker_columns=None,
    markers_by="row",
    layer=None,
    use_raw=False,
    cell2loc_suffix="_prop",
    library_key="library_id",
    library_id="TMA2-140",
    lineage_title="PTB",
    lineage_col_title="cell2location\nresults",
    marker_title_prefix="Marker",
    cmap="YlOrRd",
    shape="square",
    size=1,
    outline=False,
    img_for_genes=False,
    figsize_per_panel=(2.55, 2.45),
    wspace=0.02,
    hspace=0.06,
    row_label_x=-0.16,
    col_title_y=1.12,
    row_label_fontsize=16,
    col_title_fontsize=14,
    panel_title_fontsize=11,
    title_overrides=None,
    lineage_title_overrides=None,
    add_scale_bars=True,
    scale_bar_um=100,
    scale_bar_kwargs=None,
    italicize_gene_titles=True,
    save_dir=None,
    save_name=None,
    dpi=600,
):
    if title_overrides is None:
        title_overrides = {}
    if lineage_title_overrides is None:
        lineage_title_overrides = {}
    if scale_bar_kwargs is None:
        scale_bar_kwargs = {}
    scale_bar_kwargs = dict(scale_bar_kwargs)

    if marker_columns is None:
        marker_columns = max(len(celltype_markers_dict.get(ct, [])) for ct in objs_plot)

    if markers_by not in {"row", "column"}:
        raise ValueError("markers_by must be either 'row' or 'column'")

    if markers_by == "row":
        nrows = len(objs_plot)
        ncols = marker_columns + (1 if include_cell2location else 0)
    else:
        nrows = marker_columns + (1 if include_cell2location else 0)
        ncols = len(objs_plot)

    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(figsize_per_panel[0] * ncols, figsize_per_panel[1] * nrows),
        squeeze=False,
        gridspec_kw={"wspace": wspace, "hspace": hspace},
    )
    fig.patch.set_facecolor("white")

    if markers_by == "row":
        col_titles = []
        if include_cell2location:
            col_titles.append(lineage_col_title)
        col_titles.extend([f"{marker_title_prefix} {idx + 1}" for idx in range(marker_columns)])
        row_labels = [title_overrides.get(obj, obj) for obj in objs_plot]
    else:
        col_titles = [title_overrides.get(obj, obj) for obj in objs_plot]
        row_labels = []
        if include_cell2location:
            row_labels.append(lineage_col_title)
        row_labels.extend([f"{marker_title_prefix} {idx + 1}" for idx in range(marker_columns)])

    row_to_genes = {
        ct: (list(celltype_markers_dict.get(ct, [])) + [None] * marker_columns)[:marker_columns]
        for ct in objs_plot
    }

    if markers_by == "row":
        panel_grid = []
        for ct in objs_plot:
            panels = []
            if include_cell2location:
                panels.append(
                    {
                        "feature": f"{ct}{cell2loc_suffix}",
                        "panel_title": lineage_title_overrides.get(ct, lineage_title),
                        "img": True,
                        "is_cell2location": True,
                    }
                )
            panels.extend(
                [
                    {
                        "feature": gene,
                        "panel_title": title_overrides.get(gene, gene) if gene is not None else "",
                        "img": img_for_genes,
                        "is_cell2location": False,
                    }
                    for gene in row_to_genes[ct]
                ]
            )
            panel_grid.append(panels)
    else:
        panel_grid = []
        for row_idx in range(nrows):
            panels = []
            for ct in objs_plot:
                if include_cell2location and row_idx == 0:
                    panels.append(
                        {
                            "feature": f"{ct}{cell2loc_suffix}",
                            "panel_title": lineage_title_overrides.get(ct, lineage_title),
                            "img": True,
                            "is_cell2location": True,
                        }
                    )
                else:
                    gene_idx = row_idx - 1 if include_cell2location else row_idx
                    gene = row_to_genes[ct][gene_idx]
                    panels.append(
                        {
                            "feature": gene,
                            "panel_title": title_overrides.get(gene, gene) if gene is not None else "",
                            "img": img_for_genes,
                            "is_cell2location": False,
                        }
                    )
            panel_grid.append(panels)

    figure_titles = []
    figure_col_titles = []
    figure_row_labels = []

    for r in range(nrows):
        for c in range(ncols):
            panel = panel_grid[r][c]
            ax = axes[r, c]
            feature = panel["feature"]

            if feature is None:
                ax.axis("off")
                continue

            scatter_kwargs = dict(
                adata=adata,
                color=feature,
                img=panel["img"],
                library_key=library_key,
                library_id=library_id,
                cmap=cmap,
                shape=shape,
                size=size,
                title=None,
                outline=outline,
                crop_coord=crop_coord,
                ax=ax,
            )

            if panel["is_cell2location"]:
                sq.pl.spatial_scatter(**scatter_kwargs)
            else:
                if layer is not None:
                    scatter_kwargs["layer"] = layer
                    scatter_kwargs["use_raw"] = False
                else:
                    scatter_kwargs["use_raw"] = use_raw
                sq.pl.spatial_scatter(**scatter_kwargs)

            if add_scale_bars:
                scale_bar_px = get_scale_bar_px(
                    adata[adata.obs[library_key] == library_id].copy(),
                    scale_bar_um=scale_bar_um,
                )
                scale_bar_kwargs.setdefault("bbox_transform", ax.transAxes)
                add_scale_bar(ax, scale_bar_px, zorder=30, **scale_bar_kwargs)

            ax.set_xlabel("")
            ax.set_ylabel("")
            ax.set_xticks([])
            ax.set_yticks([])
            for spine in ax.spines.values():
                spine.set_visible(True)
                spine.set_color("black")
                spine.set_linewidth(1)

            if r == 0:
                figure_col_titles.append((ax, col_titles[c]))

            ax.set_title("")

            if panel["panel_title"]:
                figure_titles.append(
                    (
                        ax,
                        panel["panel_title"],
                        panel["is_cell2location"],
                    )
                )

    for r, row_label in enumerate(row_labels):
        figure_row_labels.append((axes[r, 0], row_label))

    fig.canvas.draw()

    for ax, title, is_cell2location in figure_titles:
        bbox = ax.get_position()
        title_kwargs = dict(
            x=bbox.x0 + bbox.width / 2,
            y=bbox.y1 + 0.01,
            ha="center",
            va="bottom",
            fontsize=panel_title_fontsize,
        )
        if is_cell2location or not italicize_gene_titles:
            fig.text(s=title, **title_kwargs)
        else:
            fig.text(
                s=title,
                fontproperties=get_publication_fontproperties(
                    style="italic",
                    size=panel_title_fontsize,
                ),
                **title_kwargs,
            )

    for ax, title in figure_col_titles:
        bbox = ax.get_position()
        fig.text(
            bbox.x0 + bbox.width / 2,
            bbox.y1 + 0.06,
            title,
            ha="center",
            va="bottom",
            fontsize=col_title_fontsize,
            fontweight="bold",
        )

    for ax, row_label in figure_row_labels:
        bbox = ax.get_position()
        fig.text(
            bbox.x0 + row_label_x * bbox.width,
            bbox.y0 + bbox.height / 2,
            row_label,
            ha="right",
            va="center",
            fontsize=row_label_fontsize,
            fontweight="bold",
        )

    if save_dir is not None:
        os.makedirs(save_dir, exist_ok=True)

        if save_name is None:
            save_name = f"spatial_lineage_markers_{library_id}.pdf"
        elif not save_name.lower().endswith(".pdf"):
            save_name = f"{save_name}.pdf"

        save_path = os.path.join(save_dir, save_name)
        fig.savefig(save_path, format="pdf", dpi=dpi, bbox_inches="tight")
        print(f"Saved figure to: {save_path}")

    return fig, axes
