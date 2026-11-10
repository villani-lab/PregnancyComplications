#%%
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib as mpl
import squidpy as sq
import os


def plot_cellint_vs_genescore(
    adata,
    cell_int,
    gene_score,
    gene_score_color,
    cell_int_color=None,
    dominant_type_key="dominant_parent_type",
    dominant_type_colors_key="dominant_parent_type_colors",
    scale_quantile=0.99,
    alpha=0.9,
    size=35,
    figsize=(6, 5),
    title=None,
    show_legends=True,
    show_mix_legend=True,
    mix_legend_label=None,
    save_to_adata=False,
    save_prefix=None,
    save_pdf=False,
    save_dir="misc_files",
    save_name=None,
    dpi = 300
):
    categories = pd.Categorical(adata.obs[dominant_type_key]).categories
    category_to_color = dict(zip(categories, adata.uns[dominant_type_colors_key]))

    if cell_int_color is None:
        if cell_int not in category_to_color:
            raise ValueError(
                f"{cell_int} not found in {dominant_type_key}; pass cell_int_color manually."
            )
        cell_int_color = category_to_color[cell_int]

    x = adata.obs[cell_int].to_numpy(dtype=float)
    y = adata.obs[gene_score].to_numpy(dtype=float)

    x_max = np.quantile(x[np.isfinite(x)], scale_quantile)
    y_pos = y[np.isfinite(y) & (y > 0)]
    y_max = np.quantile(y_pos, scale_quantile) if len(y_pos) > 0 else 1.0

    if x_max <= 0:
        x_max = 1.0
    if y_max <= 0:
        y_max = 1.0

    x_scaled = np.clip(np.nan_to_num(x, nan=0.0) / x_max, 0, 1)
    y_scaled = np.clip(np.nan_to_num(np.maximum(y, 0), nan=0.0) / y_max, 0, 1)

    cell_rgb = np.array(mpl.colors.to_rgb(cell_int_color))
    gene_rgb = np.array(mpl.colors.to_rgb(gene_score_color))
    white_rgb = np.array([1.0, 1.0, 1.0])

    weights = np.column_stack([x_scaled**2, y_scaled**2])
    weight_sum = weights.sum(axis=1, keepdims=True)
    safe_weight_sum = np.where(weight_sum == 0, 1, weight_sum)
    mixed_hue = (
        weights[:, [0]] * cell_rgb +
        weights[:, [1]] * gene_rgb
    ) / safe_weight_sum
    intensity = np.maximum(x_scaled, y_scaled)[:, None]
    point_rgb = white_rgb * (1 - intensity) + mixed_hue * intensity
    point_rgb[weight_sum.ravel() == 0] = white_rgb

    mix_intensity = intensity.ravel()
    mix_legend_color = mpl.colors.to_hex((cell_rgb + gene_rgb) / 2)

    if save_to_adata:
        if save_prefix is None:
            save_prefix = f"{gene_score}_{cell_int}s"
        adata.obs[save_prefix] = mix_intensity
        adata.obsm[f"{save_prefix}_rgb"] = point_rgb.copy()
        adata.uns[f"{save_prefix}_color_params"] = {
            "cell_int": cell_int,
            "gene_score": gene_score,
            "cell_int_color": cell_int_color,
            "gene_score_color": gene_score_color,
            "mix_legend_label": mix_legend_label if mix_legend_label is not None else save_prefix,
            "scale_quantile": scale_quantile,
            "cell_int_vmax": float(x_max),
            "gene_score_vmax": float(y_max),
            "mix_legend_color": mix_legend_color,
        }

    fig, ax = plt.subplots(figsize=figsize)
    ax.scatter(x, y, c=point_rgb, alpha=alpha, s=size)
    ax.set_xlabel(cell_int)
    ax.set_ylabel(gene_score)
    plot_title = title if title is not None else f"{cell_int} vs {gene_score}"
    ax.set_title(plot_title)

    if show_legends:
        fig.subplots_adjust(right=0.76)
        cell_cmap = mpl.colors.LinearSegmentedColormap.from_list(
            f"{cell_int}_white_to_color",
            ["#FFFFFF", cell_int_color],
        )
        gene_cmap = mpl.colors.LinearSegmentedColormap.from_list(
            f"{gene_score}_white_to_color",
            ["#FFFFFF", gene_score_color],
        )
        mix_cmap = mpl.colors.LinearSegmentedColormap.from_list(
            f"{cell_int}_{gene_score}_mix",
            ["#FFFFFF", mix_legend_color],
        )

        cax1 = fig.add_axes([0.88, 0.68, 0.03, 0.18])
        cax2 = fig.add_axes([0.88, 0.40, 0.03, 0.18])
        cax3 = fig.add_axes([0.88, 0.12, 0.03, 0.18]) if show_mix_legend else None

        sm1 = mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=0, vmax=x_max), cmap=cell_cmap)
        sm2 = mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=0, vmax=y_max), cmap=gene_cmap)

        cb1 = fig.colorbar(sm1, cax=cax1)
        cb1.ax.set_title(cell_int, pad=6)
        cb1.set_ticks([0, x_max])
        cb1.set_ticklabels(["0", f"{x_max:.2g}"])

        cb2 = fig.colorbar(sm2, cax=cax2)
        cb2.ax.set_title(gene_score, pad=6)
        cb2.set_ticks([0, y_max])
        cb2.set_ticklabels(["0", f"{y_max:.2g}"])

        if show_mix_legend:
            mix_label = mix_legend_label if mix_legend_label is not None else f"{gene_score}_{cell_int}s"
            sm3 = mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=0, vmax=1), cmap=mix_cmap)
            cb3 = fig.colorbar(sm3, cax=cax3)
            cb3.ax.set_title(mix_label, pad=6)
            cb3.set_ticks([0, 1])
            cb3.set_ticklabels(["0", "1"])

    if save_pdf:
        if save_name is None:
            save_name = f"{cell_int}_vs_{gene_score}.pdf"
        save_path = os.path.join(save_dir, save_name)
        fig.savefig(save_path, bbox_inches="tight",dpi = dpi, pad_inches=0.01)
        print(f"Saved to: {save_path}")

    plt.show()


def plot_saved_spatial_mix(
    adata,
    mix_key,
    library_key="library_id",
    library_id=None,
    crop_coord=None,
    figsize=(3.5, 2.5),
    shape="square",
    size=1,
    square_scale=26,
    outline=False,
    img=True,
    title=None,
    alpha=0.9,
    show_legends=True,
    show_mix_legend=True,
    ax=None,
):
    rgb_key = f"{mix_key}_rgb"
    params_key = f"{mix_key}_color_params"
    if rgb_key not in adata.obsm:
        raise ValueError(f"{rgb_key} not found in adata.obsm. Run plot_cellint_vs_genescore(..., save_to_adata=True) first.")
    if params_key not in adata.uns:
        raise ValueError(f"{params_key} not found in adata.uns. Run plot_cellint_vs_genescore(..., save_to_adata=True) first.")

    params = adata.uns[params_key]

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    sq.pl.spatial_scatter(
        adata,
        color=None,
        library_key=library_key,
        library_id=library_id,
        crop_coord=crop_coord,
        shape=shape,
        size=size,
        outline=outline,
        img=img,
        ax=ax,
    )

    coords = np.asarray(adata.obsm["spatial"])
    point_rgb = np.asarray(adata.obsm[rgb_key])
    point_alpha = np.asarray(adata.obs[mix_key].to_numpy(dtype=float)) * alpha
    point_rgba = np.column_stack([point_rgb, point_alpha])
    marker = "s" if shape == "square" else "o"
    marker_size = size * square_scale if shape == "square" else size * 28

    ax.scatter(
        coords[:, 0],
        coords[:, 1],
        c=point_rgba,
        s=marker_size,
        marker=marker,
        linewidths=0,
    )

    if crop_coord is not None:
        x0, y0, x1, y1 = crop_coord
        ax.set_xlim(x0, x1)
        ax.set_ylim(y1, y0)

    if title is not None:
        ax.set_title(title)

    if show_legends:
        fig.subplots_adjust(right=0.76)
        cell_cmap = mpl.colors.LinearSegmentedColormap.from_list(
            f"{params['cell_int']}_white_to_color",
            ["#FFFFFF", params["cell_int_color"]],
        )
        gene_cmap = mpl.colors.LinearSegmentedColormap.from_list(
            f"{params['gene_score']}_white_to_color",
            ["#FFFFFF", params["gene_score_color"]],
        )
        mix_cmap = mpl.colors.LinearSegmentedColormap.from_list(
            f"{mix_key}_mix",
            ["#FFFFFF", params["mix_legend_color"]],
        )

        cax1 = fig.add_axes([0.88, 0.68, 0.03, 0.18])
        cax2 = fig.add_axes([0.88, 0.40, 0.03, 0.18])
        cax3 = fig.add_axes([0.88, 0.12, 0.03, 0.18]) if show_mix_legend else None

        sm1 = mpl.cm.ScalarMappable(
            norm=mpl.colors.Normalize(vmin=0, vmax=params["cell_int_vmax"]),
            cmap=cell_cmap,
        )
        sm2 = mpl.cm.ScalarMappable(
            norm=mpl.colors.Normalize(vmin=0, vmax=params["gene_score_vmax"]),
            cmap=gene_cmap,
        )

        cb1 = fig.colorbar(sm1, cax=cax1)
        cb1.ax.set_title(params["cell_int"], pad=6)
        cb1.set_ticks([0, params["cell_int_vmax"]])
        cb1.set_ticklabels(["0", f"{params['cell_int_vmax']:.2g}"])

        cb2 = fig.colorbar(sm2, cax=cax2)
        cb2.ax.set_title(params["gene_score"], pad=6)
        cb2.set_ticks([0, params["gene_score_vmax"]])
        cb2.set_ticklabels(["0", f"{params['gene_score_vmax']:.2g}"])

        if show_mix_legend:
            sm3 = mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=0, vmax=1), cmap=mix_cmap)
            cb3 = fig.colorbar(sm3, cax=cax3)
            cb3.ax.set_title(params["mix_legend_label"], pad=6)
            cb3.set_ticks([0, 1])
            cb3.set_ticklabels(["0", "1"])

    return fig, ax


def create_projected_mix_score(
    adata,
    cell_int,
    gene_score,
    output_key="ISG_mMac_V2",
    scale_quantile=0.99,
    line_start=(0.0, 0.0),
    line_end=(1.0, 1.0),
    clip_negative_gene_score=True,
    save_params=True,
):
    x = adata.obs[cell_int].to_numpy(dtype=float)
    y = adata.obs[gene_score].to_numpy(dtype=float)

    x_valid = x[np.isfinite(x)]
    x_max = np.quantile(x_valid, scale_quantile) if len(x_valid) > 0 else 1.0
    if x_max <= 0:
        x_max = 1.0

    if clip_negative_gene_score:
        y_for_scale = np.maximum(y, 0)
    else:
        y_for_scale = y.copy()

    y_valid = y_for_scale[np.isfinite(y_for_scale)]
    y_max = np.quantile(y_valid, scale_quantile) if len(y_valid) > 0 else 1.0
    if y_max <= 0:
        y_max = 1.0

    x_scaled = np.clip(np.nan_to_num(x, nan=0.0) / x_max, 0, 1)
    y_scaled = np.clip(np.nan_to_num(y_for_scale, nan=0.0) / y_max, 0, 1)

    points = np.column_stack([x_scaled, y_scaled])
    start = np.asarray(line_start, dtype=float)
    end = np.asarray(line_end, dtype=float)
    direction = end - start
    direction_norm_sq = np.dot(direction, direction)
    if direction_norm_sq == 0:
        raise ValueError("line_start and line_end must define a non-zero line.")

    projection = ((points - start) @ direction) / direction_norm_sq
    projection = np.clip(projection, 0, 1)

    adata.obs[output_key] = projection

    if save_params:
        adata.uns[f"{output_key}_params"] = {
            "cell_int": cell_int,
            "gene_score": gene_score,
            "scale_quantile": scale_quantile,
            "cell_int_vmax": float(x_max),
            "gene_score_vmax": float(y_max),
            "line_start": tuple(line_start),
            "line_end": tuple(line_end),
            "clip_negative_gene_score": clip_negative_gene_score,
        }

    return projection

def plot_same_fov_lineages(
    adata,
    cell_colors,
    library_key="library_id",
    library_id=None,
    crop_coord=None,
    figsize=(3.5, 2.5),
    shape="square",
    size=1,
    alpha=0.7,
    outline=False,
    img=True,
    dominant_type_key="dominant_parent_type",
    dominant_type_colors_key="dominant_parent_type_colors",
    spot_scale=28,
    title=None,
    show_legends=True,
    legend_height="6%",
    legend_width="22%",
    legend_y=-0.22,
    ax=None,
):
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes

    categories = pd.Categorical(adata.obs[dominant_type_key]).categories
    category_to_color = dict(zip(categories, adata.uns[dominant_type_colors_key]))

    missing = [cell_type for cell_type in cell_colors if cell_type not in category_to_color]
    if missing:
        raise ValueError(
            f"Missing colors for {missing} in adata.uns['{dominant_type_colors_key}']"
        )

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    sq.pl.spatial_scatter(
        adata,
        color=None,
        library_key=library_key,
        library_id=library_id,
        crop_coord=crop_coord,
        shape=shape,
        size=size,
        outline=outline,
        img=img,
        ax=ax,
    )

    coords = adata.obsm["spatial"]
    marker = "s" if shape == "square" else "o"
    norm = mpl.colors.Normalize(vmin=0, vmax=1)
    cmaps = {}
    value_stack = []

    for cell_type in cell_colors:
        values = adata.obs[cell_type].to_numpy(dtype=float)
        values = np.nan_to_num(values, nan=0.0)
        values = np.clip(values, 0, 1)
        value_stack.append(values)
        cmaps[cell_type] = mpl.colors.LinearSegmentedColormap.from_list(
            f"{cell_type}_white_to_color",
            ["#FFFFFF", category_to_color[cell_type]],
        )

    value_stack = np.column_stack(value_stack)
    n_spots, n_colors = value_stack.shape
    rgba_stack = np.zeros((n_spots, n_colors, 4), dtype=float)

    for i, cell_type in enumerate(cell_colors):
        rgba_stack[:, i, :] = cmaps[cell_type](norm(value_stack[:, i]))
        rgba_stack[:, i, 3] = value_stack[:, i] * alpha

    weights = value_stack**2
    weight_sum = weights.sum(axis=1, keepdims=True)
    safe_weight_sum = np.where(weight_sum == 0, 1, weight_sum)
    blended_rgb = (rgba_stack[:, :, :3] * weights[:, :, None]).sum(axis=1) / safe_weight_sum
    blended_alpha = rgba_stack[:, :, 3].max(axis=1)

    no_signal = weight_sum.ravel() == 0
    blended_rgb[no_signal] = 1.0
    blended_alpha[no_signal] = 0.0

    blended_rgba = np.column_stack([blended_rgb, blended_alpha])
    valid = blended_alpha > 0

    ax.scatter(
        coords[valid, 0],
        coords[valid, 1],
        c=blended_rgba[valid],
        s=size * spot_scale,
        marker=marker,
        linewidths=0,
    )

    if crop_coord is not None:
        x0, y0, x1, y1 = crop_coord
        ax.set_xlim(x0, x1)
        ax.set_ylim(y1, y0)

    if title is not None:
        ax.set_title(title)

    if show_legends:
        n_legends = len(cell_colors)
        for i, cell_type in enumerate(cell_colors):
            x0 = i / max(n_legends, 1)
            cax = inset_axes(
                ax,
                width=legend_width,
                height=legend_height,
                loc="lower left",
                bbox_to_anchor=(x0, legend_y, 1, 1),
                bbox_transform=ax.transAxes,
                borderpad=0,
            )
            sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmaps[cell_type])
            cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
            cb.ax.set_title(cell_type, pad=2)
            cb.set_ticks([0, 1])
            cb.set_ticklabels(["0", "1"])

    return fig, ax

def signature_score_per_cell_min_detected(
    data,
    gene_set,
    score_title,
    min_detected_genes=1,
    expr_threshold=0,
    fill_value=np.nan,   # use 0 if you want background plotted as zero
):
    orig_len = len(gene_set)
    gene_set = [gene for gene in gene_set if gene in data.var_names]
    print(f"{len(gene_set)}/{orig_len} of the gene set genes are measured")

    dat = data[:, gene_set].X
    if hasattr(dat, "toarray"):
        dat = dat.toarray()

    # bins/cells with at least N genes detected
    detected_mask = (dat > expr_threshold).sum(axis=1) >= min_detected_genes

    mean = dat.mean(axis=0)
    std = dat.std(axis=0)
    std[std == 0] = 0.0001

    with np.errstate(divide="ignore", invalid="ignore"):
        z = (dat - mean) / std

    z[z < -5] = -5
    z[z > 5] = 5

    scores = z.mean(axis=1)

    # only keep score for bins with at least one detected ISG
    filtered_scores = np.full(data.n_obs, fill_value, dtype=float)
    filtered_scores[detected_mask] = scores[detected_mask]

    data.obs[score_title] = filtered_scores
    data.obs[f"{score_title}_n_detected"] = (dat > expr_threshold).sum(axis=1)
    data.obs[f"{score_title}_detected"] = detected_mask.astype(int)