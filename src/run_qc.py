# src/run_qc.py
# HydroQC Mini: quick QC report for bathymetry/DTM GeoTIFFs
#
# Outputs:
# - outputs/summary.json
# - outputs/report.md
# - outputs/figures/depth.png
# - outputs/figures/holes.png
# - outputs/figures/outliers.png
# - outputs/figures/roughness.png
# - outputs/figures/slope_deg.png
#
# Usage (positional input):
#   python .\src\run_qc.py .\data_raw\D5_2024.tif --out .\outputs
#
# Usage (optional --input also supported):
#   python .\src\run_qc.py --input .\data_raw\D5_2024.tif --out .\outputs

from __future__ import annotations

import argparse
import json
import math
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.warp import calculate_default_transform, reproject
from rasterio.crs import CRS

import matplotlib.pyplot as plt

# Optional, but strongly recommended for fast local-window ops and connected components
try:
    from scipy.ndimage import uniform_filter, label
    _HAS_SCIPY = True
except Exception:
    _HAS_SCIPY = False


@dataclass
class QCParams:
    outlier_mad_z: float = 6.0
    roughness_window: int = 9
    roughness_min_valid_frac: float = 0.6
    roughness_clip_percentile: float = 99.0
    max_plot_dim: int = 2000


def _ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def _pick_utm_crs(lon: float, lat: float) -> CRS:
    # UTM zone from longitude; EPSG 326xx for north, 327xx for south
    zone = int((lon + 180.0) / 6.0) + 1
    if lat >= 0:
        epsg = 32600 + zone
    else:
        epsg = 32700 + zone
    return CRS.from_epsg(epsg)


def _scaled_transform(transform: rasterio.Affine, src_w: int, src_h: int, dst_w: int, dst_h: int) -> rasterio.Affine:
    # When using rasterio.read(... out_shape=...), pixel size changes by scale factors.
    # This updates the affine transform accordingly.
    scale_x = src_w / float(dst_w)
    scale_y = src_h / float(dst_h)
    return transform * rasterio.Affine.scale(scale_x, scale_y)


def _robust_outlier_mask(x: np.ndarray, z_thresh: float) -> np.ndarray:
    """
    Robust outlier mask using MAD-based z-score:
      z = 0.6745 * (x - median) / MAD
    Returns mask (True = outlier). NaNs are treated as not-outlier here (mask separately).
    """
    vals = x[np.isfinite(x)]
    if vals.size < 50:
        return np.zeros_like(x, dtype=bool)

    med = np.median(vals)
    mad = np.median(np.abs(vals - med))
    if mad <= 1e-12:
        return np.zeros_like(x, dtype=bool)

    z = 0.6745 * (x - med) / mad
    m = np.isfinite(z) & (np.abs(z) > z_thresh)
    return m


def _online_stats_update(count: int, mean: float, M2: float, new_vals: np.ndarray) -> Tuple[int, float, float]:
    # Welford update for mean/std over a 1D array
    for v in new_vals:
        count += 1
        delta = v - mean
        mean += delta / count
        delta2 = v - mean
        M2 += delta * delta2
    return count, mean, M2


def _compute_global_stats(ds: rasterio.DatasetReader) -> Dict[str, Any]:
    """
    Compute stats (valid %, min/max/mean/std) in a memory-safe way using block windows.
    Valid is finite values (not NaN, not +/-inf).
    """
    band = 1
    total = ds.width * ds.height

    valid_count = 0
    vmin = np.inf
    vmax = -np.inf

    count = 0
    mean = 0.0
    M2 = 0.0

    # Use dataset block windows when tiled, otherwise iterate with a reasonable window grid
    try:
        windows = list(ds.block_windows(band))
    except Exception:
        windows = []
        step = 1024
        for row_off in range(0, ds.height, step):
            for col_off in range(0, ds.width, step):
                h = min(step, ds.height - row_off)
                w = min(step, ds.width - col_off)
                windows.append(((row_off // step, col_off // step), rasterio.windows.Window(col_off, row_off, w, h)))

    for _, w in windows:
        arr = ds.read(band, window=w, masked=False).astype(np.float64, copy=False)
        finite = np.isfinite(arr)
        if not np.any(finite):
            continue

        vals = arr[finite]
        valid_count += int(vals.size)

        vmin = min(vmin, float(np.min(vals)))
        vmax = max(vmax, float(np.max(vals)))

        count, mean, M2 = _online_stats_update(count, mean, M2, vals.ravel())

    if valid_count == 0:
        return {
            "valid_percent": 0.0,
            "min": None,
            "max": None,
            "mean": None,
            "std": None,
        }

    var = (M2 / (count - 1)) if count > 1 else 0.0
    std = float(math.sqrt(max(var, 0.0)))

    return {
        "valid_percent": 100.0 * (valid_count / float(total)),
        "min": float(vmin),
        "max": float(vmax),
        "mean": float(mean),
        "std": std,
    }


def _read_downsample(ds: rasterio.DatasetReader, max_dim: int) -> Tuple[np.ndarray, rasterio.Affine]:
    """
    Read a downsampled array for plotting and spatial QC layers.
    Keeps aspect ratio and uses average resampling.
    Returns (array float32 with NaNs, transform for the downsample grid).
    """
    band = 1
    src_w, src_h = ds.width, ds.height
    scale = min(1.0, max_dim / float(max(src_w, src_h)))
    dst_w = max(2, int(round(src_w * scale)))
    dst_h = max(2, int(round(src_h * scale)))

    arr = ds.read(
        band,
        out_shape=(dst_h, dst_w),
        resampling=Resampling.average,
        masked=False,
    ).astype(np.float32)

    # Normalize nodata to NaN
    nod = ds.nodata
    if nod is not None and np.isfinite(nod):
        arr = np.where(arr == nod, np.nan, arr)

    # Some EMODnet files store NaN nodata already; keep it.
    # Also guard inf.
    arr = np.where(np.isfinite(arr), arr, np.nan)

    ds_transform = ds.transform
    dst_transform = _scaled_transform(ds_transform, src_w, src_h, dst_w, dst_h)
    return arr, dst_transform


def _connected_missing_regions(missing_mask: np.ndarray) -> Tuple[int, int]:
    """
    Returns (num_regions, largest_region_pixels).
    Runs on downsampled boolean mask, True = missing.
    """
    if not np.any(missing_mask):
        return 0, 0

    if _HAS_SCIPY:
        structure = np.array([[0, 1, 0],
                              [1, 1, 1],
                              [0, 1, 0]], dtype=np.int8)  # 4-connected
        lab, n = label(missing_mask.astype(np.uint8), structure=structure)
        if n == 0:
            return 0, 0
        counts = np.bincount(lab.ravel())
        # counts[0] is background
        largest = int(np.max(counts[1:])) if counts.size > 1 else 0
        return int(n), largest

    # Fallback: very simple (and slower) BFS on downsample.
    # Works if your downsample is not huge (<= ~2000x2000).
    h, w = missing_mask.shape
    visited = np.zeros_like(missing_mask, dtype=bool)
    regions = 0
    largest = 0

    for y in range(h):
        for x in range(w):
            if not missing_mask[y, x] or visited[y, x]:
                continue
            regions += 1
            stack = [(y, x)]
            visited[y, x] = True
            size = 0
            while stack:
                cy, cx = stack.pop()
                size += 1
                for ny, nx in ((cy - 1, cx), (cy + 1, cx), (cy, cx - 1), (cy, cx + 1)):
                    if 0 <= ny < h and 0 <= nx < w and missing_mask[ny, nx] and not visited[ny, nx]:
                        visited[ny, nx] = True
                        stack.append((ny, nx))
            largest = max(largest, size)

    return regions, largest


def _reproject_to_utm(
    src: np.ndarray,
    src_crs: CRS,
    src_transform: rasterio.Affine,
    bounds_lonlat: Tuple[float, float, float, float],
) -> Tuple[np.ndarray, rasterio.Affine, CRS]:
    """
    Reproject a float array (with NaNs) to an auto-selected UTM CRS.
    """
    left, bottom, right, top = bounds_lonlat
    lon_c = (left + right) / 2.0
    lat_c = (bottom + top) / 2.0
    dst_crs = _pick_utm_crs(lon_c, lat_c)

    dst_transform, dst_w, dst_h = calculate_default_transform(
        src_crs, dst_crs, src.shape[1], src.shape[0], left, bottom, right, top
    )

    dst = np.full((dst_h, dst_w), np.nan, dtype=np.float32)

    # rasterio.reproject handles NaN nodata if we pass src_nodata/dst_nodata as NaN
    reproject(
        source=src.astype(np.float32),
        destination=dst,
        src_transform=src_transform,
        src_crs=src_crs,
        dst_transform=dst_transform,
        dst_crs=dst_crs,
        src_nodata=np.nan,
        dst_nodata=np.nan,
        resampling=Resampling.bilinear,
    )
    dst = np.where(np.isfinite(dst), dst, np.nan)
    return dst, dst_transform, dst_crs


def _slope_degrees(z: np.ndarray, transform: rasterio.Affine) -> np.ndarray:
    """
    Compute slope (degrees) using central differences on a projected grid (meters).
    Masks pixels where required neighbors are NaN.
    """
    dx = abs(transform.a)
    dy = abs(transform.e)

    slope = np.full_like(z, np.nan, dtype=np.float32)
    if z.shape[0] < 3 or z.shape[1] < 3:
        return slope

    # central differences
    zc = z[1:-1, 1:-1]
    zx1 = z[1:-1, 2:]
    zx0 = z[1:-1, :-2]
    zy1 = z[2:, 1:-1]
    zy0 = z[:-2, 1:-1]

    ok = np.isfinite(zc) & np.isfinite(zx1) & np.isfinite(zx0) & np.isfinite(zy1) & np.isfinite(zy0)
    dzdx = (zx1 - zx0) / (2.0 * dx)
    dzdy = (zy1 - zy0) / (2.0 * dy)

    s = np.arctan(np.sqrt(dzdx * dzdx + dzdy * dzdy))
    sdeg = np.degrees(s).astype(np.float32)

    inner = slope[1:-1, 1:-1]
    inner[ok] = sdeg[ok]
    slope[1:-1, 1:-1] = inner
    return slope


def _roughness_abs_dev(z: np.ndarray, window: int, min_valid_frac: float) -> np.ndarray:
    """
    Roughness = |z - local_mean(z)| using a square window.
    Uses sum/count via uniform_filter if scipy is available.
    """
    if window < 3:
        window = 3
    if window % 2 == 0:
        window += 1

    rough = np.full_like(z, np.nan, dtype=np.float32)
    if not _HAS_SCIPY:
        # Fallback: compute roughness only on a small grid would be needed.
        # Raise a clear error instead of silently producing nonsense.
        raise RuntimeError("scipy is required for roughness (pip install scipy).")

    valid = np.isfinite(z).astype(np.float32)
    z0 = np.where(np.isfinite(z), z, 0.0).astype(np.float32)

    # local sum and local count
    local_sum = uniform_filter(z0, size=window, mode="nearest") * (window * window)
    local_cnt = uniform_filter(valid, size=window, mode="nearest") * (window * window)

    # avoid div-by-zero
    mean = np.where(local_cnt > 0, local_sum / local_cnt, np.nan)
    dev = np.abs(z - mean)

    min_cnt = min_valid_frac * (window * window)
    dev = np.where(local_cnt >= min_cnt, dev, np.nan).astype(np.float32)
    return dev


def _robust_vmin_vmax(arr: np.ndarray, lo: float = 2.0, hi: float = 98.0) -> Tuple[float, float]:
    vals = arr[np.isfinite(arr)]
    if vals.size < 50:
        return (float(np.nanmin(arr)), float(np.nanmax(arr)))
    vmin = float(np.percentile(vals, lo))
    vmax = float(np.percentile(vals, hi))
    if not np.isfinite(vmin) or not np.isfinite(vmax) or vmin == vmax:
        vmin = float(np.nanmin(arr))
        vmax = float(np.nanmax(arr))
    return vmin, vmax


def _save_fig(path: Path) -> None:
    plt.tight_layout()
    plt.savefig(path, dpi=160, bbox_inches="tight")
    plt.close()


def main() -> None:
    p = argparse.ArgumentParser(description="HydroQC Mini: quick QC report for bathymetry/DTM rasters")
    # Support both positional input and --input
    p.add_argument("input_pos", nargs="?", help="Path to input GeoTIFF (positional).")
    p.add_argument("--input", dest="input_opt", help="Path to input GeoTIFF (optional flag).")

    p.add_argument("--out", default="outputs", help="Output folder.")
    p.add_argument("--outlier-mad-z", type=float, default=6.0, help="MAD z-threshold for outliers (default: 6).")
    p.add_argument("--roughness-window", type=int, default=9, help="Window size for roughness (odd int, default: 9).")
    p.add_argument("--roughness-min-valid-frac", type=float, default=0.6, help="Min valid fraction in window (default: 0.6).")
    p.add_argument("--roughness-clip-percentile", type=float, default=99.0, help="Clip roughness for plotting (default: 99).")
    p.add_argument("--max-plot-dim", type=int, default=2000, help="Max plot dimension for downsample (default: 2000).")

    args = p.parse_args()

    in_path = args.input_opt or args.input_pos
    if not in_path:
        raise SystemExit("No input provided. Use positional input or --input.")

    params = QCParams(
        outlier_mad_z=float(args.outlier_mad_z),
        roughness_window=int(args.roughness_window),
        roughness_min_valid_frac=float(args.roughness_min_valid_frac),
        roughness_clip_percentile=float(args.roughness_clip_percentile),
        max_plot_dim=int(args.max_plot_dim),
    )

    in_path = Path(in_path)
    out_dir = Path(args.out)
    fig_dir = out_dir / "figures"
    _ensure_dir(out_dir)
    _ensure_dir(fig_dir)

    with rasterio.open(in_path) as ds:
        crs = ds.crs
        if crs is None:
            raise RuntimeError("Input raster has no CRS.")

        bounds = ds.bounds  # in raster CRS
        # Expect lon/lat if EPSG:4326, but keep generic text
        bounds_tuple = (float(bounds.left), float(bounds.bottom), float(bounds.right), float(bounds.top))

        # Global stats (streaming)
        gstats = _compute_global_stats(ds)

        # Downsampled array for plots and spatial QC
        depth_ds, depth_transform = _read_downsample(ds, params.max_plot_dim)

        missing = ~np.isfinite(depth_ds)
        outliers = _robust_outlier_mask(depth_ds, params.outlier_mad_z)

        # Missing connected components on downsample
        n_regions, largest_region = _connected_missing_regions(missing)

        # Depth plot (mask missing + outliers)
        depth_plot = depth_ds.copy()
        depth_plot[missing] = np.nan
        depth_plot[outliers] = np.nan

        dvmin, dvmax = _robust_vmin_vmax(depth_plot, 2.0, 98.0)

        plt.figure(figsize=(10, 7))
        plt.title("Depth (masked NaN + robust outliers)")
        im = plt.imshow(depth_plot, vmin=dvmin, vmax=dvmax)
        plt.colorbar(im, fraction=0.03, pad=0.04)
        _save_fig(fig_dir / "depth.png")

        # NoData mask figure
        plt.figure(figsize=(10, 7))
        plt.title("NoData mask (1 = missing)")
        im = plt.imshow(missing.astype(np.float32), vmin=0.0, vmax=1.0, cmap="gray")
        plt.colorbar(im, fraction=0.03, pad=0.04)
        _save_fig(fig_dir / "holes.png")

        # Outlier mask figure
        plt.figure(figsize=(10, 7))
        plt.title("Outlier mask (1 = outlier)")
        im = plt.imshow(outliers.astype(np.float32), vmin=0.0, vmax=1.0, cmap="gray")
        plt.colorbar(im, fraction=0.03, pad=0.04)
        _save_fig(fig_dir / "outliers.png")

        # Slope + roughness computed after reprojection to UTM (meters)
        utm_arr, utm_transform, utm_crs = _reproject_to_utm(
            depth_plot,  # already masks missing/outliers to NaN
            src_crs=crs,
            src_transform=depth_transform,
            bounds_lonlat=bounds_tuple,
        )

        slope = _slope_degrees(utm_arr, utm_transform)

        # Roughness
        rough = _roughness_abs_dev(utm_arr, params.roughness_window, params.roughness_min_valid_frac)

        # Clip roughness for nicer visualization (keep NaNs)
        rvals = rough[np.isfinite(rough)]
        if rvals.size > 50:
            rclip = float(np.percentile(rvals, params.roughness_clip_percentile))
            rough_plot = np.clip(rough, 0.0, rclip)
        else:
            rough_plot = rough

        # Plot slope
        svmin, svmax = _robust_vmin_vmax(slope, 2.0, 98.0)
        plt.figure(figsize=(10, 7))
        plt.title("Slope (degrees) after UTM reprojection")
        im = plt.imshow(slope, vmin=max(0.0, svmin), vmax=svmax)
        plt.colorbar(im, fraction=0.03, pad=0.04)
        _save_fig(fig_dir / "slope_deg.png")

        # Plot roughness
        rvmin, rvmax = _robust_vmin_vmax(rough_plot, 2.0, 98.0)
        plt.figure(figsize=(10, 7))
        plt.title("Roughness (abs deviation from local mean)")
        im = plt.imshow(rough_plot, vmin=max(0.0, rvmin), vmax=rvmax)
        plt.colorbar(im, fraction=0.03, pad=0.04)
        _save_fig(fig_dir / "roughness.png")

        # Summary JSON
        summary: Dict[str, Any] = {
            "width": int(ds.width),
            "height": int(ds.height),
            "crs": str(crs),
            "bounds": [bounds_tuple[0], bounds_tuple[1], bounds_tuple[2], bounds_tuple[3]],
            "dtype": str(ds.dtypes[0]),
            "nodata": "nan" if ds.nodata is None or (isinstance(ds.nodata, float) and np.isnan(ds.nodata)) else ds.nodata,
            "scale": 1.0,
            "offset": 0.0,
            "valid_percent": float(gstats["valid_percent"]),
            "min": gstats["min"],
            "max": gstats["max"],
            "mean": gstats["mean"],
            "std": gstats["std"],
            "missing_regions_downsample": int(n_regions),
            "largest_missing_region_pixels_downsample": int(largest_region),
            "downsample_shape": [int(depth_ds.shape[0]), int(depth_ds.shape[1])],
            "params": {
                "outlier_mad_z": params.outlier_mad_z,
                "roughness_window": params.roughness_window,
                "roughness_min_valid_frac": params.roughness_min_valid_frac,
                "roughness_clip_percentile": params.roughness_clip_percentile,
                "max_plot_dim": params.max_plot_dim,
            },
            "utm_crs": str(utm_crs),
        }

    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    # Report.md
    report_md = f"""# HydroQC Mini Report

## Input
- File: `{in_path.name}`
- CRS: `{summary["crs"]}`
- Raster size: {summary["width"]} x {summary["height"]}
- Bounds (in CRS units): ({summary["bounds"][0]}, {summary["bounds"][1]}, {summary["bounds"][2]}, {summary["bounds"][3]})
- Data type: {summary["dtype"]}
- NoData: {summary["nodata"]}

## Parameters
- Outlier MAD z-threshold: {summary["params"]["outlier_mad_z"]}
- Roughness window: {summary["params"]["roughness_window"]}
- Roughness min valid fraction: {summary["params"]["roughness_min_valid_frac"]}
- Roughness clip percentile (plot): {summary["params"]["roughness_clip_percentile"]}
- Max plot dimension (downsample): {summary["params"]["max_plot_dim"]}
- UTM CRS used for slope/roughness: `{summary["utm_crs"]}`

## Quick QC metrics (global, streamed)
- Valid pixels: {summary["valid_percent"]:.2f}%
- Min depth: {summary["min"]}
- Max depth: {summary["max"]}
- Mean depth: {summary["mean"]}
- Std depth: {summary["std"]}

## Missing data (downsample view)
- Connected missing regions (count): {summary["missing_regions_downsample"]}
- Largest missing region (pixels, downsample): {summary["largest_missing_region_pixels_downsample"]}
- Downsample shape: {summary["downsample_shape"][1]} x {summary["downsample_shape"][0]}

## Figures
![Depth](figures/depth.png)

![NoData mask](figures/holes.png)

![Outliers](figures/outliers.png)

![Roughness](figures/roughness.png)

![Slope](figures/slope_deg.png)
"""
    (out_dir / "report.md").write_text(report_md, encoding="utf-8")


if __name__ == "__main__":
    main()
