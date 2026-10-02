"""Plot suitability and bus-level renewable profile diagnostics."""

import argparse
from pathlib import Path

import geopandas as gpd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import pandas as pd
import rasterio
import xarray as xr
import yaml
from matplotlib import colors
from rasterio.enums import Resampling
from rasterio.plot import plotting_extent
from rasterio.windows import from_bounds


ROOT = Path(__file__).resolve().parents[1]
PALETTES = {"solar": "YlGn", "onwind": "PuBuGn"}


def project_path(path):
    path = Path(path)
    return path if path.is_absolute() else ROOT / path


def bus_label(name):
    if name.endswith(("_q1", "_q2")):
        site, quality_class = name.rsplit("_", 1)
        return f"{bus_label(site)} {quality_class.upper()}"
    if "_site_" in name:
        island, site = name.split("_site_", 1)
        return f"M{island.removeprefix('microgrid_')}-S{int(site) + 1}"
    return name


def profile_metrics(profile_path, regions):
    with xr.open_dataset(profile_path) as dataset:
        if not {"bus", "time"}.issubset(dataset.profile.dims):
            raise ValueError(f"Missing bus/time profile dimensions in {profile_path}")
        buses = dataset.bus.to_index().astype(str)
        region_buses = pd.Index(regions["name"].astype(str))
        sites = dataset.site.to_index().astype(str) if "site" in dataset.coords else buses
        if buses.empty or not buses.is_unique or not region_buses.is_unique or not set(sites).issubset(region_buses):
            raise ValueError(f"Profile buses must be a nonempty subset of region IDs: {profile_path}")

        times = pd.DatetimeIndex(dataset.time.values)
        if len(times) < 2:
            raise ValueError("At least two profile timestamps are needed")
        intervals = np.diff(times.asi8) / 3.6e12
        if not np.allclose(intervals, intervals[0]) or intervals[0] <= 0:
            raise ValueError("Profile timestamps must have a regular positive interval")

        values = dataset.profile.transpose("time", "bus").to_numpy()
        capacity = dataset.p_nom_max.to_numpy()
        if not np.isfinite(values).all() or not np.isfinite(capacity).all():
            raise ValueError(f"Non-finite profile or capacity values in {profile_path}")
        metrics = pd.DataFrame(
            {
                "p_nom_max": capacity,
                "full_load_hours": values.sum(axis=0) * intervals[0],
                "capacity_factor": values.mean(axis=0) * 100,
                "weight": dataset.weight.to_numpy() if "weight" in dataset else capacity,
                "site": sites,
            },
            index=buses,
        )
        metrics.index.name = "bus"
        return metrics, times[0], times[-1]


def draw_bus_labels(ax, regions):
    if {"x", "y"}.issubset(regions.columns):
        positions = gpd.GeoSeries(gpd.points_from_xy(regions.x, regions.y), crs=4326).to_crs(regions.crs)
    else:
        positions = regions.geometry.representative_point()
    for (_, region), point in zip(regions.iterrows(), positions):
        ax.text(
            point.x, point.y, bus_label(str(region["name"])), ha="center", va="center",
            fontsize=7, fontweight="bold", color="#132d31",
            bbox={"boxstyle": "round,pad=0.15", "facecolor": "white", "edgecolor": "none", "alpha": 0.85},
        )


def draw_suitability(ax, raster_path, regions, threshold, palette):
    with rasterio.open(raster_path) as raster:
        if raster.crs is None or raster.count != 1:
            raise ValueError(f"Suitability raster needs a CRS and one band: {raster_path}")
        projected = regions.to_crs(raster.crs)
        west, south, east, north = projected.total_bounds
        padding = 0.06 * max(east - west, north - south)
        window = from_bounds(west - padding, south - padding, east + padding, north + padding,
                             raster.transform)
        height = min(max(round(window.height), 1), 1400)
        width = min(max(round(window.width), 1), 1400)
        scores = raster.read(1, window=window, out_shape=(height, width), boundless=True,
                             resampling=Resampling.nearest, masked=True)
        transform = raster.window_transform(window) * raster.transform.scale(window.width / width, window.height / height)
        bounds = plotting_extent(scores, transform)
        valid = ~np.ma.getmaskarray(scores) & np.isfinite(scores.data)
        eligible = valid & (scores.data >= threshold)
        if not eligible.any():
            raise ValueError(f"No raster pixels meet the suitability cutoff: {raster_path}")

        ax.set_facecolor("#e5e9e7")
        ax.imshow(np.where(eligible, scores.data, np.nan), extent=bounds, origin="upper",
                  cmap=palette, vmin=threshold, vmax=max(float(scores.data[eligible].max()), threshold + 0.1),
                  interpolation="nearest")
        projected.boundary.plot(ax=ax, color="#1c4345", linewidth=1.1)
        draw_bus_labels(ax, projected)
        ax.set_xlim(west - padding, east + padding)
        ax.set_ylim(south - padding, north + padding)
        ax.legend(handles=[patches.Patch(facecolor="#e5e9e7", label="Below cutoff / no score"),
                           patches.Patch(facecolor=plt.get_cmap(palette)(0.8), label="Eligible")],
                  loc="lower left", frameon=True, fontsize=9)
        ax.set_title(f"Suitable pixels (score >= {threshold:g})", loc="left", fontsize=13, fontweight="bold")
        ax.set_axis_off()
        return projected


def draw_metric_map(ax, regions, metrics, key, title, unit, palette):
    projected = regions.copy()
    if key == "p_nom_max":
        by_site = metrics.groupby("site")[key].sum()
    else:
        weighted = (metrics[key] * metrics["weight"]).groupby(metrics["site"]).sum()
        by_site = weighted / metrics.groupby("site")["weight"].sum()
    projected[key] = projected["name"].astype(str).map(by_site)
    maximum = max(float(metrics[key].max()), 1e-9)
    norm = colors.Normalize(vmin=0, vmax=maximum)
    projected.plot(ax=ax, column=key, cmap=palette, norm=norm, edgecolor="#25494a",
                   linewidth=0.9, missing_kwds={"color": "#e5e9e7"})
    draw_bus_labels(ax, projected)
    scalar = plt.cm.ScalarMappable(norm=norm, cmap=palette)
    ax.figure.colorbar(scalar, ax=ax, fraction=0.035, pad=0.02, label=unit)
    ax.set_title(title, loc="left", fontsize=12, fontweight="bold")
    ax.set_axis_off()


def draw_sorted_bars(ax, metrics, key, title, unit, palette):
    ordered = metrics[key].sort_values(ascending=False)
    maximum = max(float(ordered.max()), 1e-9)
    shades = plt.get_cmap(palette)(colors.Normalize(vmin=0, vmax=maximum)(ordered.to_numpy()))
    ax.barh([bus_label(bus) for bus in ordered.index], ordered.to_numpy(), color=shades, height=0.62)
    ax.invert_yaxis()
    ax.set_xlim(0, maximum * 1.19)
    for row, value in enumerate(ordered):
        ax.text(value + maximum * 0.025, row, f"{value:,.0f}" if key != "capacity_factor" else f"{value:.1f}",
                va="center", color="#25494a", fontsize=10)
    ax.set_xlabel(unit)
    ax.set_title(title, loc="left", fontsize=12, fontweight="bold")
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.grid(axis="x", color="#d7dedb", linewidth=0.7)
    ax.set_axisbelow(True)


def plot_technology(technology, config, regions_path, profiles_dir, output_dir):
    settings = config["renewable"][technology]["suitability"]
    threshold = float(settings["threshold"])
    if not np.isfinite(threshold):
        raise ValueError(f"Suitability threshold must be finite: {technology}")
    raster_path = project_path(settings["path"])
    regions = gpd.read_file(regions_path)
    if regions.empty or regions.crs is None or "name" not in regions:
        raise ValueError("Bus regions need a CRS and unique names")
    profile_path = profiles_dir / f"profile_{technology}.nc"
    metrics, first, last = profile_metrics(profile_path, regions)
    palette = PALETTES[technology]

    fig = plt.figure(figsize=(15, 19), facecolor="#f7f9f7", layout="constrained")
    grid = fig.add_gridspec(4, 2, height_ratios=[1.65, 1, 1, 1], hspace=0.12, wspace=0.17)
    raster_ax = fig.add_subplot(grid[0, :])
    projected = draw_suitability(raster_ax, raster_path, regions, threshold, palette)
    diagnostics = [
        ("p_nom_max", "Maximum installable capacity", "MW"),
        ("full_load_hours", "Full-load hours", "h over profile period"),
        ("capacity_factor", "Mean capacity factor", "%"),
    ]
    for row, (key, title, unit) in enumerate(diagnostics, start=1):
        draw_metric_map(fig.add_subplot(grid[row, 0]), projected, metrics, key, title, unit, palette)
        draw_sorted_bars(fig.add_subplot(grid[row, 1]), metrics, key, title, unit, palette)

    fig.suptitle(f"{technology.upper()}  |  Suitability and renewable potential", x=0.05, ha="left",
                 fontsize=20, fontweight="bold", color="#173e40")
    fig.text(0.05, 0.005, f"Bus profiles  {first:%d %b %Y} to {last:%d %b %Y}  |  Gray pixels are ineligible",
             fontsize=10, color="#556e6a")
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"{technology}_suitability_profiles.png"
    fig.savefig(output_path, dpi=180, facecolor=fig.get_facecolor())
    plt.close(fig)
    print(output_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=ROOT / "configs/config.hygalapagos.yaml")
    parser.add_argument("--regions", type=Path, help="Override the scenario's renewable bus regions")
    parser.add_argument("--profiles-dir", type=Path, default=ROOT / "resources/renewable_profiles")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "resources/renewable_profiles/plots")
    parser.add_argument("--technology", choices=["solar", "onwind"], nargs="+", default=["solar", "onwind"])
    args = parser.parse_args()
    with args.config.open() as config_file:
        scenario = yaml.safe_load(config_file)
    regions_path = args.regions
    if regions_path is None:
        site_mode = scenario.get("mode") == "green_field" and scenario.get("generation_sites", {}).get("enabled", False)
        regions_path = ROOT / (
            "resources/shapes/candidate_sites.geojson" if site_mode
            else "resources/bus_regions/regions_onshore.geojson"
        )
    for technology in args.technology:
        plot_technology(technology, scenario, regions_path, args.profiles_dir, args.output_dir)