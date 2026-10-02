import logging

import geopandas as gpd
import numpy as np
import rasterio
from rasterio.features import geometry_mask
from rasterio.warp import Resampling, reproject
from rasterio.windows import Window, from_bounds
from scipy.ndimage import distance_transform_edt, label
from shapely.geometry import MultiPoint, Point
from shapely import voronoi_polygons

from _helpers_dist import configure_logging, sets_path_to_root

logger = logging.getLogger(__name__)


def build_candidate_sites(regions, clusters, suitability, settings):
    if not 1 <= settings["min_sites"] <= settings["max_sites"]:
        raise ValueError("Candidate site minimum must be positive and no greater than maximum")
    if any(settings[key] <= 0 for key in (
        "min_spacing_km", "min_eligible_area_km2", "max_connection_km"
    )):
        raise ValueError("Candidate site area, spacing and connection range must be positive")
    regions = regions.to_crs("EPSG:4326")
    candidates = []
    with rasterio.open(suitability["onwind"]["path"]) as grid:
        if grid.crs is None or grid.res[0] <= 0 or grid.res[1] <= 0:
            raise ValueError("Wind suitability raster needs a metric CRS and valid resolution")
        metric_regions = regions.to_crs(grid.crs)
        metric_clusters = clusters.to_crs(grid.crs)
        pixel_area = abs(grid.transform.a * grid.transform.e)
        minimum_pixels = max(1, int(np.ceil(settings["min_eligible_area_km2"] * 1e6 / pixel_area)))

        for region in metric_regions.itertuples():
            island = region.name_microgrid.removesuffix("_gen_bus")
            load_points = metric_clusters.loc[
                metric_clusters.name_microgrid == island, "geometry"
            ]
            if load_points.empty:
                raise ValueError(f"No load clusters found for {island}")
            window = from_bounds(*region.geometry.bounds, transform=grid.transform)
            window = window.round_offsets().round_lengths().intersection(
                Window(0, 0, grid.width, grid.height)
            )
            transform = grid.window_transform(window)
            shape = (int(window.height), int(window.width))
            inside = geometry_mask(
                [region.geometry], out_shape=shape, transform=transform, invert=True
            )
            masks = {}
            scores = {}
            for technology, spec in suitability.items():
                with rasterio.open(spec["path"]) as source:
                    if source.count != 1 or source.crs is None:
                        raise ValueError(f"Invalid suitability raster for {technology}")
                    values = np.full(shape, np.nan, dtype=np.float32)
                    reproject(
                        source=rasterio.band(source, 1),
                        destination=values,
                        src_nodata=source.nodata,
                        dst_nodata=np.nan,
                        dst_transform=transform,
                        dst_crs=grid.crs,
                        resampling=Resampling.nearest,
                    )
                    masks[technology] = inside & np.isfinite(values) & (
                        values >= spec["threshold"]
                    )
                    scores[technology] = values

            eligible = np.logical_or.reduce(list(masks.values()))
            components, count = label(eligible, structure=np.ones((3, 3)))
            if not count:
                raise ValueError(f"No eligible generation area in {island}")

            rows, cols = np.indices(shape)
            east = transform.c + (cols + 0.5) * transform.a
            north = transform.f + (rows + 0.5) * transform.e
            distance_to_load = np.minimum.reduce(
                [np.hypot(east - point.x, north - point.y) for point in load_points]
            )
            permitted = distance_to_load <= settings["max_connection_km"] * 1000
            interior = distance_transform_edt(eligible) * abs(grid.res[0])
            score = np.zeros(shape, dtype=float)
            for technology, mask in masks.items():
                if mask.any():
                    threshold = suitability[technology]["threshold"]
                    upper = np.percentile(scores[technology][mask], 95)
                    quality = np.clip(
                        (scores[technology] - threshold) / max(upper - threshold, 1e-6),
                        0, 1,
                    )
                    score += np.where(mask, 0.5 + 0.5 * quality, 0) / mask.sum()
            score *= (0.5 + np.minimum(interior / 1000, 1)) / (
                1 + distance_to_load / 20000
            )
            score[~permitted] = -np.inf

            sizes = np.bincount(components.ravel())
            ranked = []
            for component in range(1, count + 1):
                if sizes[component] < minimum_pixels:
                    continue
                component_score = np.where(components == component, score, -np.inf)
                best = np.unravel_index(np.argmax(component_score), shape)
                if np.isfinite(component_score[best]):
                    ranked.append((float(component_score[best]), best))
            ranked.sort(key=lambda item: (-item[0], item[1]))
            chosen = []
            spacing = settings["min_spacing_km"] * 1000
            for _, location in ranked:
                if all(
                    np.hypot(east[location] - east[other], north[location] - north[other]) >= spacing
                    for other in chosen
                ):
                    chosen.append(location)
                if len(chosen) == settings["max_sites"]:
                    break

            while chosen and len(chosen) < settings["min_sites"]:
                separation = np.minimum.reduce(
                    [np.hypot(east - east[point], north - north[point]) for point in chosen]
                )
                extra_score = np.where(
                    eligible & permitted & (separation >= spacing),
                    score * separation / spacing,
                    -np.inf,
                )
                next_location = np.unravel_index(np.argmax(extra_score), shape)
                if not np.isfinite(extra_score[next_location]):
                    break
                chosen.append(next_location)

            if not chosen:
                raise ValueError(f"No viable generation sites within connection range of {island}")
            points = [
                (float(east[location]), float(north[location])) for location in chosen
            ]
            cells = (
                list(voronoi_polygons(MultiPoint(points), extend_to=region.geometry.envelope).geoms)
                if len(points) > 1
                else [region.geometry]
            )
            for index, point in enumerate(points):
                cell = next(polygon for polygon in cells if polygon.covers(Point(point)))
                candidates.append(
                    {
                        "name": f"{island}_site_{index}",
                        "name_microgrid": island,
                        "geometry": cell.intersection(region.geometry),
                        "point": point,
                    }
                )
            logger.info("Selected %s candidate generation sites for %s", len(points), island)

    result = gpd.GeoDataFrame(candidates, geometry="geometry", crs=grid.crs)
    points = gpd.GeoSeries.from_xy(
        [entry["point"][0] for entry in candidates],
        [entry["point"][1] for entry in candidates],
        crs=grid.crs,
    ).to_crs("EPSG:4326")
    result["x"] = points.x.to_numpy()
    result["y"] = points.y.to_numpy()
    return result.drop(columns="point").to_crs("EPSG:4326")


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers_dist import mock_snakemake

        snakemake = mock_snakemake("dist_build_candidate_sites")
        sets_path_to_root("pypsa-distribution")

    configure_logging(snakemake)
    sites = build_candidate_sites(
        gpd.read_file(snakemake.input.regions),
        gpd.read_file(snakemake.input.clusters),
        snakemake.params.suitability,
        snakemake.params.settings,
    )
    sites.to_file(snakemake.output[0], driver="GeoJSON")