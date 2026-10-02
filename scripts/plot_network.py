import argparse
from pathlib import Path

import cartopy.crs as ccrs
import geopandas as gpd
import matplotlib.pyplot as plt
import pandas as pd
import pypsa
from matplotlib.lines import Line2D
from pypsa.plot import add_legend_circles, add_legend_patches


ROOT = Path(__file__).resolve().parents[1]
COLORS = {
    "onwind": "#179C7D",
    "solar": "#F6A800",
    "diesel": "#ce6757",
    "battery": "#468ebd",
    "default": "#0065BD",
    "da20": "#39C1CD",
    "clow": "#A8508C",
    "da20clow": "#7C154D",
    "min": "#179C7D",
    "max": "#7C154D",
    "downstream_low": "#179C7D",
    "downstream_high": "#7C154D",
    "production": "#0065BD",
    "transport": "#39C1CD",
    "storage": "#A8508C",
    "distribution": "#F6A800",
}


def plot_network(
    network_path, buildings_path, output_path, basemap=False,
    country_borders_path=ROOT / "pypsa-earth/data/gadm/gadm41_ECU/gadm41_ECU.gpkg",
):
    network = pypsa.Network(network_path)
    buildings = gpd.read_file(buildings_path)
    if buildings.empty or buildings.crs is None:
        raise ValueError("Building polygons must have geometries and a CRS")

    ac_buses = network.buses.index[network.buses.carrier == "AC"]
    if ac_buses.empty or network.buses.loc[ac_buses, ["x", "y"]].isna().any().any():
        raise ValueError("The network needs AC buses with geographic coordinates")

    generators = network.generators.loc[
        network.generators.bus.isin(ac_buses)
        & (network.generators.p_nom_opt > 0)
        & (network.generators.p_nom_opt < 1e5),
        ["bus", "carrier", "p_nom_opt"],
    ]
    storage = network.storage_units.loc[
        network.storage_units.bus.isin(ac_buses)
        & (network.storage_units.p_nom_opt > 0),
        ["bus", "carrier", "p_nom_opt"],
    ].copy()
    storage.loc[storage.carrier.isin(["lithium", "lead acid"]), "carrier"] = "battery"
    capacities = pd.concat([generators, storage]).groupby(["bus", "carrier"])["p_nom_opt"].sum()
    if capacities.empty:
        raise ValueError("No positive generation or storage capacity at AC buses")
    unknown = capacities.index.get_level_values("carrier").difference(COLORS)
    if not unknown.empty:
        raise ValueError(f"Missing plot colors for carriers: {', '.join(unknown)}")

    buildings = buildings.to_crs(4326)
    building_west, building_south, building_east, building_north = buildings.total_bounds
    bus_coords = network.buses.loc[ac_buses, ["x", "y"]]
    west, south = min(building_west, bus_coords.x.min()), min(building_south, bus_coords.y.min())
    east, north = max(building_east, bus_coords.x.max()), max(building_north, bus_coords.y.max())
    pad_x, pad_y = (east - west) * 0.12, (north - south) * 0.12
    bounds = [west - pad_x, east + pad_x, south - pad_y, north + pad_y]
    projection = ccrs.Mercator(central_longitude=float(network.buses.loc[ac_buses, "x"].mean()))
    projected_buildings = buildings.to_crs(projection.to_string())
    country = gpd.read_file(country_borders_path, layer="ADM_ADM_0")
    if country.empty or country.crs is None:
        raise ValueError("Country border layer must contain geometries with a CRS")
    # Clip the boundary, not the polygon: clipping polygons creates false lines at the map edges.
    country_border = country.to_crs(4326).boundary.clip((bounds[0], bounds[2], bounds[1], bounds[3]))
    if country_border.empty:
        raise ValueError("Country border does not intersect the network map extent")
    projected_border = country_border.to_crs(projection.to_string())

    display_network = network.copy()
    display_network.remove("Link", display_network.links.index)
    display_network.remove("Bus", display_network.buses.index.difference(ac_buses))
    max_capacity = capacities.groupby(level="bus").sum().max()
    reference_area = 3.14159 * (0.025 * (north - south)) ** 2
    bus_sizes = capacities / max_capacity * reference_area
    line_widths = display_network.lines.s_nom_opt.clip(lower=0) / 20

    fig, ax = plt.subplots(figsize=(15, 8), subplot_kw={"projection": projection})
    try:
        ax.set_extent(bounds, crs=ccrs.PlateCarree())
        if basemap:
            import contextily as ctx

            ctx.add_basemap(ax, crs=projection.to_string(), source=ctx.providers.OpenStreetMap.HOT, alpha=0.6)
        display_network.plot.map(
            ax=ax, projection=projection, geomap=True, geomap_color=False,
            bus_size=bus_sizes, bus_color=COLORS, bus_alpha=0.9,
            line_width=line_widths, line_color="#303d40",
            boundaries=bounds,
        )
        projected_buildings.plot(
            ax=ax, column="cluster_id", cmap="tab20", alpha=0.7,
            edgecolor="#777777", linewidth=0.15, legend=False, aspect="auto", zorder=1,
        )
        projected_border.plot(
            ax=ax, color="#526b72", linewidth=0.9, aspect="auto", zorder=1.5,
        )

        unbuilt = display_network.lines.loc[
            display_network.lines.index.str.endswith("_connection")
            & (display_network.lines.s_nom_opt <= 1e-5)
        ]
        for line in unbuilt.itertuples():
            start, end = display_network.buses.loc[[line.bus0, line.bus1], ["x", "y"]].itertuples(index=False)
            ax.plot([start.x, end.x], [start.y, end.y], transform=ccrs.PlateCarree(),
                    color="#8c9694", linewidth=0.75, linestyle="--", zorder=2)
        pending_sites = display_network.buses.loc[unbuilt.bus0]
        ax.scatter(pending_sites.x, pending_sites.y, transform=ccrs.PlateCarree(),
                   s=25, facecolors="white", edgecolors="#4a6262", linewidths=1,
                   zorder=5)

        carriers = capacities.index.get_level_values("carrier").unique()
        add_legend_patches(
            ax, [COLORS[carrier] for carrier in carriers], list(carriers),
            legend_kw={"loc": "upper right", "title": "Installed capacity (MW)"},
        )
        ax.add_artist(ax.get_legend())
        add_legend_circles(
            ax, [reference_area, reference_area / 4],
            [f"{max_capacity:.0f} MW", f"{max_capacity / 4:.0f} MW"],
            patch_kw={"facecolor": "#728284", "edgecolor": "#394849", "alpha": 0.9},
            legend_kw={"loc": "lower right", "title": "Bus total"},
        )
        ax.add_artist(ax.get_legend())
        if not pending_sites.empty:
            ax.legend(handles=[
                Line2D([], [], color="#8c9694", linestyle="--", label="Candidate connection"),
                Line2D([], [], marker="o", linestyle="None", markerfacecolor="white",
                       markeredgecolor="#4a6262", label="Unbuilt site"),
            ], loc="lower left", title="Generation sites")

        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=250, bbox_inches="tight")
    finally:
        plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot installed AC network capacity over building clusters")
    parser.add_argument("--network", type=Path, default=ROOT / "networks/results/elec.nc")
    parser.add_argument("--buildings", type=Path, default=ROOT / "resources/buildings/cluster_with_buildings.geojson")
    parser.add_argument("--output", type=Path, default=ROOT / "resources/network_plot.png")
    parser.add_argument("--country-borders", type=Path, default=ROOT / "pypsa-earth/data/gadm/gadm41_ECU/gadm41_ECU.gpkg")
    parser.add_argument("--basemap", action="store_true", help="Download OpenStreetMap tiles")
    args = parser.parse_args()
    plot_network(args.network, args.buildings, args.output, args.basemap, args.country_borders)