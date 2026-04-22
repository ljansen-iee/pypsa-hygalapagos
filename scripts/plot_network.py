import cartopy.crs as ccrs
import contextily as ctx
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pypsa

# ── Load data ─────────────────────────────────────────────────────────────────
n = pypsa.Network("../networks/results/elec.nc")
buildings = gpd.read_file(
    "../resources/buildings/cluster_with_buildings.geojson"
)

# ── Projection (Mercator centred on network) ───────────────────────────────────
crs = ccrs.Mercator(central_longitude=n.buses.x.mean())

# Reproject buildings into the axes' Mercator CRS
buildings_merc = buildings.to_crs(crs.to_string())

# ── Carrier colours ───────────────────────────────────────────────────────────
carrier_colors = {
    "onwind":  "#74c476",
    "solar":   "#fdd835",
    "diesel":  "#e57373",
    "load":    "#bdbdbd",
    "lithium": "#64b5f6",
    "line":    "#444444",
}

# ── Bus sizes: total p_nom_opt per bus (generators + storage) ─────────────────
# Drop load-shedding / slack generators (p_nom_opt > 1e5 are placeholders)
real_gens = n.generators[n.generators["p_nom_opt"] <= 1e5]
gen_cap = (
    real_gens.groupby("bus")["p_nom_opt"].sum()
    .reindex(n.buses.index, fill_value=0)
)
sto_cap = (
    n.storage_units.groupby("bus")["p_nom_opt"].sum()
    .reindex(n.buses.index, fill_value=0)
)
total_cap = gen_cap + sto_cap
# bus_size must be in degrees² — PyPSA internally multiplies by area_factor²
# (~111km/deg)² to convert to Mercator metres² when a projection is used.
# Largest bus → circle radius = 2.5 % of the map's N-S extent in degrees.
_dy = n.buses.y.max() - n.buses.y.min()   # degrees
_r_max = 0.025 * _dy                       # degrees
_area_max = 3.14159 * _r_max ** 2          # degrees²  (~1.5 km radius after projection)
bus_sizes = total_cap / total_cap.max() * _area_max

# ── Line widths: s_nom_opt ────────────────────────────────────────────────────
line_width = (n.lines["s_nom_opt"] / n.lines["s_nom_opt"].max() * 2).clip(lower=0.3)

# ── Figure ────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(15, 15), subplot_kw={"projection": crs})

# # Set axes extent from buildings bounds so contextily infers the correct zoom.
# # A small padding (5 %) is added on all sides.
_bounds = buildings_merc.total_bounds  # [minx, miny, maxx, maxy] in metres
_pad_x = (_bounds[2] - _bounds[0]) * 0.12
_pad_y = (_bounds[3] - _bounds[1]) * 0.12
ax.set_extent(
    [_bounds[0] - _pad_x, _bounds[2] + _pad_x,
     _bounds[1] - _pad_y, _bounds[3] + _pad_y],
    crs=crs,
)

# ── Contextily basemap ────────────────────────────────────────────────────────
contextily_opts = {
    "add_basemap": True,
    "basemap_source": ctx.providers.OpenStreetMap.HOT,  # alternatives: Esri.WorldImagery (satellite), OpenStreetMap.Mapnik
    "basemap_alpha": 0.6,
    "basemap_zoom": "auto",  # extent is now set → auto resolves to a valid level
}
if contextily_opts["add_basemap"]:
    ctx.add_basemap(
        ax,
        crs=crs.to_string(),
        source=contextily_opts["basemap_source"],
        alpha=contextily_opts["basemap_alpha"],
        zoom=contextily_opts["basemap_zoom"],
    )



with plt.rc_context({"patch.linewidth": 0.3}):
    n.plot(
        ax=ax,
        projection=crs,
        geomap_color=False,
        # geomap_resolution="10m",
        bus_size=bus_sizes,
        bus_color=carrier_colors["diesel"],
        bus_alpha=0.85,
        line_widths=line_width,
        line_colors=carrier_colors["line"],
        line_alpha=0.8,
    )

buildings_merc = buildings.to_crs(crs.to_string())
# Background: buildings – reproject first, then plot with transform=crs
# (same pattern as regions.plot() in the pypsa-earth example)
buildings_merc.plot(
    ax=ax,
    column="cluster_id",
    cmap="tab20",
    alpha=0.95,
    edgecolor="orange",
    linewidth=0.8,
    legend=False,
    aspect="auto",
    transform=crs,
    zorder=1,
)
# ── Legend ────────────────────────────────────────────────────────────────────
legend_handles = [
    mpatches.Patch(color=carrier_colors["onwind"],  label="Wind capacity"),
    mpatches.Patch(color=carrier_colors["solar"],   label="Solar capacity"),
    mpatches.Patch(color=carrier_colors["diesel"],  label="Diesel / Bus node"),
    mpatches.Patch(color=carrier_colors["lithium"], label="Battery storage"),
    mpatches.Patch(color=carrier_colors["line"],    alpha=0.8, label="Lines (width ∝ s_nom_opt)"),
]
ax.legend(handles=legend_handles, loc="upper right", fontsize=8, framealpha=0.9)

plt.tight_layout()
plt.savefig("../resources/network_plot.png", dpi=500, bbox_inches="tight")
plt.show()


