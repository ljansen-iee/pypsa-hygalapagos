# -*- coding: utf-8 -*-
"""
Estimates the population and the electric load of each microgrid.

Relevant Settings
-----------------

.. code:: yaml

    microgrids_list:
        microgridX: 
          lon_min:
          lon_max: 
          lat_min: 
          lat_max: 
    load:
        scaling_factor:

Inputs
------
- ``data/sample_profile.csv``: a load profile, which will be scaled through a scaling_factor to obtain the per person load

Outputs
-------
- ``resources/shapes/microgrid_shapes.geojson``: a geojson file of the shape of each microgrid,
- ``resources/masked_files/masked_file_{i+1}.tif``,
- ``resources/demand/microgrid_load_{i+1}.csv``: the electric load of the microgid,

Description
-----------
The rule :mod:`build_demand` contains functions that are used to create a shape file of the microgrid, to mask a raster with the shape file and to estimate 
the population. Then the population is multiplied for the per person load and the microgrid load is then obtained. The process applies to all the microgrids specified in config.yaml.
"""

import json
import logging
import os
import shutil

import geopandas as gpd
import numpy as np
import pandas as pd
import pypsa
import rasterio
import rasterio.mask
import requests
from _helpers_dist import (
    configure_logging,
    sets_path_to_root,
    two_2_three_digits_country,
)

_logger = logging.getLogger(__name__)
_logger.setLevel(logging.INFO)


def get_WorldPop_data(
    country_code,
    year,
    update=False,
    out_logging=False,
    size_min=300,
):
    """
    Download tiff file for each country code using the standard method from worldpop datastore with 1kmx1km resolution.

    Parameters
    ----------
    country_code : str
        Two letter country codes of the downloaded files.
        Files downloaded from https://data.worldpop.org/ datasets WorldPop UN adjusted
    year : int
        Year of the data to download
    update : bool
        Update = true, forces re-download of files
    size_min : int
        Minimum size of each file to download
    Returns
    -------
    WorldPop_inputfile : str
        Path of the file
    """

    three_digits_code = two_2_three_digits_country(country_code)

    if out_logging:
        _logger.info("Get WorldPop datasets")

    if country_code == "XK":
        WorldPop_filename = f"srb_ppp_{year}_UNadj_constrained.tif"
        WorldPop_urls = [
            f"https://data.worldpop.org/GIS/Population/Global_2000_2020_Constrained/2020/BSGM/SRB/{WorldPop_filename}",
            f"https://data.worldpop.org/GIS/Population/Global_2000_2020_Constrained/2020/maxar_v1/SRB/{WorldPop_filename}",
        ]
    else:
        WorldPop_filename = (
            f"{three_digits_code.lower()}_ppp_{year}_UNadj_constrained.tif"
        )
        # Urls used to possibly download the file
        WorldPop_urls = [
            f"https://data.worldpop.org/GIS/Population/Global_2000_2020_Constrained/2020/BSGM/{two_2_three_digits_country(country_code).upper()}/{WorldPop_filename}",
            f"https://data.worldpop.org/GIS/Population/Global_2000_2020_Constrained/2020/maxar_v1/{two_2_three_digits_country(country_code).upper()}/{WorldPop_filename}",
        ]

    WorldPop_inputfile = os.path.join(
        os.getcwd(),
        "pypsa-earth",
        "data",
        "WorldPop",
        WorldPop_filename,
    )  # Input filepath tif

    if not os.path.exists(WorldPop_inputfile) or update is True:
        if out_logging:
            _logger.warning(
                f"{WorldPop_filename} does not exist, downloading to {WorldPop_inputfile}"
            )
        #  create data/osm directory
        os.makedirs(os.path.dirname(WorldPop_inputfile), exist_ok=True)

        loaded = False
        for WorldPop_url in WorldPop_urls:
            with requests.get(WorldPop_url, stream=True) as r:
                with open(WorldPop_inputfile, "wb") as f:
                    if float(r.headers["Content-length"]) > size_min:
                        shutil.copyfileobj(r.raw, f)
                        loaded = True
                        break
        if not loaded:
            _logger.error(f"Impossible to download {WorldPop_filename}")

    return WorldPop_inputfile, WorldPop_filename


def estimate_microgrid_population(raster_path, shapes_path, output_file):
    """
    Estimates the population within each microgrid by using raster data and shape geometries.
    The function processes population density raster data and calculates the total population
    for each microgrid by masking the raster data using the corresponding geometries from a
    GeoJSON file. The population estimates are saved as a CSV file.

    Parameters
    ----------
    raster_path : str
        Path to the population density raster file (GeoTIFF format).
    shapes_path : str
        Path to the GeoJSON file containing the microgrid geometries.
    output_file : str
        Path to the CSV file where the population estimates will be saved.
    Returns
    -------
    pd.DataFrame
        A DataFrame containing the names of microgrids and their corresponding population estimates.
    """
    # Dictionary to store the population data for each microgrid
    population_data = {}
    # Load the GeoJSON file containing microgrid geometries
    shapes = gpd.read_file(shapes_path)
    # Iterate through each microgrid geometry
    for i, shape in shapes.iterrows():
        name = shape["name_microgrid"]  # Extract the name of the microgrid
        # Open the raster file and mask it using the microgrid geometry
        with rasterio.open(raster_path) as src:
            # Mask the raster data to only include the area within the microgrid
            masked, out_transform = rasterio.mask.mask(src, [shape.geometry], crop=True)
            # Update the raster metadata for the masked area
            out_meta = src.meta.copy()
            out_meta.update(
                {
                    "driver": "GTiff",
                    "height": masked.shape[1],
                    "width": masked.shape[2],
                    "transform": out_transform,
                }
            )
        # Calculate the total population within the microgrid by summing non-negative raster values
        pop_microgrid = masked[masked >= 0].sum()
        population_data[name] = pop_microgrid
    # Convert the population data dictionary to a DataFrame
    population_df = pd.DataFrame(
        list(population_data.items()), columns=["Microgrid_Name", "Population"]
    )

    return population_df


def calculate_load(
    p,
    raster_path,
    shapes_path,
    sample_profile,
    output_file,
    input_path,
    microgrids_list,
    start_date,
    end_date,
    inclusive,
):
    """
    Calculate the microgrid demand based on a load profile provided as input,
    appropriately scaled according to the population calculated for each cluster.
    The output includes a time-indexed DataFrame containing the load for each bus in the microgrid
    and is saved as a CSV file.

    Parameters
    ----------
    n : object
        PyPSA network object containing snapshots.
    p : int or float
        Scaling factor for the per-unit load.
    raster_path : str
        Path to the raster file containing population density data.
    shapes_path : str
        Path to the GeoJSON file containing the geometries of the microgrids.
    sample_profile : str
        Path to the CSV file containing the sample load profile.
    output_file : str
        Path where the resulting load profile CSV file will be saved.
    input_path : str
        Path to the CSV file containing building classifications.
    microgrids_list : dict
        Dictionary with microgrid names as keys and their cluster information as values.
    start_date : str
        Start date for filtering the time series data.
    end_date : str
        End date for filtering the time series data.
    inclusive : str
        Specifies whether the filtering is inclusive of the start or end date. Possible values: "left" or "right".

    Returns
    -------
    pd.DataFrame
        DataFrame containing the calculated load profile for all microgrids.

    """
    # Estimate the population for the two microgrids
    pop_microgrid = estimate_microgrid_population(raster_path, shapes_path, output_file)
    # Load the building classification data
    building_class = pd.read_csv(input_path)
    # Dictionary to store the load profiles for each microgrid
    microgrid_dataframes = {}
    # Load the sample load profile. It is a repeating annual shape (8760
    # hourly values), not tied to any specific calendar year, so it is
    # aligned to the requested snapshots by (month, day, hour) rather than
    # by literal date -- that way `snapshots.start/end` in the config can be
    # set to any year, not just the year the profile happens to be indexed as.
    df = pd.read_csv(sample_profile)
    per_unit_load_template = (df["0"] / p).to_numpy()  # Scale by the provided factor `p`
    template_index = pd.date_range(
        start="2001-01-01", periods=len(per_unit_load_template), freq="h"
    )
    per_unit_by_key = pd.Series(
        per_unit_load_template,
        index=[(t.month, t.day, t.hour) for t in template_index],
    )

    # Generate the snapshots range for filtering
    snapshots_range = pd.date_range(
        start=start_date, end=end_date, freq="h", inclusive=inclusive
    )

    # Reuse Feb 28's shape for Feb 29 in leap years, since the template
    # itself only covers a 365-day (non-leap) year.
    lookup_keys = [
        (t.month, 28 if (t.month, t.day) == (2, 29) else t.day, t.hour)
        for t in snapshots_range
    ]
    per_unit_load = per_unit_by_key.reindex(lookup_keys).to_numpy()

    # Loop over each microgrid
    for grid_name in microgrids_list.keys():
        microgrid_buildings = building_class.loc[
            building_class["name_microgrid"] == grid_name
        ]
        if microgrid_buildings.empty or microgrid_buildings["cluster_id"].isna().any():
            raise ValueError(
                f"Buildings in {grid_name} are missing assigned distribution buses; "
                "check brown-field voltage filtering and nearest-node assignment."
            )

        building_for_cluster = (
            microgrid_buildings.groupby("cluster_id")["count"]
            .sum()
            .astype(float)
        )
        total_buildings = building_for_cluster.sum()
        if total_buildings <= 0:
            raise ValueError(f"{grid_name} has no buildings with a positive count.")

        # Retrieve the population for the current microgrid
        population_rows = pop_microgrid.loc[
            pop_microgrid["Microgrid_Name"] == grid_name, "Population"
        ]
        if population_rows.empty:
            raise ValueError(f"No WorldPop population value found for {grid_name}.")
        pop_for_microgrid = float(population_rows.iloc[0])
        if not np.isfinite(pop_for_microgrid) or pop_for_microgrid <= 0:
            raise ValueError(
                f"{grid_name} has buildings but no positive WorldPop population; "
                "refusing to emit zero demand."
            )

        # Calculate the population per building and per cluster
        population_per_building = pop_for_microgrid / total_buildings
        population_per_cluster = building_for_cluster * population_per_building
        bus_names = [
            f"{grid_name}_bus_{int(bus_id) if float(bus_id).is_integer() else bus_id}"
            for bus_id in population_per_cluster.index
        ]
        load_per_cluster = pd.DataFrame(
            np.outer(
                per_unit_load,
                population_per_cluster.to_numpy(dtype=float),
            ),
            columns=bus_names,
            index=snapshots_range,
        )
        # Add the DataFrame for the microgrid to the dictionary
        microgrid_dataframes[grid_name] = load_per_cluster

    # Concatenate all microgrid DataFrames horizontally
    all_load_per_cluster = pd.concat(microgrid_dataframes.values(), axis=1)
    if all_load_per_cluster.empty or (all_load_per_cluster.sum(axis=0) <= 0).any():
        raise ValueError("Generated load profiles contain an empty or zero-demand bus.")
    if not np.isfinite(all_load_per_cluster.to_numpy(dtype=float)).all():
        raise ValueError("Generated load profiles contain non-finite values.")

    # Save the cumulative results to a CSV file with time index as the first column
    all_load_per_cluster.to_csv(output_file, index_label="Time")
    return all_load_per_cluster


def calculate_load_ramp(
    input_file_buildings,
    p,
    raster_path,
    shapes_path,
    sample_profile,
    output_file,
    input_file_profile_tier1,
    input_file_profile_tier2,
    input_file_profile_tier3,
    input_file_profile_tier4,
    input_file_profile_tier5,
    output_path_csv,
    tier_percent,
    date_start,
    date_end,
    inclusive,
    microgrid_list,
    std,
):
    # Upload of buildings and data demand for each tier
    cleaned_buildings = gpd.read_file(input_file_buildings)
    demand_files = [
        input_file_profile_tier1,
        input_file_profile_tier2,
        input_file_profile_tier3,
        input_file_profile_tier4,
        input_file_profile_tier5,
    ]

    mean_demand_tier_df = pd.DataFrame(
        {
            f"tier_{i+1}": pd.read_excel(file)["mean"]
            for i, file in enumerate(demand_files)
        }
    )
    std_demand_tier_df = pd.DataFrame(
        {
            f"tier_{i+1}": pd.read_excel(file)["std"]
            for i, file in enumerate(demand_files)
        }
    )
    mean_demand_tier_df.insert(0, "tier_0", np.zeros(len(mean_demand_tier_df)))
    std_demand_tier_df.insert(0, "tier_0", np.zeros(len(mean_demand_tier_df)))
    mean_demand_tier_df.index = pd.date_range(
        "00:00:00", periods=len(mean_demand_tier_df), freq="H"
    ).time
    std_demand_tier_df.index = pd.date_range(
        "00:00:00", periods=len(mean_demand_tier_df), freq="H"
    ).time

    pop = estimate_microgrid_population(raster_path, shapes_path, output_file)

    all_microgrid_loads = pd.DataFrame()

    for grid_name, grid_data in microgrid_list.items():
        microgrid_buildings = cleaned_buildings[
            cleaned_buildings["name_microgrid"] == grid_name
        ]
        # Calculate the population density for the current microgrid based only on house buildings
        house = microgrid_buildings[microgrid_buildings["tags_building"] == "house"]
        pop_microgrid = pop.loc[
            pop["Microgrid_Name"] == grid_name, "Population"
        ].values[0]
        density = pop_microgrid / house["area_m2"].sum()

        # Calculate population per cluster
        grouped_buildings = microgrid_buildings.groupby("cluster_id")
        clusters = np.sort(microgrid_buildings["cluster_id"].unique())
        house_area_for_cluster = [
            grouped_buildings.get_group(cluster)[
                grouped_buildings.get_group(cluster)["tags_building"] == "house"
            ]["area_m2"].sum()
            for cluster in clusters
        ]
        population_df = pd.DataFrame(
            {"cluster": clusters, "house_area_for_cluster": house_area_for_cluster}
        ).set_index("cluster")
        population_df["people_for_cluster"] = (
            population_df["house_area_for_cluster"] * density
        ).round()
        tier_pop_df = pd.DataFrame(
            np.outer(population_df["people_for_cluster"], tier_percent),
            index=population_df.index.astype(int),
        )

        if inclusive == "left":
            date_range = pd.date_range(start=date_start, end=date_end, freq="D")[:-1]
        else:
            date_range = pd.date_range(start=date_start, end=date_end, freq="D")

        mean_demand_tier_df_extended = pd.concat(
            [mean_demand_tier_df] * len(date_range), ignore_index=True
        )
        std_demand_tier_df_extended = pd.concat(
            [std_demand_tier_df] * len(date_range), ignore_index=True
        )

        # Calculate load for each cluster and tier
        if std == "on":
            result_dict = {}
            for k, pop_cluster in tier_pop_df.iterrows():
                load_df = pd.DataFrame()
                for j, n_person in enumerate(pop_cluster / 7):  # Scale by family size
                    mean_load = mean_demand_tier_df_extended.iloc[:, j] * n_person
                    std_load = np.random.normal(
                        mean_demand_tier_df_extended.iloc[:, j],
                        std_demand_tier_df_extended.iloc[:, j],
                    ) * np.sqrt(n_person)
                    total_load = (mean_load + std_load) / 1e6
                    load_df[f"tier_{j}"] = total_load
                result_dict[f"{grid_name}_bus_{k}"] = load_df
        elif std == "off":
            result_dict = {}
            for k, pop_cluster in tier_pop_df.iterrows():
                load_df = pd.DataFrame()
                for j, n_person in enumerate(pop_cluster / 7):  # Scale by family size
                    mean_load = mean_demand_tier_df_extended.iloc[:, j] * n_person
                    total_load = (mean_load) / 1e6
                    load_df[f"tier_{j}"] = total_load
                result_dict[f"{grid_name}_bus_{k}"] = load_df

        # Aggregate total load per cluster
        tot_result_dict = {
            f"{k}": df.sum(axis=1).rename(f"{k}") for k, df in result_dict.items()
        }
        tot_loads_df = pd.concat(tot_result_dict.values(), axis=1)
        if inclusive == "left":
            date_range_tot = pd.date_range(start=date_start, end=date_end, freq="H")[
                :-1
            ]
        else:
            date_range_tot = pd.date_range(start=date_start, end=date_end, freq="H")
        tot_loads_df.index = date_range_tot

        # Replace zero values with a small value just for avoid problem with plotting
        small_value = 1e-26
        tot_loads_df.loc[:, (tot_loads_df == 0).all()] = small_value

        all_microgrid_loads = pd.concat([all_microgrid_loads, tot_loads_df], axis=1)

    all_microgrid_loads.to_csv(output_path_csv)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers_dist import mock_snakemake

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake("build_demand")
        sets_path_to_root("pypsa-distribution")

    configure_logging(snakemake)

    sample_profile = snakemake.input["sample_profile"]
    tier_percent = snakemake.params.tier["tier_percent"]
    date_start = snakemake.params.snapshots["start"]
    date_end = snakemake.params.snapshots["end"]
    inclusive = snakemake.params.snapshots["inclusive"]
    microgrids_list = snakemake.config["microgrids_list"]
    build_demand_model = snakemake.params.build_demand_model["type"]
    std = snakemake.params.build_demand_model["std"]

    assert (
        len(snakemake.config["countries"]) == 1
    ), "Error: only a country shall be specified"

    worldpop_path, worldpop_flname = get_WorldPop_data(
        snakemake.config["countries"][
            0
        ],  # TODO: this needs fix to generalize the countries
        snakemake.config["build_shape_options"]["year"],
        False,
    )

    estimate_microgrid_population(
        worldpop_path,
        snakemake.input["microgrid_shapes"],
        snakemake.output["electric_load"],
    )
    if build_demand_model == "From_file":
        calculate_load(
            snakemake.config["load"]["scaling_factor"],
            worldpop_path,
            snakemake.input["microgrid_shapes"],
            sample_profile,
            snakemake.output["electric_load"],
            snakemake.input["building_csv"],
            microgrids_list,
            date_start,
            date_end,
            inclusive,
        )

    elif build_demand_model == "Ramp":
        calculate_load_ramp(
            snakemake.input["clusters_with_buildings"],
            snakemake.config["load"]["scaling_factor"],
            worldpop_path,
            snakemake.input["microgrid_shapes"],
            sample_profile,
            snakemake.output["electric_load"],
            snakemake.input["profile_Tier1"],
            snakemake.input["profile_Tier2"],
            snakemake.input["profile_Tier3"],
            snakemake.input["profile_Tier4"],
            snakemake.input["profile_Tier5"],
            snakemake.output["electric_load"],
            tier_percent,
            date_start,
            date_end,
            inclusive,
            microgrids_list,
            std,
        )
