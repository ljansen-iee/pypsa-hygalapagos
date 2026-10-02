# -*- coding: utf-8 -*-
"""
Solves linear optimal power flow for a network iteratively.

.. code:: yaml

    solving:
        tmpdir:
        options:
            formulation:
            clip_p_max_pu:
            load_shedding:
            noisy_costs:
            nhours:
            min_iterations:
            max_iterations:
            skip_iterations:
            track_iterations:
        
Inputs
------
- ``networks/elec.nc``

Outputs
-------
- ``networks/results/networks/elec.nc``: Solved PyPSA network including optimisation results
  
Description
-----------
Total annual system costs are minimised with PyPSA. The full formulation of the
linear optimal power flow (plus investment planning
is provided in the
`documentation of PyPSA <https://pypsa.readthedocs.io/en/latest/optimal_power_flow.html#linear-optimal-power-flow>`_.
The optimization is based on the ``pyomo=False`` setting in the :func:`network.lopf` function.
Additionally, some extra constraints specified in :mod:`prepare_network` are added.
"""

import os

import numpy as np
import pandas as pd
import pypsa
from pypsa.descriptors import get_switchable_as_dense as get_as_dense
from _helpers_dist import configure_logging, pe_helpers, sets_path_to_root


def load_hydrogen_costs(tech_costs, cost_config, Nyears=1):
    """Read and annualize one PyPSA-Earth technology-cost file."""
    return pe_helpers.prepare_costs(
        tech_costs,
        cost_config,
        cost_config["output_currency"],
        cost_config["fill_values"],
        Nyears,
        cost_config["default_exchange_rate"],
        cost_config["future_exchange_rate_strategy"],
        cost_config["custom_future_exchange_rate"],
    )


def add_hydrogen(n, costs):
    """Add local hydrogen production, reconversion, and storage."""
    ac_buses = n.buses.index[n.buses.carrier == "AC"]
    if ac_buses.empty:
        raise ValueError("Cannot add hydrogen without AC buses")

    for carrier in ("H2", "H2 Electrolysis", "H2 Fuel Cell", "H2 Store Tank"):
        if carrier not in n.carriers.index:
            n.add("Carrier", carrier)

    n.madd(
        "Bus",
        ac_buses + " H2",
        location=ac_buses,
        carrier="H2",
        x=n.buses.loc[ac_buses, "x"].to_numpy(),
        y=n.buses.loc[ac_buses, "y"].to_numpy(),
    )
    n.madd(
        "Link",
        ac_buses + " H2 Electrolysis",
        bus0=ac_buses,
        bus1=ac_buses + " H2",
        p_nom_extendable=True,
        carrier="H2 Electrolysis",
        efficiency=1 / costs.at["Alkaline electrolyzer large size", "electricity-input"],
        capital_cost=costs.at["Alkaline electrolyzer large size", "fixed"],
        lifetime=costs.at["Alkaline electrolyzer large size", "lifetime"],
    )
    n.madd(
        "Link",
        ac_buses + " H2 turbine",
        bus0=ac_buses + " H2",
        bus1=ac_buses,
        p_nom_extendable=True,
        carrier="H2 Fuel Cell",
        efficiency=costs.at["OCGT", "efficiency"],
        capital_cost=costs.at["OCGT", "fixed"] * costs.at["OCGT", "efficiency"],
        marginal_cost=costs.at["OCGT", "VOM"],
        lifetime=costs.at["OCGT", "lifetime"],
    )
    n.madd(
        "Store",
        ac_buses + " H2 Store Tank",
        bus=ac_buses + " H2",
        e_nom_extendable=True,
        e_cyclic=True,
        carrier="H2 Store Tank",
        capital_cost=costs.at["Hydrogen-store", "fixed"],
        lifetime=costs.at["Hydrogen-store", "lifetime"],
    )

    return n


def add_methanol(n, costs):
    """Add local methanol assets and a flat 4,000 MWh/year market load."""
    ac_buses = n.buses.index[n.buses.carrier == "AC"]
    for carrier in ("methanol", "methanolisation", "CO2 feedstock"):
        if carrier not in n.carriers.index:
            n.add("Carrier", carrier)

    n.madd(
        "Bus",
        ac_buses + " methanol",
        location=ac_buses,
        carrier="methanol",
        unit="MWh_LHV",
        x=n.buses.loc[ac_buses, "x"].to_numpy(),
        y=n.buses.loc[ac_buses, "y"].to_numpy(),
    )
    methanol_energy_density = 5.54 * 791 * 1e-3  # MWh/m3
    n.madd(
        "Store",
        ac_buses + " methanol Store",
        bus=ac_buses + " methanol",
        e_nom_extendable=True,
        e_cyclic=True,
        carrier="methanol",
        capital_cost=(
            costs.at["General liquid hydrocarbon storage (product)", "fixed"]
            / methanol_energy_density
        ),
    )
    n.madd(
        "Bus",
        ac_buses + " CO2 feedstock",
        location=ac_buses,
        carrier="CO2 feedstock",
        unit="tCO2",
        x=n.buses.loc[ac_buses, "x"].to_numpy(),
        y=n.buses.loc[ac_buses, "y"].to_numpy(),
    )
    n.madd(
        "Generator",
        ac_buses + " CO2 supply",
        bus=ac_buses + " CO2 feedstock",
        carrier="CO2 feedstock",
        p_nom=np.inf,
        marginal_cost=80.0,
    )

    hydrogen_input = costs.at["methanolisation", "hydrogen-input"]
    efficiency = 1 / hydrogen_input
    n.madd(
        "Link",
        ac_buses + " methanolisation",
        bus0=ac_buses + " H2",
        bus1=ac_buses + " methanol",
        bus2=ac_buses + " CO2 feedstock",
        bus3=ac_buses,
        carrier="methanolisation",
        p_nom_extendable=True,
        efficiency=efficiency,
        efficiency2=-costs.at["methanol", "CO2 intensity"],
        efficiency3=-costs.at["methanolisation", "electricity-input"] / hydrogen_input,
        capital_cost=costs.at["methanolisation", "fixed"] / hydrogen_input,
        marginal_cost=costs.at["methanolisation", "VOM"] / hydrogen_input,
        lifetime=costs.at["methanolisation", "lifetime"],
    )

    hours = n.snapshot_weightings.generators.sum()
    if not np.isfinite(hours) or hours <= 0:
        raise ValueError("Methanol demand requires positive snapshot generator weightings")
    n.madd(
        "Load",
        ac_buses + " methanol market",
        bus=ac_buses + " methanol",
        carrier="methanol",
        p_set=4000 / (hours * len(ac_buses)),
    )
    return n


def prepare_network(n, solve_opts):

    if "clip_p_max_pu" in solve_opts:
        for df in (n.generators_t.p_max_pu, n.storage_units_t.inflow):
            df.where(df > solve_opts["clip_p_max_pu"], other=0.0, inplace=True)

    def get_load_shedding_capacity(n, safety_margin=1.2):
        """
        Calculate required load shedding p_nom per bus based on the
        maximum aggregated load observed in any snapshot.

        Parameters
        ----------
        n : pypsa.Network
            The PyPSA network
        safety_margin : float, default 1.2
            Safety factor to apply to the maximum load

        Returns
        -------
        pd.Series
            Required p_nom per bus for load shedding.
        """

        load_profiles = get_as_dense(n, "Load", "p_set")
        load_by_bus = load_profiles.T.groupby(n.loads.bus).sum().T
        co2_buses = n.buses.index[n.buses.carrier.isin(["co2", "co2 stored"])]
        load_by_bus = load_by_bus.drop(columns=co2_buses, errors="ignore")

        # Load shedding can cover positive demand only. Negative-only buses, such
        # as process-emission or CO2 buses, must not receive negative capacities.
        load_shedding_p_nom = (
            load_by_bus.max(axis=0).clip(lower=0.0) * safety_margin
        )

        return load_shedding_p_nom.reindex(n.buses.index, fill_value=0.0)


    if solve_opts.get("load_shedding"):
        required_p_nom = get_load_shedding_capacity(n, safety_margin=1.2)
        n.add("Carrier", "load shedding", color="#dd2e23", nice_name="Load shedding")

        load_shedding_buses = n.buses.index[
            ~n.buses.carrier.isin(["co2", "co2 stored"])
        ]
        n.madd(
            "Generator",
            load_shedding_buses,
            " load shedding",
            bus=load_shedding_buses,
            carrier="load shedding",
            sign=1,
            marginal_cost=solve_opts.get("load_shedding") * 1000,
            p_nom=required_p_nom.reindex(load_shedding_buses, fill_value=0.5e6),
        )
    return n


def fixing_missing_carriers(n):
    """
    Fixing missing carriers for each component in the network
    to prevent pypsa.consistency warning.
    """

    components = ["lines", "generators", "buses", "storage_units"]
    all_carriers = set()

    for comp in components:
        df = getattr(n, comp)
        if "carrier" in df.columns:
            carriers = df["carrier"].dropna().unique()
            all_carriers.update(carriers)

    existing_carriers = set(n.carriers.index)
    missing_carriers = all_carriers - existing_carriers

    for carrier in missing_carriers:
        n.add("Carrier", carrier)

    return n


def solve_network(n, solver_name, capacity_targets=None, **solver_options):
    capacity_targets = capacity_targets or {}

    def constrain_capacity(network, snapshots):
        capacities = network.model["Generator-p_nom"]
        for carrier, target in capacity_targets.items():
            generators = network.generators.index[
                (network.generators.carrier == carrier)
                & network.generators.p_nom_extendable
            ]
            fixed = network.generators.loc[
                (network.generators.carrier == carrier)
                & ~network.generators.p_nom_extendable,
                "p_nom",
            ].sum()
            if generators.empty:
                raise ValueError(f"No extendable generators for capacity target: {carrier}")
            network.model.add_constraints(
                capacities.loc[generators].sum() == target - fixed,
                name=f"{carrier}_capacity_target",
            )

    status, condition = n.optimize(
        solver_name=solver_name,
        extra_functionality=constrain_capacity if capacity_targets else None,
        **solver_options,
    )
    if (status, condition) != ("ok", "optimal"):
        raise RuntimeError(f"Network optimization failed: {status}, {condition}")

    return n


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers_dist import mock_snakemake

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        snakemake = mock_snakemake("solve_network")
        sets_path_to_root("pypsa-distribution")

    configure_logging(snakemake)

    solver_options = snakemake.config["solving"]["solver"].copy()
    solver_name = solver_options.pop("name")

    n = pypsa.Network(snakemake.input[0])
    Nyears = n.snapshot_weightings.objective.sum() / 8760.0
    cost_config = snakemake.config["costs"].copy()
    costs = load_hydrogen_costs(
        snakemake.input["tech_costs"],
        cost_config,
        Nyears,
    )
    n = add_hydrogen(n, costs)
    n = add_methanol(n, costs)

    n = prepare_network(n, snakemake.config["solving"]["options"])
    n = fixing_missing_carriers(n)

    n = solve_network(
        n, solver_name,
        capacity_targets=snakemake.config["solving"].get("capacity_targets"),
        **solver_options,
    )

    n.export_to_netcdf(snakemake.output[0])
