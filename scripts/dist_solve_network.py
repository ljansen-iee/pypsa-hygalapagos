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
from _helpers_dist import configure_logging, sets_path_to_root


def _calculate_annuity(lifetime, discount_rate):
    if isinstance(discount_rate, pd.Series):
        return pd.Series(1.0 / lifetime, index=discount_rate.index).where(
            discount_rate == 0,
            discount_rate
            / (1.0 - 1.0 / (1.0 + discount_rate) ** lifetime),
        )
    if discount_rate > 0:
        return discount_rate / (1.0 - 1.0 / (1.0 + discount_rate) ** lifetime)
    return 1.0 / lifetime


def load_hydrogen_costs(tech_costs, cost_config, Nyears=1):
    costs = pd.read_csv(
        tech_costs,
        index_col=["technology", "year", "parameter"],
    )

    costs.loc[costs.unit.str.contains("/kW", na=False), "value"] *= 1e3
    costs.loc[costs.unit.str.contains("USD", na=False), "value"] *= cost_config[
        "USD2013_to_EUR2013"
    ]

    costs = (
        costs.loc[pd.IndexSlice[:, cost_config["year"], :], "value"]
        .unstack(level=2)
        .groupby("technology")
        .sum(min_count=1)
    )
    costs = costs.fillna(
        {
            "discount rate": cost_config["discountrate"],
            "FOM": 0,
            "VOM": 0,
            "fuel": 0,
            "efficiency": 1,
            "lifetime": 25,
        }
    )
    costs["capital_cost"] = (
        _calculate_annuity(costs["lifetime"], costs["discount rate"])
        + costs["FOM"] / 100.0
    ) * costs["investment"] * Nyears
    costs["marginal_cost"] = costs["VOM"] + costs["fuel"] / costs["efficiency"]

    return costs


def _find_cost_technology(costs, candidates):
    for candidate in candidates:
        if candidate in costs.index:
            return candidate
    raise KeyError(f"None of the cost technologies are available: {candidates}")


def add_hydrogen(n, costs):
    """Add local hydrogen production, reconversion, and storage."""
    ac_buses = n.buses.index[n.buses.carrier == "AC"]
    if ac_buses.empty:
        raise ValueError("Cannot add hydrogen without AC buses")

    for carrier in ("H2", "H2 Electrolysis", "H2 Fuel Cell", "H2 Store Tank"):
        if carrier not in n.carriers.index:
            n.add("Carrier", carrier)

    h2_buses = ac_buses + " H2"
    for ac_bus, h2_bus in zip(ac_buses, h2_buses):
        n.add(
            "Bus",
            h2_bus,
            carrier="H2",
            x=n.buses.at[ac_bus, "x"],
            y=n.buses.at[ac_bus, "y"],
        )

        n.add(
            "Link",
            f"{ac_bus} H2 Electrolysis",
            bus0=ac_bus,
            bus1=h2_bus,
            p_nom_extendable=True,
            carrier="H2 Electrolysis",
            efficiency=costs.at["electrolysis", "efficiency"],
            capital_cost=costs.at["electrolysis", "capital_cost"],
            lifetime=costs.at["electrolysis", "lifetime"],
        )

        n.add(
            "Link",
            f"{h2_bus} H2 Fuel Cell",
            bus0=h2_bus,
            bus1=ac_bus,
            p_nom_extendable=True,
            carrier="H2 Fuel Cell",
            efficiency=costs.at["fuel cell", "efficiency"],
            capital_cost=costs.at["fuel cell", "capital_cost"]
            * costs.at["fuel cell", "efficiency"],
            lifetime=costs.at["fuel cell", "lifetime"],
        )

    storage_technology = _find_cost_technology(
        costs,
        (
            "hydrogen storage tank type 1 including compressor",
            "hydrogen storage tank",
            "hydrogen storage",
        ),
    )
    for h2_bus in h2_buses:
        n.add(
            "Store",
            f"{h2_bus} Store Tank",
            bus=h2_bus,
            e_nom_extendable=True,
            e_cyclic=True,
            carrier="H2 Store Tank",
            capital_cost=costs.at[storage_technology, "capital_cost"],
        )

    return n


def prepare_network(n, solve_opts):

    if "clip_p_max_pu" in solve_opts:
        for df in (n.generators_t.p_max_pu, n.storage_units_t.inflow):
            df.where(df > solve_opts["clip_p_max_pu"], other=0.0, inplace=True)

    load_shedding = solve_opts.get("load_shedding")
    if load_shedding:
        n.add("Carrier", "Load")
        buses_i = n.buses.query("carrier == 'AC'").index
        if not np.isscalar(load_shedding):
            load_shedding = 8e3  # Eur/kWh
        # intersect between macroeconomic and surveybased
        # willingness to pay
        # http://journal.frontiersin.org/article/10.3389/fenrg.2015.00055/full)
        # 1e2 is practical relevant, 8e3 good for debugging
        n.madd(
            "Generator",
            buses_i,
            " load",
            bus=buses_i,
            carrier="load",
            sign=1e-3,  # Adjust sign to measure p and p_nom in kW instead of MW
            marginal_cost=load_shedding,
            p_nom=1e9,  # kW
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


def solve_network(n, solver_name, **solver_options):
    n.optimize(solver_name=solver_name)

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
    costs = load_hydrogen_costs(
        snakemake.input["tech_costs"],
        snakemake.config["costs"],
        Nyears,
    )
    n = add_hydrogen(n, costs)

    n = prepare_network(n, snakemake.config["solving"]["options"])
    n = fixing_missing_carriers(n)

    n = solve_network(n, solver_name, **solver_options)

    n.export_to_netcdf(snakemake.output[0])
