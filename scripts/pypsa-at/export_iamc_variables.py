"""
Export variables in IAMC data model for all regions.

IAMC variable naming convention:
Category|Subcategory|Specification

Notes
-----
https://docs.ece.iiasa.ac.at/iamc.html
https://pyam-iamc.readthedocs.io/en/stable/
"""

import logging

import pandas as pd
from pyam import IamDataFrame
from pypsa.statistics import port_efficiency

from evals.constants import DataModel as DM
from evals.constants import TradeTypes
from evals.fileio import read_networks
from evals.statistic import collect_myopic_statistics
from evals.utils import (
    calculate_input_share,
    filter_by,
    filter_for_carrier_connected_to,
    insert_index_level,
    rename_aggregate,
)
from scripts._helpers import configure_logging, mock_snakemake

logger = logging.getLogger()

IDX = [DM.YEAR, DM.LOCATION, "unit"]


def _extract(ds: pd.Series, **filter_kwargs) -> pd.Series:
    """Extract and group filter results."""
    results = filter_by(ds, **filter_kwargs)
    ds.drop(results.index, inplace=True)
    return results.groupby(IDX).sum()
    # return filter_by(ds, drop=True, **filter_kwargs).groupby(IDX).sum()


def _get_traded_energy(n, var, bus_carrier, direction, subcat):
    """
    Calculate the trade statistics.

    Parameters
    ----------
    n
    var
    bus_carrier
    direction
    subcat

    Returns
    -------
    :
    """
    for scope in (TradeTypes.DOMESTIC, TradeTypes.FOREIGN):
        trade = n.statistics.trade_energy(
            bus_carrier=bus_carrier, direction=direction, scope=scope
        )
        trade = trade[trade.gt(0)].groupby("location").sum()
        var[f"Primary Energy|{subcat}|{direction.title()} {scope.title()}"] = trade


def _sum_by_subcategory(var, subcat):
    return (
        pd.concat([var[v] for v in var.keys() if f"|{subcat}|" in v]).groupby(IDX).sum()
    )


def primary_oil(var: dict) -> dict:
    """
    Calculate the amounts of oil entering a region.

    Assuming that regional oil production is consumed in the same region, the netted
    energy balance becomes the primary energy (import) and final energy (export) of
    the region. It is important to note that this is not the same as fossil oil imports,
    but all oil imports including renewable oil production from other regions.
    Lacking an oil network that connects regions, we assume that all oil is consumed
    locally and only oil production surplus becomes exported as final energy.


    There are a few caveats here:
    1. the EU oil bus mixes fossil oil and renewable liquids
    2. It is not possible to distinguish oil types anymore
       once they enter the EU bus
    3. Regional transformation Links are connected to the EU bus
    4. Anything other than "oil import" is considered "liquids"
    5. liquids may, or may not be renewable, depending on the input
       of the transformation technology
    6. We assume, that all liquids produced in a region are consumed in that region
       (This is why the energy_balance is calculated instead of the supply)
    6. the same consumption split (oil vs liquids) is used for all regions,
       because the regional split cannot be calculated.
    7. All oil imported to EU is fossil oil

    Parameters
    ----------
    var
        The collection to append new variables in.

    Returns
    -------
    :
        The updated variables' collection.
    """
    # # netted_liquids = n.statistics.energy_balance(
    # #     groupby="location", bus_carrier="oil", comps="Link"
    # # ).drop("EU", errors="ignore")
    # netted_liquids = collect_myopic_statistics(networks, "energy_balance", bus_carrier="oil", comps="Link").drop("EU", level="location")
    # liquids_consumption = netted_liquids[netted_liquids.le(0)].abs()

    # increase fossil oil demand by this factor to account for refining losses
    _refining_efficiencies = set()
    for n in networks.values():
        _refining_efficiencies.add(
            port_efficiency(n, "Link", "1").filter(like="oil refining").item()
        )
    if len(_refining_efficiencies) != 1:
        raise NotImplementedError("Multiple efficiencies not supported.")
    else:
        oil_refining_eff = _refining_efficiencies.pop()

    # assuming that all local oil production. Let's not filter_by components, to capture all but Stores.
    production = (
        filter_by(SUPPLY, bus_carrier="oil")
        .drop("Store", level="component")
        .groupby(IDX)
        .sum()
    )
    consumption = (
        filter_by(DEMAND, bus_carrier="oil")
        .drop("Store", level="component")
        .groupby(IDX)
        .sum()
    )
    regional_oil_imports = consumption.sub(production, fill_value=0).clip(lower=0)

    # assert primary oil and imports are equal per year
    eu_oil_import = filter_by(SUPPLY, carrier=["oil primary", "import oil"])
    try:
        pd.testing.assert_series_equal(
            regional_oil_imports.groupby("year").sum(),
            eu_oil_import.groupby("year").sum(),
            check_names=False,
        )
    except AssertionError:
        logger.warning("Oil amounts mismatch.")

    var["Primary Energy|Oil|Fossil Oil"] = (
        regional_oil_imports / oil_refining_eff
    )  # total amounts before refining
    var["Primary Energy|Oil|Refining Losses"] = regional_oil_imports * (
        1 - oil_refining_eff
    )
    var["Primary Energy|Oil|Global Import"] = _extract(SUPPLY, carrier="import oil")
    var["Primary Energy|Oil|Primary"] = _extract(SUPPLY, carrier="oil primary")
    var["Primary Energy|Oil"] = regional_oil_imports

    # var["Primary Energy|Oil|Fossil"] = oil  # including losses
    # var["Primary Energy|Oil|Liquids"] = liquids  # this is secondary energy
    # var["Primary Energy|Oil|Global Import"] = _extract(SUPPLY, carrier="import oil")
    # var["Primary Energy|Oil"] = oil + liquids
    # var["Primary Energy|Oil|Refining Losses"] = oil * (1 - oil_refining_eff)

    # # calculate the split of fossil oil generation to liquids production and assume
    # # the same split for all regions
    # total_oil_import = n.statistics.supply(bus_carrier="oil primary").item()
    # total_liquids_production = n.statistics.supply(
    #     groupby="bus_carrier", bus_carrier="oil", comps="Link"
    # ).item()
    # fossil_share = total_oil_import / total_liquids_production
    # oil = liquids_consumption * fossil_share / oil_refining_eff
    #
    # liquids = liquids_consumption * (1 - fossil_share)
    #
    # var["Primary Energy|Oil|Fossil"] = oil  # including losses
    # var["Primary Energy|Oil|Liquids"] = liquids  # this is secondary energy
    # var["Primary Energy|Oil|Global Import"] = _extract(SUPPLY, carrier="import oil")
    # var["Primary Energy|Oil"] = oil + liquids
    # var["Primary Energy|Oil|Refining Losses"] = oil * (1 - oil_refining_eff)

    return var


def primary_gas(var) -> dict:
    var["Primary Energy|Gas|Import Foreign"] = _extract(
        IMPORT_FOREIGN, bus_carrier="gas"
    )
    var["Primary Energy|Gas|Import Domestic"] = _extract(
        IMPORT_DOMESTIC, bus_carrier="gas"
    )
    var["Primary Energy|Gas|Production"] = _extract(
        SUPPLY, carrier="gas", bus_carrier="gas", component="Generator"
    )
    var["Primary Energy|Gas|Import Global"] = _extract(
        SUPPLY,
        carrier="import gas",
        bus_carrier="gas",
        component="Generator",
    )
    var["Primary Energy|Gas"] = _sum_by_subcategory(var, "Gas")

    return var


def primary_waste(var: dict) -> dict:
    """

    Parameters
    ----------
    var

    Returns
    -------
    :
    """
    bus_carrier = "municipal solid waste"
    var["Primary Energy|Waste|Import Foreign"] = _extract(
        IMPORT_FOREIGN, bus_carrier=bus_carrier
    )
    var["Primary Energy|Waste|Import Domestic"] = _extract(
        IMPORT_DOMESTIC, bus_carrier=bus_carrier
    )
    var["Primary Energy|Waste|Solid"] = _extract(
        SUPPLY, bus_carrier=bus_carrier, component="Generator"
    )
    var["Primary Energy|Waste"] = _sum_by_subcategory(var, "Waste")

    return var


def primary_coal(var: dict) -> dict:
    """
    Calculate the amounts of coal consumed in a region.

    Coal is not produced by any Link, therefore it's safe to assume
    all Link withdrawal is imported fossil coal or lignite.

    Parameters
    ----------
    var
        The collection to append new variables in.

    Returns
    -------
    :
        The updated variables' collection.
    """
    var["Primary Energy|Coal|Hard"] = (
        filter_by(DEMAND, bus_carrier="coal", component="Link").groupby(IDX).sum()
    )
    var["Primary Energy|Coal|Lignite"] = (
        filter_by(DEMAND, bus_carrier="lignite", component="Link").groupby(IDX).sum()
    )
    var["Primary Energy|Coal"] = _sum_by_subcategory(var, "Coal")

    # remove EU coal generators from the to-do list
    coal_generators = filter_by(
        SUPPLY, bus_carrier=["coal", "lignite"], component="Generator"
    )
    SUPPLY.drop(coal_generators, inplace=True)

    return var


def primary_hydrogen(var: dict) -> dict:
    """
    Calculate the amounts of hydrogen imported into a region.

    There are global import of Hydrogen, found in `Generator`
    components, and various types of H2 pipelines, that bring
    Hydrogen into regions.

    Parameters
    ----------
    var
        The collection to append new variables in.

    Returns
    -------
    :
        The updated variables' collection.
    """
    bus_carrier = "H2"
    var["Primary Energy|Hydrogen|Import Foreign"] = _extract(
        IMPORT_FOREIGN, bus_carrier=bus_carrier
    )
    var["Primary Energy|Hydrogen|Import Domestic"] = _extract(
        IMPORT_DOMESTIC, bus_carrier=bus_carrier
    )
    var["Primary Energy|Hydrogen|Import Global"] = _extract(
        SUPPLY, carrier="import H2", bus_carrier=bus_carrier, component="Generator"
    )

    var["Primary Energy|Hydrogen"] = _sum_by_subcategory(var, "Hydrogen")

    return var


def primary_biomass(var: dict) -> dict:
    """
    Calculate the amounts of biomass generated in a region.

    Parameters
    ----------
    var
        The collection to append new variables in.

    Returns
    -------
    :
        The updated variables' collection.
    """
    bus_carrier = "solid biomass"
    var["Primary Energy|Biomass|Import Foreign"] = _extract(
        IMPORT_FOREIGN, bus_carrier=bus_carrier
    )
    var["Primary Energy|Biomass|Import Domestic"] = _extract(
        IMPORT_DOMESTIC, bus_carrier=bus_carrier
    )

    var["Primary Energy|Biomass|Solid"] = _extract(
        SUPPLY, bus_carrier=bus_carrier, component="Generator"
    )
    var["Primary Energy|Biomass|Biogas"] = _extract(
        SUPPLY, bus_carrier="biogas", component="Generator"
    )
    var["Primary Energy|Biomass"] = _sum_by_subcategory(var, "Biomass")

    return var


def primary_hydro(var: dict) -> dict:
    """
    Calculate the hydropower generated per region.

    The inflow split is kept, in case PHS receives an update
    that supplies inflow amounts to it.

    Parameters
    ----------
    var
        The collection to append new variables in.

    Returns
    -------
    :
        The updated variables' collection.
    """
    var["Primary Energy|Hydro|PHS"] = _extract(
        SUPPLY, carrier="PHS", component="StorageUnit"
    )
    var["Primary Energy|Hydro|Reservoir"] = _extract(
        SUPPLY, carrier="hydro", component="StorageUnit"
    )
    var["Primary Energy|Hydro|Run-of-River"] = _extract(
        SUPPLY, carrier="ror", component="Generator"
    )
    var["Primary Energy|Hydro"] = _sum_by_subcategory(var, "Hydro")

    return var


def primary_solar(var: dict) -> dict:
    """
    Calculate solar energy generated per region.

    Parameters
    ----------
    n
        A solved network.
    var
        The collection to append new variables in.

    Returns
    -------
    :
        The updated variables' collection.
    """
    var["Primary Energy|Solar|Utility"] = _extract(
        SUPPLY, carrier="solar", component="Generator"
    )
    var["Primary Energy|Solar|HSAT"] = _extract(
        SUPPLY, carrier="solar-hsat", component="Generator"
    )
    var["Primary Energy|Solar|Rooftop"] = _extract(
        SUPPLY, carrier="solar rooftop", component="Generator"
    )
    var["Primary Energy|Solar"] = _sum_by_subcategory(var, "Solar")

    return var


def primary_liquids(var: dict) -> dict:
    var["Primary Energy|Liquids|Unsustainable Bioliquids"] = _extract(
        SUPPLY, carrier="import oil"
    )

    return var


def primary_nuclear(var: dict) -> dict:
    """
    Calculate the uranium demand for nuclear power plants per region.

    Parameters
    ----------
    var
        The collection to append new variables in.

    Returns
    -------
    :
        The updated variables' collection.
    """
    var["Primary Energy|Nuclear|Uranium"] = (
        filter_by(DEMAND, carrier="uranium", component="Link").groupby(IDX).sum()
    )
    # var["Primary Energy|Nuclear|Electricity"] = _extract(SUPPLY, carrier="nuclear", component="Link")  # is secondary energy
    # remove EU uranium generators from the to-do list
    uranium_generators = filter_by(
        SUPPLY, bus_carrier=["coal", "lignite"], component="Generator"
    )
    SUPPLY.drop(uranium_generators, inplace=True)

    return var


def primary_ammonia(var: dict) -> dict:
    """
    Calculate the ammonium imported per region.

    Parameters
    ----------
    n
        A solved network.
    var
        The collection to append new variables in.

    Returns
    -------
    :
        The updated variables' collection.
    """
    var["Primary Energy|Ammonium|Import"] = _extract(SUPPLY, carrier="import NH3")

    return var


def primary_wind(var: dict) -> dict:
    """
    Calculate wind energy generated per region.

    Parameters
    ----------
    n
        A solved network.
    var
        The collection to append new variables in.

    Returns
    -------
    :
        The updated variables' collection.
    """
    var["Primary Energy|Wind|Onshore"] = _extract(
        SUPPLY, carrier="onwind", component="Generator"
    )
    var["Primary Energy|Wind|Offshore"] = _extract(
        SUPPLY, carrier=["offwind-ac", "offwind-dc"], component="Generator"
    )
    var["Primary Energy|Wind"] = _sum_by_subcategory(var, "Wind")

    return var


def primary_heat(var: dict) -> dict:
    """
    Calculate heat generation and enthalpy of evaporation.

    Parameters
    ----------
    var
        The collection to append new variables in.

    Returns
    -------
    :
        The updated variables' collection.
    """
    heat_bus_carrier = ["rural heat", "urban decentral heat", "urban central heat"]
    # link_energy_balance = n.statistics.energy_balance(
    #     groupby=["location", "carrier", "bus_carrier"],
    #     comps="Link",
    # )
    enthalpy_heat = (
        collect_myopic_statistics(networks, comps="Link", statistic="energy_balance")
        .pipe(filter_for_carrier_connected_to, heat_bus_carrier)
        .drop(["co2", "co2 stored"], level=DM.BUS_CARRIER)
        .pipe(calculate_input_share, heat_bus_carrier)
        .pipe(filter_by, bus_carrier=["ambient heat", "latent heat"])
        .pipe(insert_index_level, "MWh_th", "unit")
    )

    # # for every heat bus, calculate the amounts of supply for heat
    # to_concat = []
    # for bc in heat_bus_carrier:
    #     p = (
    #         link_energy_balance.pipe(filter_for_carrier_connected_to, bc)
    #         # CO2 supply are CO2 emissions that do not help heat production
    #         .drop(["co2", "co2 stored"], level=DM.BUS_CARRIER)
    #         .pipe(calculate_input_share, bc)
    #         .pipe(filter_by, bus_carrier=["ambient heat", "latent heat"])
    #     )
    #     p.attrs["unit"] = "MWh_th"
    #     to_concat.append(p)
    #
    # heat_supply = pd.concat(to_concat)
    var["Primary Energy|Heat|Latent"] = (
        filter_by(enthalpy_heat, bus_carrier="latent heat").groupby(IDX).sum()
    )
    var["Primary Energy|Heat|Ambient"] = (
        filter_by(enthalpy_heat, bus_carrier="ambient heat").groupby(IDX).sum()
    )

    solar_thermal_carr = [
        c for c in SUPPLY.index.unique("carrier") if "solar thermal" in c
    ]
    var["Primary Energy|Heat|Solar"] = _extract(
        SUPPLY, carrier=solar_thermal_carr, component="Generator"
    )
    var["Primary Energy|Heat|Geothermal"] = _extract(SUPPLY, carrier="geothermal heat")
    var["Primary Energy|Heat"] = _sum_by_subcategory(var, "Heat")

    return var


def get_renewable_generation(var: dict) -> dict:
    """

    Parameters
    ----------
    var

    Returns
    -------
    :
    """
    # solar
    # wind
    # biomass
    # heat?
    #


def combine_variables(var: dict) -> pd.Series:
    """
    Combine variables into a single dataframe.

    Parameters
    ----------
    var

    Returns
    -------
    :
    """
    ds = pd.concat({k: v for k, v in var.items() if not v.empty})
    ds.index = ds.index.rename({None: "Variable"})

    return ds


def collect_system_cost() -> pd.Series:
    """
    Extract total energy system cost per region.

    Parameters
    ----------
    n
        The pypsa network instance.

    Returns
    -------
    :
        Total CAPEX plus OPEX per model region in billions EUR (2020).
    """
    # Nodal OPEX and nodal CAPEX in billion EUR2020
    var = {}
    # CAPEX and OPEX units are incorrect and need to be updated
    unit = "billion EUR2020"
    # CAPEX and OPEX are not used anywhere else, hence they are local
    myopic_opex = (
        collect_myopic_statistics(networks, "opex", **kwargs)
        .pipe(rename_aggregate, unit, level="unit")
        .div(1e9)
    )
    myopic_capex = (
        collect_myopic_statistics(networks, "capex", **kwargs)
        .pipe(rename_aggregate, unit, level="unit")
        .div(1e9)
    )

    # FixMe: Why are there negative values in OPEX?

    # rewrite using global myopic metrics
    var["System Costs|OPEX"] = myopic_opex.groupby(IDX).sum()
    var["System Costs|CAPEX"] = myopic_capex.groupby(IDX).sum()
    var["System Costs"] = var["System Costs|CAPEX"] + var["System Costs|OPEX"]

    # todo: enable, or test later
    # assert vars["System Costs"].mul(1e9).sum() == n.objective, "Total system costs do not match the optimization result."
    # ds = pd.concat(var)
    # df.index = ds.index = ds.index.rename({None: "Variable"})

    return combine_variables(var)


def collect_primary_energy() -> pd.Series:
    """
    Extract all primary energy variables from the networks.

    In general, primary energy is the supply side of
    network components. If a component has both, supply and demand,
    the balance is returned and a warning is raised.

    Variables for primary energy follow the naming scheme:
    Primary Energy|<Category>|<SubCategory>

    Returns
    -------
    :
        The primary energy for all regions and years.
    """
    var = {}
    primary_gas(var)
    primary_oil(var)
    primary_hydrogen(var)
    primary_waste(var)
    primary_coal(var)
    primary_biomass(var)
    primary_hydro(var)
    primary_solar(var)
    primary_nuclear(var)
    primary_ammonia(var)
    primary_wind(var)
    primary_heat(var)

    return combine_variables(var)


def collect_secondary_energy(n) -> pd.DataFrame:
    """Extract all secondary energy variables from the networks."""


def collect_final_energy(n) -> pd.DataFrame:
    """Extract all final energy variables from the networks."""


if __name__ == "__main__":
    if "snakemake" not in globals():
        snakemake = mock_snakemake(
            "export_iamc_variables",
            run="KN2045_Mix",
        )
    configure_logging(snakemake)

    # template = pd.read_excel(snakemake.input.template, sheet_name="data")

    # want to use 2050 first, to catch all technologies
    # during development and debugging sessions
    networks = read_networks(sorted(snakemake.input.networks, reverse=True))

    groupby = ["location", "carrier", "bus_carrier", "unit"]
    kwargs = {
        "aggregate_components": False,
        "drop_zeros": False,
        "drop_unit": False,
    }
    myopic_energy_balance = collect_myopic_statistics(
        networks, "energy_balance", **kwargs
    )
    # calculate all statistics and process them to IAMC data model.
    SUPPLY = collect_myopic_statistics(networks, "supply", groupby=groupby, **kwargs)
    DEMAND = collect_myopic_statistics(
        networks, "withdrawal", groupby=groupby, **kwargs
    )
    IMPORT_FOREIGN = collect_myopic_statistics(
        networks, "trade_energy", scope=TradeTypes.FOREIGN, direction="import", **kwargs
    )
    EXPORT_FOREIGN = collect_myopic_statistics(
        networks, "trade_energy", scope=TradeTypes.FOREIGN, direction="export", **kwargs
    )
    IMPORT_DOMESTIC = collect_myopic_statistics(
        networks,
        "trade_energy",
        scope=TradeTypes.DOMESTIC,
        direction="import",
        **kwargs,
    )
    EXPORT_DOMESTIC = collect_myopic_statistics(
        networks,
        "trade_energy",
        scope=TradeTypes.DOMESTIC,
        direction="export",
        **kwargs,
    )
    # Idea: calculate all once and extract from global mutable series
    # to ensure nothing is forgotten. The global series however is a risk.

    # iamc_variables = []
    # for year, n in networks.items():
    system_cost = collect_system_cost()
    primary_energy = collect_primary_energy()
    # iamc_variables.append(collect_system_cost())
    # iamc_variables.append(collect_primary_energy(n))
    # iamc_variables.append(collect_secondary_energy(n))
    # iamc_variables.append(collect_final_energy(n))

    # df = pd.concat(iamc_variables)
    df = pd.concat([system_cost, primary_energy])

    df = insert_index_level(df, "PyPSA-AT", "model")
    df = insert_index_level(df, snakemake.wildcards.run, "scenario")
    df.index = df.index.rename({"location": "region"})

    # # drop all values for dummy template

    iamc = IamDataFrame(df)

    meta = pd.Series(
        {
            "Model": "PyPSA-AT",
            "Commit SHA": "",
            "Repository": "https://gitlab.aggm.at/philip.worschischek/pypsa-at",
            "Scenario": snakemake.wildcards.run,
            "Quality Assessment": "draft",
            "Release for publication": "no",
        }
    ).to_frame("value")
    meta.index.name = "key"
    with pd.ExcelWriter(snakemake.output.exported_variables) as writer:
        iamc.to_excel(writer, sheet_name="data", index=False)
        meta.to_excel(writer, sheet_name="meta", index=False)
