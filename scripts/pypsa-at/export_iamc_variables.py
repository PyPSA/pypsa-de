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

from evals.constants import BusCarrier, TradeTypes
from evals.constants import DataModel as DM
from evals.fileio import read_networks
from evals.statistic import collect_myopic_statistics
from evals.utils import (
    calculate_input_share,
    filter_by,
    filter_for_carrier_connected_to,
    get_transmission_techs,
    insert_index_level,
    rename_aggregate,
)
from scripts._helpers import configure_logging, mock_snakemake

logger = logging.getLogger(__file__)

IDX = [DM.YEAR, DM.LOCATION, "unit"]


def _extract(ds: pd.Series, **filter_kwargs) -> pd.Series:
    """Extract and group filter results."""
    results = filter_by(ds, **filter_kwargs)
    ds.drop(results.index, inplace=True)
    return results.groupby(IDX).sum()


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


def _sum_variables_by_prefix(var, prefix):
    return (
        pd.concat([var[v] for v in var.keys() if v.startswith(prefix)])
        .groupby(IDX)
        .sum()
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

    # assuming that all local oil production is consumed locally.
    # Let's not filter_by components, to capture all but Stores.
    production = (
        filter_by(SUPPLY, bus_carrier=["oil"])  # , "unsustainable bioliquids"
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
    regional_oil_deficit = consumption.add(production, fill_value=0).clip(lower=0)

    # assert primary oil and imports are equal to regional demands
    eu_oil_import = filter_by(SUPPLY, carrier=["oil primary", "import oil"])
    try:
        pd.testing.assert_series_equal(
            regional_oil_deficit.groupby("year").sum(),
            eu_oil_import.groupby("year").sum(),
            check_names=False,
        )
    except AssertionError:
        logger.warning("Oil amounts mismatch.")

    var["Primary Energy|Oil|Fossil Oil"] = (
        regional_oil_deficit / oil_refining_eff
    )  # total amounts before refining
    var["Primary Energy|Oil|Refining Losses"] = regional_oil_deficit * (
        1 - oil_refining_eff
    )
    var["Primary Energy|Oil|Global Import"] = _extract(
        SUPPLY, carrier="import oil"
    )  # EU
    var["Primary Energy|Oil|Primary"] = _extract(SUPPLY, carrier="oil primary")  # EU
    var["Primary Energy|Oil|Refining"] = _extract(SUPPLY, carrier="oil refining")  # EU

    # unsustainable bioliquids have regional bus generators
    var["Primary Energy|Oil|Unsustainable Bioliquids"] = _extract(
        SUPPLY, carrier="unsustainable bioliquids", component="Generator"
    )

    var["Primary Energy|Oil"] = var["Primary Energy|Oil|Fossil Oil"].add(
        var["Primary Energy|Oil|Unsustainable Bioliquids"], fill_value=0
    )
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
    var["Primary Energy|Gas"] = _sum_variables_by_prefix(var, "Primary Energy|Gas")

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
    var["Primary Energy|Waste"] = _sum_variables_by_prefix(var, "Primary Energy|Waste")

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
    ).mul(-1)
    var["Primary Energy|Coal|Lignite"] = (
        filter_by(DEMAND, bus_carrier="lignite", component="Link").groupby(IDX).sum()
    ).mul(-1)
    var["Primary Energy|Coal"] = _sum_variables_by_prefix(var, "Primary Energy|Coal")

    # remove EU coal generators from the to-do list
    coal_generators = filter_by(
        SUPPLY, bus_carrier=["coal", "lignite"], component="Generator"
    )
    SUPPLY.drop(coal_generators.index, inplace=True)

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
    var["Primary Energy|H2|Import Foreign"] = _extract(
        IMPORT_FOREIGN, bus_carrier=bus_carrier
    )
    var["Primary Energy|H2|Import Domestic"] = _extract(
        IMPORT_DOMESTIC, bus_carrier=bus_carrier
    )
    var["Primary Energy|H2|Import Global"] = _extract(
        SUPPLY, carrier="import H2", bus_carrier=bus_carrier, component="Generator"
    )

    var["Primary Energy|H2"] = _sum_variables_by_prefix(var, "Primary Energy|H2")

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
    var["Primary Energy|Biomass"] = _sum_variables_by_prefix(
        var, "Primary Energy|Biomass"
    )

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
    var["Primary Energy|Hydro"] = _sum_variables_by_prefix(var, "Primary Energy|Hydro")

    return var


def primary_solar(var: dict) -> dict:
    """
    Calculate solar energy generated per region.

    Parameters
    ----------
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
    var["Primary Energy|Solar"] = _sum_variables_by_prefix(var, "Primary Energy|Solar")

    return var


def primary_liquids(var: dict) -> dict:
    """
    Calculate the amounts of liquid fuels generated per region.

    Parameters
    ----------
        The collection to append new variables in.

    Returns
    -------
    :
        The updated variables' collection.
    """
    var["Primary Energy|Liquids|Unsustainable Bioliquids"] = _extract(
        SUPPLY, carrier="unsustainable bioliquids", component="Generator"
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
        filter_by(DEMAND, carrier="uranium", component="Link")
        .groupby(IDX)
        .sum()
        .mul(-1)
    )
    uranium_generators = filter_by(SUPPLY, bus_carrier="uranium", component="Generator")
    SUPPLY.drop(uranium_generators.index, inplace=True)

    return var


def primary_ammonia(var: dict) -> dict:
    """
    Calculate the ammonium imported per region.

    Parameters
    ----------
    var
        The collection to append new variables in.

    Returns
    -------
    :
        The updated variables' collection.
    """
    # there is no regional ammonium demand
    var["Primary Energy|Ammonium|Import"] = _extract(SUPPLY, carrier="import NH3")

    return var


def primary_methanol(var: dict) -> dict:
    """
    Calculate methanol imported per region.

    Parameters
    ----------
    var
        The collection to append new variables in.

    Returns
    -------
    :
        The updated variables' collection.
    """
    # todo: same logic as oil for regional demands.
    # regional_demand = filter_by(DEMAND, bus_carrier="methanol", component="Link").
    # regional_supply = filter_by(SUPPLY, bus_carrier="methanol", component="Link").sum()
    # regional_supply + regional_demand
    # filter_by(SUPPLY, carrier="import methanol", component="Link").sum()
    var["Primary Energy|Methanol|Import"] = _extract(SUPPLY, carrier="import methanol")

    return var


def primary_wind(var: dict) -> dict:
    """
    Calculate wind energy generated per region.

    Parameters
    ----------
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
    var["Primary Energy|Wind"] = _sum_variables_by_prefix(var, "Primary Energy|Wind")

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
        LINK_BALANCE.pipe(filter_for_carrier_connected_to, heat_bus_carrier)
        .drop(["co2", "co2 stored"], level=DM.BUS_CARRIER)
        .pipe(calculate_input_share, heat_bus_carrier)
        .pipe(filter_by, bus_carrier=["ambient heat", "latent heat"])
        .pipe(insert_index_level, "MWh_th", "unit")
    )

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
    var["Primary Energy|Heat"] = _sum_variables_by_prefix(var, "Primary Energy|Heat")

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
    primary_liquids(var)
    primary_methanol(var)

    assert filter_by(SUPPLY, component="Generator").empty

    return combine_variables(var)


def secondary_electricity_supply(var: dict) -> dict:
    """

    Parameters
    ----------
    var

    Returns
    -------
    :
    """
    prefix = "Secondary Energy|Electricity"
    bc = ["AC", "low voltage"]

    var[f"{prefix}|Gas"] = _extract(
        SUPPLY,
        carrier=[
            "CCGT",
            "OCGT",
            "urban central gas CHP",
        ],
        bus_carrier=bc,
    )
    var[f"{prefix}|Oil"] = _extract(
        SUPPLY, carrier="urban central oil CHP", bus_carrier=bc
    )
    var[f"{prefix}|Coal|Hard Coal"] = _extract(
        SUPPLY, carrier=["coal", "urban central coal CHP"], bus_carrier=bc
    )
    var[f"{prefix}|Coal|Lignite"] = _extract(
        SUPPLY, carrier=["lignite", "urban central lignite CHP"], bus_carrier=bc
    )
    var[f"{prefix}|H2"] = _extract(
        SUPPLY,
        carrier=["urban central H2 CHP", "urban central H2 retrofit CHP"],
        bus_carrier=bc,
    )
    var[f"{prefix}|Solid Biomass"] = _extract(
        SUPPLY,
        carrier=["solid biomass", "urban central solid biomass CHP"],
        bus_carrier=bc,
    )
    var[f"{prefix}|Nuclear"] = _extract(SUPPLY, carrier="nuclear", bus_carrier=bc)
    var[f"{prefix}|Waste|w/o CC"] = _extract(
        SUPPLY, carrier="waste CHP", bus_carrier=bc
    )
    var[f"{prefix}|Waste|w CC"] = _extract(
        SUPPLY, carrier="waste CHP CC", bus_carrier=bc
    )

    var[prefix] = _sum_variables_by_prefix(var, prefix)

    # subcategory aggregation must happen after prefix aggregations, or the
    # sum will be distorted
    var[f"{prefix}|Coal"] = _sum_variables_by_prefix(var, f"{prefix}|Coal")
    var[f"{prefix}|Waste"] = _sum_variables_by_prefix(var, f"{prefix}|Waste")

    # distribution grid losses are no supply, but we deal with it now remove all
    # electricity from the global supply statistic
    var[f"{prefix}|Distribution Grid Losses"] = _extract(
        SUPPLY, carrier="electricity distribution grid"
    ) + _extract(DEMAND, carrier="electricity distribution grid")  # negative values

    assert filter_by(SUPPLY, bus_carrier=bc, component="Link").empty

    return var


def secondary_gas_supply(var: dict) -> dict:
    """

    Parameters
    ----------
    var

    Returns
    -------
    :
    """
    # Secondary Energy|<output bus_carrier>|from <input bus_carrier>|<subcategory>
    prefix = "Secondary Energy|Gas"
    bc = "gas"
    var[f"{prefix}|Biogas|w/o CC"] = _extract(
        SUPPLY, carrier="biogas to gas", bus_carrier=bc
    )
    var[f"{prefix}|Biogas|w CC"] = _extract(
        SUPPLY, carrier="biogas to gas CC", bus_carrier=bc
    )
    # _extract(DEMAND, carrier=["biogas to gas", "biogas to gas CC"], bus_carrier="biogas")

    var[f"{prefix}|Solid Biomass|w/o CC"] = _extract(
        SUPPLY, carrier="BioSNG", bus_carrier=bc
    )
    var[f"{prefix}|Solid Biomass|w CC"] = _extract(
        SUPPLY, carrier="BioSNG CC", bus_carrier=bc
    )
    # _extract(DEMAND, carrier=["BioSNG", "BioSNG CC"], bus_carrier="solid biomass")

    var[f"{prefix}|Sabatier"] = _extract(SUPPLY, carrier="Sabatier", bus_carrier=bc)

    var[prefix] = _sum_variables_by_prefix(var, prefix)

    assert filter_by(SUPPLY, bus_carrier=bc, component="Link").empty

    # move to Carbon function
    # var["Secondary Carbon|Biogas|atmosphere"] = _extract(DEMAND, carrier="biogas to gas", bus_carrier="co2")  # CO2 credit
    # var["Secondary Carbon|Biogas|stored"] = _extract(SUPPLY, carrier="biogas to gas CC", bus_carrier="co2")  # carbon capture

    return var


def secondary_hydrogen_supply(var: dict) -> dict:
    """

    Parameters
    ----------
    var

    Returns
    -------
    :
    """
    prefix = "Secondary Energy|H2"
    bc = "H2"
    var[f"{prefix}|Electricity"] = _extract(
        SUPPLY, carrier="H2 Electrolysis", bus_carrier=bc
    )
    # var[f"{prefix}|Electricity|Losses"] = _extract(SUPPLY, carrier="H2 Electrolysis", bus_carrier=bc)
    var[f"{prefix}|Gas|SMR w/o CC"] = _extract(SUPPLY, carrier="SMR", bus_carrier=bc)
    var[f"{prefix}|Gas|SMR w CC"] = _extract(SUPPLY, carrier="SMR CC", bus_carrier=bc)
    var[f"{prefix}|Methanol|w/o CC"] = _extract(
        SUPPLY, carrier="Methanol steam reforming", bus_carrier=bc
    )
    var[f"{prefix}|Methanol|w CC"] = _extract(
        SUPPLY, carrier="Methanol steam reforming CC", bus_carrier=bc
    )

    var[prefix] = _sum_variables_by_prefix(var, prefix)

    assert filter_by(SUPPLY, bus_carrier=bc, component="Link").empty

    return var


def secondary_methanol_supply(var: dict) -> dict:
    """

    Parameters
    ----------
    var

    Returns
    -------
    :
    """
    prefix = "Secondary Energy|Methanol"
    bc = "methanol"

    # need to distinguish between methanol and heat output
    methanolisation_inputs_for_methanol = (
        LINK_BALANCE.drop(["co2", "co2 stored"], level="bus_carrier")
        .pipe(filter_for_carrier_connected_to, bc)
        .pipe(calculate_input_share, bc)
        .pipe(filter_by, carrier="methanolisation")
    )
    var[f"{prefix}|H2"] = (
        filter_by(methanolisation_inputs_for_methanol, bus_carrier="H2")
        .pipe(insert_index_level, "MWh_LHV", "unit")
        .groupby(IDX)
        .sum()
    )
    var[f"{prefix}|AC"] = (
        filter_by(methanolisation_inputs_for_methanol, bus_carrier="AC")
        .pipe(insert_index_level, "MWh_el", "unit")
        .groupby(IDX)
        .sum()
    )
    _extract(SUPPLY, carrier="methanolisation", bus_carrier=bc)

    var[prefix] = _sum_variables_by_prefix(var, prefix)

    assert filter_by(SUPPLY, bus_carrier=bc, component="Link").empty
    # biomass-to-methanol?

    return var


def secondary_oil(var: dict) -> dict:
    """

    Parameters
    ----------
    var

    Returns
    -------
    :
    """
    prefix = "Secondary Energy|Oil"
    bc = "oil"

    var[f"{prefix}|Solid Biomass|Biomass2Liquids|w/o CC"] = _extract(
        SUPPLY, carrier="biomass to liquid", bus_carrier=bc
    )
    var[f"{prefix}|Solid Biomass|Biomass2Liquids|w CC"] = _extract(
        SUPPLY, carrier="biomass to liquid CC", bus_carrier=bc
    )

    # electrobiofuels has 2 inputs: solid biomass and H2 and one output
    electrobiofuels_inputs_for_oil = (
        LINK_BALANCE.drop(["co2", "co2 stored"], level="bus_carrier")
        .pipe(filter_for_carrier_connected_to, bc)
        .pipe(calculate_input_share, bc)
        .pipe(filter_by, carrier="electrobiofuels")
    )
    var[f"{prefix}|H2|Electrobiofuels"] = (
        filter_by(electrobiofuels_inputs_for_oil, bus_carrier="H2")
        .pipe(insert_index_level, "MWh_LHV", "unit")
        .groupby(IDX)
        .sum()
    )
    var[f"{prefix}|AC|Electrobiofuels"] = (
        filter_by(electrobiofuels_inputs_for_oil, bus_carrier="solid biomass")
        .pipe(insert_index_level, "MWh_LHV", "unit")
        .groupby(IDX)
        .sum()
    )
    _extract(SUPPLY, carrier="electrobiofuels", bus_carrier=bc)

    var[f"{prefix}|H2|Fischer-Tropsch"] = _extract(
        SUPPLY, carrier="Fischer-Tropsch", bus_carrier=bc
    )

    var[f"{prefix}|Unsustainable Bioliquids"] = _extract(
        SUPPLY, carrier="unsustainable bioliquids", bus_carrier=bc
    )

    assert filter_by(SUPPLY, bus_carrier=bc, component="Link").empty

    var[prefix] = _sum_variables_by_prefix(var, prefix)
    var[f"{prefix}|H2"] = _sum_variables_by_prefix(var, f"{prefix}|H2")
    var[f"{prefix}|AC"] = _sum_variables_by_prefix(var, f"{prefix}|AC")
    var[f"{prefix}|Solid Biomass"] = _sum_variables_by_prefix(
        var, f"{prefix}|Solid Biomass"
    )

    return var


def secondary_ammonia(var: dict) -> dict:
    """

    Parameters
    ----------
    var

    Returns
    -------
    :
    """
    prefix = "Secondary Energy|NH3"
    bc = "NH3"

    haber_bosch_input_for_nh3 = (
        LINK_BALANCE.drop(["co2", "co2 stored"], level="bus_carrier")
        .pipe(filter_for_carrier_connected_to, bc)
        .pipe(calculate_input_share, bc)
        .pipe(filter_by, carrier="Haber-Bosch")
        .pipe(insert_index_level, "n/a", "unit")
    )
    var[f"{prefix}|AC"] = (
        filter_by(haber_bosch_input_for_nh3, carrier="Haber-Bosch", bus_carrier="AC")
        .groupby(IDX)
        .sum()
    )
    var[f"{prefix}|H2"] = (
        filter_by(haber_bosch_input_for_nh3, carrier="Haber-Bosch", bus_carrier="H2")
        .groupby(IDX)
        .sum()
    )
    _extract(SUPPLY, carrier="Haber-Bosch", bus_carrier=bc)

    return var


def secondary_heat(var: dict) -> dict:
    """

    Parameters
    ----------
    var

    Returns
    -------
    :
    """
    prefix = "Secondary Energy|Heat"
    bc = BusCarrier.heat_buses()

    var[f"{prefix}|Solid Biomass|Boiler"] = _extract(
        SUPPLY,
        carrier=["rural biomass boiler", "urban decentral biomass boiler"],
        bus_carrier=bc,
    )
    var[f"{prefix}|Solid Biomass|CHP"] = _extract(
        SUPPLY, carrier="urban central solid biomass CHP", bus_carrier=bc
    )

    var[f"{prefix}|AC|Ground Heat Pump"] = _extract(
        SUPPLY, carrier="rural ground heat pump", bus_carrier=bc
    )
    var[f"{prefix}|AC|Air Heat Pump"] = _extract(
        SUPPLY,
        carrier=[
            "urban decentral air heat pump",
            "rural air heat pump",
            "urban central air heat pump",
        ],
        bus_carrier=bc,
    )
    var[f"{prefix}|Gas|Boiler"] = _extract(
        SUPPLY,
        carrier=[
            "rural gas boiler",
            "urban central gas boiler",
            "urban decentral gas boiler",
        ],
        bus_carrier=bc,
    )

    var[f"{prefix}|AC|Resistive Heater"] = _extract(
        SUPPLY,
        carrier=[
            "rural resistive heater",
            "urban decentral resistive heater",
            "urban central resistive heater",
        ],
        bus_carrier=bc,
    )

    var[f"{prefix}|Oil|Boiler"] = _extract(
        SUPPLY,
        carrier=["rural oil boiler", "urban decentral oil boiler"],
        bus_carrier=bc,
    )

    var[f"{prefix}|H2|CHP"] = _extract(
        SUPPLY,
        carrier=["urban central H2 CHP", "urban central H2 retrofit CHP"],
        bus_carrier=bc,
    )

    var[f"{prefix}|Gas|CHP"] = _extract(
        SUPPLY, carrier="urban central gas CHP", bus_carrier=bc
    )

    var[f"{prefix}|Oil|CHP"] = _extract(
        SUPPLY,
        carrier="urban central oil CHP",
        bus_carrier=bc,
    )

    var[f"{prefix}|Coal|CHP"] = _extract(
        SUPPLY,
        carrier=["urban central coal CHP", "urban central lignite CHP"],
        bus_carrier=bc,
    )

    var[f"{prefix}|Waste|CHP w/o CC"] = _extract(
        SUPPLY, carrier="waste CHP", bus_carrier=bc
    )
    var[f"{prefix}|Waste|CHP w CC"] = _extract(
        SUPPLY, carrier="waste CHP CC", bus_carrier=bc
    )

    # Excess heat from technologies
    var[f"{prefix}|H2|Sabatier"] = _extract(SUPPLY, carrier="Sabatier", bus_carrier=bc)
    var[f"{prefix}|H2|Fischer-Tropsch"] = _extract(
        SUPPLY, carrier="Fischer-Tropsch", bus_carrier=bc
    )
    var[f"{prefix}|AC|Electrolysis"] = _extract(
        SUPPLY, carrier="H2 Electrolysis", bus_carrier=bc
    )

    # methanolisation and Haber-Bosch have multiple inputs
    methanolisation_input_for_heat = (
        LINK_BALANCE.drop(["co2", "co2 stored"], level="bus_carrier")
        .pipe(filter_for_carrier_connected_to, bc)
        .pipe(calculate_input_share, bc)
        .pipe(filter_by, carrier=["methanolisation", "Haber-Bosch"])
        .pipe(insert_index_level, "MWh_LHV", "unit")
    )
    var[f"{prefix}|AC|Methanolisation"] = (
        filter_by(
            methanolisation_input_for_heat, carrier="methanolisation", bus_carrier="AC"
        )
        .groupby(IDX)
        .sum()
    )
    var[f"{prefix}|H2|Methanolisation"] = (
        filter_by(
            methanolisation_input_for_heat, carrier="methanolisation", bus_carrier="H2"
        )
        .groupby(IDX)
        .sum()
    )
    _extract(SUPPLY, carrier="methanolisation", bus_carrier=bc)

    var[f"{prefix}|AC|Haber-Bosch"] = (
        filter_by(
            methanolisation_input_for_heat, carrier="Haber-Bosch", bus_carrier="AC"
        )
        .groupby(IDX)
        .sum()
    )
    var[f"{prefix}|H2|Haber-Bosch"] = (
        filter_by(
            methanolisation_input_for_heat, carrier="Haber-Bosch", bus_carrier="H2"
        )
        .groupby(IDX)
        .sum()
    )
    _extract(SUPPLY, carrier="Haber-Bosch", bus_carrier=bc)

    # deal with water pit charger losses now to clear all heat buses
    var[f"{prefix}|Heat|Storage Losses"] = _extract(
        SUPPLY,
        carrier="urban central water pits discharger",
        bus_carrier="urban central heat",
    ) + _extract(
        DEMAND,
        carrier="urban central water pits charger",
        bus_carrier="urban central heat",
    )
    # drop the supply/demand at the other bus side of (dis)charger links
    _extract(
        SUPPLY,
        carrier="urban central water pits charger",
        bus_carrier="urban central water pits",
    )
    _extract(
        DEMAND,
        carrier="urban central water pits discharger",
        bus_carrier="urban central water pits",
    )

    assert filter_by(SUPPLY, bus_carrier=bc, component="Link").empty

    var[prefix] = _sum_variables_by_prefix(var, prefix)

    return var


def secondary_waste(var: dict) -> dict:
    prefix = "Secondary Energy|Waste"
    bc = "non-sequestered HVC"

    var[f"{prefix}|Waste|Naptha"] = _extract(
        SUPPLY, carrier="naphtha for industry", bus_carrier=bc
    )
    var[f"{prefix}|Waste|Municipal solid waste"] = _extract(
        SUPPLY, carrier="municipal solid waste", bus_carrier=bc
    )

    assert filter_by(SUPPLY, bus_carrier=bc, component="Link").empty

    var[prefix] = _sum_variables_by_prefix(var, prefix)

    return var


def secondary_solid_biomass_supply(var: dict) -> dict:
    var["Secondary Energy|Solid Biomass|Boiler Error"] = _extract(
        SUPPLY,
        carrier=["rural biomass boiler", "urban decentral biomass boiler"],
        bus_carrier="solid biomass",
    )

    if not var["Secondary Energy|Solid Biomass|Boiler Error"].empty:
        logger.warning(
            f"Solid biomass boilers supply to solid biomass bus. Total amount of energy supplied = {var['Secondary Energy|Solid Biomass|Boiler Error'].sum():.2f} MWh_LHV"
        )

    assert filter_by(SUPPLY, bus_carrier="solid biomass", component="Link").empty

    return var


def collect_secondary_energy() -> pd.Series:
    """Extract all secondary energy variables from the networks."""
    var = {}
    secondary_electricity_supply(var)
    secondary_gas_supply(var)
    secondary_hydrogen_supply(var)
    secondary_methanol_supply(var)
    secondary_oil(var)
    secondary_heat(var)
    secondary_ammonia(var)
    secondary_waste(var)

    # solid biomass is produced by some boilers, which is wrong of course
    # but needs to be addressed nevertheless to correct balances
    secondary_solid_biomass_supply(var)

    # Links that connect to buses with single loads. They are skipped in
    # IAMC variables, because their buses are only needed because of
    # PyPSA restrictions.
    ignore_bus_carrier = [
        "EV battery",
        "agriculture machinery oil",
        "coal for industry",
        "gas for industry",
        "industry methanol",
        "kerosene for aviation",
        "land transport oil",
        "naphtha for industry",
        "shipping methanol",
        "shipping oil",
        "solid biomass for industry",
    ]
    assert (
        filter_by(SUPPLY, component="Link")
        .drop("t_co2", level="unit", errors="ignore")
        .drop(ignore_bus_carrier, level="bus_carrier", errors="ignore")
        .empty
    )

    return combine_variables(var)


def collect_final_energy() -> pd.DataFrame:
    """Extract all final energy variables from the networks."""
    # by sector


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
    # calculate all statistics and process them to IAMC data model. The idea is to
    # calculate everything once and remove rows from the global statistic. This way
    # we make sure that nothing is counted twice or is forgotten.
    SUPPLY = collect_myopic_statistics(networks, "supply", groupby=groupby, **kwargs)
    DEMAND = collect_myopic_statistics(
        networks, "withdrawal", groupby=groupby, **kwargs
    ).mul(-1)
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

    # necessary for Links with multiple inputs
    LINK_BALANCE = collect_myopic_statistics(
        networks, comps="Link", statistic="energy_balance"
    )

    # all transmission is already in trade_energy.
    transmission_carrier = [t[1] for t in get_transmission_techs(networks)]
    SUPPLY.drop(transmission_carrier, level="carrier", errors="ignore", inplace=True)
    DEMAND.drop(transmission_carrier, level="carrier", errors="ignore", inplace=True)

    # collect transformed energy system variables
    primary_energy = collect_primary_energy()
    secondary_energy = collect_secondary_energy()
    system_cost = collect_system_cost()

    df = pd.concat([system_cost, primary_energy, secondary_energy])

    # for global_statistic in [
    #     SUPPLY,
    #     DEMAND,
    #     IMPORT_FOREIGN,
    #     EXPORT_FOREIGN,
    #     IMPORT_DOMESTIC,
    #     EXPORT_DOMESTIC,
    # ]:
    #     assert global_statistic.empty, f"Statistic not transformed: {global_statistic}."

    df = insert_index_level(df, "PyPSA-AT", "model")
    df = insert_index_level(df, snakemake.wildcards.run, "scenario")
    df.index = df.index.rename({"location": "region"})  # comply with IAMC data model

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
