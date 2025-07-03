import pandas as pd
from pypsa.statistics import port_efficiency

from evals.constants import TradeTypes
from evals.utils import filter_by


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
        pd.concat([var[v] for v in var.keys() if f"|{subcat}|" in v])
        .groupby("location")
        .sum()
    )


def primary_oil(n, var):
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
    n
        A solved network.
    var
        The collection to append new variables in.

    Returns
    -------
    :
        The updated variables' collection.
    """
    netted_liquids = n.statistics.energy_balance(
        groupby="location", bus_carrier="oil", comps="Link"
    ).drop("EU", errors="ignore")
    liquids_consumption = netted_liquids[netted_liquids.le(0)].abs()

    # increase fossil oil demand by this factor to account for refining losses
    oil_refining_efficiency = (
        port_efficiency(n, "Link", "1").filter(like="oil refining").item()
    )
    # calculate the split of fossil oil generation to liquids production and assume
    # the same split for all regions
    total_oil_import = n.statistics.supply(bus_carrier="oil primary").item()
    total_liquids_production = n.statistics.supply(
        groupby="bus_carrier", bus_carrier="oil", comps="Link"
    ).item()
    fossil_share = total_oil_import / total_liquids_production
    oil = liquids_consumption * fossil_share / oil_refining_efficiency

    liquids = liquids_consumption * (1 - fossil_share)

    var["Primary Energy|Oil|Fossil"] = oil  # including losses
    var["Primary Energy|Oil|Liquids"] = liquids
    var["Primary Energy|Oil"] = oil + liquids
    var["Primary Energy|Oil|Refining Losses"] = oil * (1 - oil_refining_efficiency)

    return var


def primary_gas(n, var):
    # for scope in (TradeTypes.DOMESTIC, TradeTypes.FOREIGN):
    #     gas_trade = n.statistics.trade_energy(
    #         bus_carrier="gas", direction="import", scope=scope
    #     )
    #     gas_trade = gas_trade[gas_trade.gt(0)].groupby("location").sum()
    #     var[f"Primary Energy|Gas|Import {scope.title()}"] = gas_trade
    _get_traded_energy(n, var, "gas", "import", "Gas")

    gas_generation = n.statistics.supply(
        groupby=["location", "carrier"], bus_carrier="gas", comps="Generator"
    )
    # var["Primary Energy|Gas|Import Foreign"] = gas_trade_foreign
    # var["Primary Energy|Gas|Import Domestic"] = gas_trade_domestic
    var["Primary Energy|Gas|Production"] = (
        filter_by(gas_generation, carrier="gas").groupby("location").sum()
    )
    var["Primary Energy|Gas|Import Global"] = (
        filter_by(gas_generation, carrier="import gas").groupby("location").sum()
    )
    # summing up does not work bc of misaligned indices. Need to concat.
    var["Primary Energy|Gas"] = _sum_by_subcategory(var, "Gas")

    return var


def primary_waste(n, var):
    """

    Parameters
    ----------
    n
    var

    Returns
    -------
    :
    """
    _get_traded_energy(n, var, "municipal solid waste", "import", "Waste")

    var["Primary Energy|Waste|Solid"] = n.statistics.supply(
        groupby="location", bus_carrier="municipal solid waste", comps="Generator"
    )
    var["Primary Energy|Waste"] = _sum_by_subcategory(var, "Waste")

    return var


def primary_coal(n, var):
    """
    Calculate the amounts of coal consumed in a region.

    Coal is not produced by any Link, therefore it's safe to assume
    all withdrawal is imported fossil coal or lignite.

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

    var["Primary Energy|Coal|Hard"] = n.statistics.withdrawal(
        groupby="location", bus_carrier="coal"
    ).drop("Store")
    var["Primary Energy|Coal|Lignite"] = n.statistics.withdrawal(
        groupby="location", bus_carrier="lignite"
    ).drop("Store")
    var["Primary Energy|Coal"] = (
        var["Primary Energy|Coal|Hard"] + var["Primary Energy|Coal|Lignite"]
    )

    return var


def primary_hydrogen(n, var):
    """
    Calculate the amounts of hydrogen imported into a region.

    There are global import of Hydrogen, found in `Generator`
    components, and various types of H2 pipelines, that bring
    Hydrogen into regions.

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

    _get_traded_energy(n, var, "H2", "import", "Hydrogen")

    var["Primary Energy|Hydrogen|Import Global"] = (
        n.statistics.supply(
            groupby=["location", "carrier"], bus_carrier="H2", comps="Generator"
        )
        .pipe(filter_by, carrier="import H2")
        .droplevel("carrier")
    )

    var["Primary Energy|Hydrogen"] = _sum_by_subcategory(var, "Hydrogen")

    return var


def primary_biomass(n, var):
    """
    Calculate the amounts of biomass generated in a region.

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

    _get_traded_energy(n, var, "solid biomass", "import", "Biomass")

    var["Primary Energy|Biomass|Solid"] = n.statistics.supply(
        groupby="location", bus_carrier="solid biomass", comps="Generator"
    )
    var["Primary Energy|Biomass|Biogas"] = n.statistics.supply(
        groupby="location", bus_carrier="biogas", comps="Generator"
    )
    var["Primary Energy|Biomass"] = _sum_by_subcategory(var, "Biomass")

    return var


def primary_hydro(n, var):
    """

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

    n.statistics.supply(
        bus_carrier="AC", groupby=["location", "carrier"], comps="StorageUnit"
    ).pipe(filter_by, location="AT32")

    hydro = n.statistics.phs_split(drop_hydro_cols=False)

    var["Primary Energy|Hydro|PHS"] = filter_by(
        hydro, carrier="PHS Dispatched Power from Inflow"
    ).droplevel(["carrier", "bus_carrier"])
    var["Primary Energy|Hydro|Reservoir"] = filter_by(
        hydro, carrier="hydro Dispatched Power from Inflow"
    ).droplevel(["carrier", "bus_carrier"])
    var["Primary Energy|Hydro|Run-of-River"] = (
        n.statistics.supply(
            groupby=["location", "carrier"], comps="Generator", bus_carrier="AC"
        )
        .pipe(filter_by, carrier="ror")
        .droplevel("carrier")
    )
    var["Primary Energy|Hydro"] = _sum_by_subcategory(var, "Hydro")

    return var


def primary_solar(n, var):
    """

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
    solar = n.statistics.supply(
        groupby=["location", "carrier"],
        comps="Generator",
        bus_carrier=["AC", "low voltage"],
    )
    var["Primary Energy|Solar|solar"] = filter_by(solar, carrier="solar").droplevel(
        "carrier"
    )
    var["Primary Energy|Solar|solar-hsat"] = filter_by(
        solar, carrier="solar-hsat"
    ).droplevel("carrier")
    var["Primary Energy|Solar|solar rooftop"] = filter_by(
        solar, carrier="solar rooftop"
    ).droplevel("carrier")
    var["Primary Energy|Solar"] = _sum_by_subcategory(var, "Solar")

    return var


def primary_nuclear(n, var):
    """

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

    return var


def primary_ammonia(n, var):
    """

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
    return var


def primary_wind(n, var):
    """

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
    return var


def primary_heat(n, var):
    """

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

    # ambient heat from heat pumps
    # latent heat from CHPs
    # Solar heat
    # Geothermal

    return var


def get_renewable_generation(var):
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
