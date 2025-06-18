import pandas as pd
from pypsa.statistics import port_efficiency

from evals.constants import TradeTypes
from evals.utils import filter_by


def primary_liquids(n, var):
    # "Primary Energy|Oil"  (= sum of all oil variables)
    # "Primary Energy|Oil|Heat"
    # "Primary Energy|Oil|Electricity"
    # "Primary Energy|Oil|Liquids" (= Oil minus rest)
    # what about Fischer-Tropsch?

    # need to reduce primary energy by renewable oil production in the same region
    # renewable_oil_production = n.statistics.supply(
    #     groupby=["location", "carrier", "bus_carrier"], comps="Link", bus_carrier="oil"
    # ).drop(
    #     ["unsustainable bioliquids", "oil refining", "import oil"],
    #     axis=0,
    #     level="carrier",
    #     errors="ignore",
    # )
    # oil_usage = (
    #     n.statistics.withdrawal(
    #         groupby=["location", "carrier", "bus_carrier"], bus_carrier="oil"
    #     )
    #     .drop("Store", axis=0, level="component", errors="ignore")
    #     .droplevel("component")
    #     .mul(-1)  # withdrawal is negative
    # )
    #
    # # assuming that regional oil production is consumed in the same region
    # oil_saldo = (
    #     oil_usage.groupby("location")
    #     .sum()
    #     .add(renewable_oil_production.groupby("location").sum(), fill_value=0)
    # )
    # oil_import = oil_saldo[oil_saldo.le(0)]
    # oil_export = oil_saldo[oil_saldo.gt(0)]
    #
    # var["Primary Energy|Oil"] = oil_import.divide(oil_refining_efficiency)
    # var["Final Energy|Oil"] = oil_export

    # FixMe: broken, I think
    oil_balance = n.statistics.energy_balance(
        groupby="location", bus_carrier="oil", comps="Link"
    ).drop("EU", errors="ignore")
    liquids_all = oil_balance[oil_balance.le(0)].abs()

    # increase fossil oil demand by this factor to account for refining losses
    oil_refining_efficiency = (
        port_efficiency(n, "Link", "1").filter(like="oil refining").item()
    )
    # calculate the split of fossil oil generation to liquids production and assume
    # the same split for all regions
    oil_fossil_eu = n.statistics.supply(bus_carrier="oil primary").item()
    renewable_liquids = n.statistics.supply(
        groupby="bus_carrier", bus_carrier="oil", comps="Link"
    ).item()
    fossil_fraction = liquids_all * oil_fossil_eu / renewable_liquids
    # increase fossil oil share by refining losses
    liquids_w_losses = (
        liquids_all - fossil_fraction + fossil_fraction / oil_refining_efficiency
    )
    # assuming that regional oil production is consumed in the same region, the netted
    # energy balance becomes the primary energy (import) and final energy (export) of
    # the region. It is important to note that this is not the same as fossil oil imports,
    # but all oil imports including renewable oil production from other regions.
    # Lacking an oil network that connects regions, we assume that all oil is consumed
    # locally and only oil production surplus becomes exported as final energy.
    var["Primary Energy|Liquids"] = liquids_w_losses
    var["Primary Energy|Liquids|oil refining losses"] = liquids_w_losses - liquids_all

    # should be
    # var["Primary Energy|Oil"]
    # var["Primary Energy|Oil|Refining Losses"]

    # TEST: oil import must be same as total sum of oil primary
    # pd.testing.assert_series_equal(
    #     oil_import, oil_energy_balance[oil_energy_balance.le(0)], check_names=False
    # )  # PASS -> its the same
    return var


def primary_gas(n, var):
    for scope in (TradeTypes.DOMESTIC, TradeTypes.FOREIGN):
        gas_trade = n.statistics.trade_energy(
            bus_carrier="gas", direction="import", scope=scope
        )
        gas_trade = gas_trade[gas_trade.gt(0)].groupby("location").sum()
        var[f"Primary Energy|Gas|Import {scope.title()}"] = gas_trade

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
    var["Primary Energy|Gas"] = (
        pd.concat(
            [
                var[f"Primary Energy|Gas|Import {TradeTypes.FOREIGN.title()}"],
                var[f"Primary Energy|Gas|Import {TradeTypes.DOMESTIC.title()}"],
                gas_generation.droplevel("carrier"),
            ]
        )
        .groupby("location")
        .sum()
    )

    return var


def primary_waste(n, var):
    # non-sequestered HVC + municipal solid waste
    # municipal solid waste is generated regionally
    # municipal solid waste supplies to HVC bus + draws from CO2 atmosphere
    # naphtha for industry withdraws from oil and supplies naphtha + waste to HVC bus + process emissions
    # waste CHP only draws from HVC bus
    # Primary Energy|Waste is only municipal solid waste Generators, the rest is secondary energy
    # plus 'municipal solid waste transport' import amounts
    # n.statistics.supply(groupby=["location", "carrier", "bus_carrier"], bus_carrier=["non-sequestered HVC", "municipal solid waste"])
    waste_generation = n.statistics.supply(
        groupby="location", bus_carrier="municipal solid waste", comps="Generator"
    )
    waste_import = n.statistics.trade_energy(
        bus_carrier="municipal solid waste",
        direction="import",
        scope=("domestic", "foreign"),
    )
    waste_import = waste_import[waste_import.gt(0)].groupby("location").sum()

    var["Primary Energy|Waste"] = waste_generation.add(waste_import, fill_value=0)
    var["Primary Energy|Waste|Import"] = waste_import

    return var


def primary_coal(n, var):
    return var


def primary_hydrogen(n, vars):
    return vars


def primary_biomass(n, vars):
    return vars


def primary_hydro(n, vars):
    return vars


def primary_solar(n, vars):
    return vars


def primary_nuclear(n, vars):
    return vars


def primary_ammonia(n, vars):
    return vars


def primary_wind(n, vars):
    return vars


def primary_heat(n, vars):
    return vars
