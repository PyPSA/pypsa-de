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
import re

import pandas as pd
from pyam import IamDataFrame
from pypsa.statistics import port_efficiency

from evals.constants import DataModel as DM
from evals.constants import TradeTypes
from evals.fileio import read_networks
from evals.statistic import collect_myopic_statistics
from evals.utils import (
    filter_by,
    get_transmission_techs,
    insert_index_level,
    rename_aggregate,
)
from scripts._helpers import configure_logging, mock_snakemake

logger = logging.getLogger(__file__)

YEAR_LOC = [DM.YEAR, DM.LOCATION]
IDX = YEAR_LOC + ["unit"]
PRIMARY = "Primary Energy"
SECONDARY = "Secondary Energy"
FINAL = "Final Energy"

BC_ALIAS = {
    "urban central heat": "Heat",
    "urban decentral heat": "Heat",
    "rural heat": "Heat",
    "gas": "Gas",
    "low voltage": "AC",
    "oil": "Oil",
    "coal": "Coal",
    "lignite": "Coal",
    "municipal solid waste": "Waste",
    "solid biomass": "Biomass",
    "methanol": "Methanol",
    "non-sequestered HVC": "Waste",
    "uranium": "Uranium",
    "unsustainable bioliquids": "Liquids",
    # load buses
    "naphtha for industry": "Oil",
    "agriculture machinery oil": "Oil",
    "coal for industry": "Coal",
    "gas for industry": "Gas",
    "kerosene for aviation": "Oil",
    "land transport oil": "Oil",
    "shipping methanol": "Methanol",
    "shipping oil": "Oil",
    "solid biomass for industry": "Biomass",
    "EV battery": "AC",
    # Store
    "urban central water pits": "Heat",
    "urban central water tanks": "Heat",
}


class WriteOnceDict(dict):
    """Prevent overwriting existing keys."""

    def __setitem__(self, key: str, value: pd.Series):
        if key in self:
            raise KeyError(f"Key {key!r} already exists.")
        super().__setitem__(key, value)


def _process_single_input_link(
    var: dict,
    supply: pd.Series,
    demand: pd.Series,
    bc_out: pd.Index,
    bc_in: str,
    technology: str,
) -> dict:
    bc_in = BC_ALIAS.get(bc_in, bc_in)

    for bc in bc_out:
        label = f"{SECONDARY}|{BC_ALIAS.get(bc, bc)}|{bc_in}|{technology}"
        branch_output = filter_by(supply, bus_carrier=bc)
        var[label] = branch_output.groupby(IDX).sum()

    demand_unit = demand.index.unique("unit").item()
    losses = supply.groupby(YEAR_LOC).sum() + demand.groupby(YEAR_LOC).sum()
    losses = insert_index_level(losses, demand_unit, "unit", pos=2).mul(-1)

    var[f"{SECONDARY}|Losses|{bc_in}|{technology}"] = losses[losses > 0]
    losses_neg = losses[losses < 0]
    if not losses_neg.empty:
        # negative losses are expected for X-to-heat technologies that
        # utilize enthalpy heat and for heat pumps
        assert technology in ("CHP", "Boiler", "Ground Heat Pump", "Air Heat Pump"), (
            f"Unknown technology with efficiencies > 1: {technology}."
        )
        var[f"{SECONDARY}|Ambient Heat|{bc_in}|{technology}"] = rename_aggregate(
            losses_neg, "MWh_LHV", level="unit"
        ).mul(-1)

    # do not count ambient heat for demand + supply + losses = 0 to stay true
    bc_out_alias = [BC_ALIAS.get(bc, bc) for bc in bc_out] + ["Losses"]
    pattern = rf"{SECONDARY}\|({'|'.join(bc_out_alias)})\|{bc_in}\|{technology}"
    total_vars = (
        pd.concat({k: v for k, v in var.items() if re.match(pattern, k)})
        .groupby(YEAR_LOC)
        .sum()
    )
    # total demand + supply + losses - ambient heat = 0
    ambient_heat = var.get(
        f"{SECONDARY}|Ambient Heat|{bc_in}|{technology}", pd.Series()
    )
    if not ambient_heat.empty:
        assert (
            total_vars.add(demand.groupby(YEAR_LOC).sum()).sub(
                ambient_heat, fill_value=0
            )
            <= 1e-5
        ).all()
    else:
        assert (total_vars.add(demand.groupby(YEAR_LOC).sum()) <= 1e-5).all()

    return var


def transform_link(
    var, carrier: str | list, technology: str, debug: bool = False
) -> dict:
    """
    Transform a Link component into supply and transformation losses.

    The Link demand is equal to supply + losses and not included in
    the output to avoid redundant data. Losses have positive signs
    and the demand bus unit.

    Parameters
    ----------
    var
    carrier
    technology
    debug

    Returns
    -------
    :
    """
    supply = filter_by(SUPPLY, carrier=carrier, component="Link")
    demand = filter_by(DEMAND, carrier=carrier, component="Link")

    if supply.empty and demand.empty:
        logger.warning(
            f"No supply or demand found for {carrier}. Skipping transformation."
        )
        return var

    bc_in = demand.index.unique("bus_carrier")
    bc_out = supply.index.unique("bus_carrier")

    if len(bc_in) > 1:
        for bus_carrier_demand in bc_in:
            demand_bc = filter_by(demand, bus_carrier=bus_carrier_demand)
            demand_share = demand_bc.sum() / demand.sum()
            supply_bc = supply * demand_share
            var = _process_single_input_link(
                var,
                supply_bc,
                demand_bc,
                bc_out,
                bus_carrier_demand,
                technology,
            )
    else:
        var = _process_single_input_link(
            var, supply, demand, bc_out, bc_in.item(), technology
        )

    # optionally skipping global statistics update is useful during development
    if debug:
        return var

    # remove from global statistic to prevent double counting
    SUPPLY.drop(supply.index, inplace=True)
    # adding demand to IAMC is redundant, because supply + losses == demand,
    # and we the demand bus_carrier is known from the variable name.
    DEMAND.drop(demand.index, inplace=True)

    return var
    # for bc in bc_out:
    #     label = f"{SECONDARY}|{BC_ALIAS.get(bc, bc)}|{bc_in}|{technology}"
    #     branch_output = filter_by(supply, bus_carrier=bc)
    #     var[label] = branch_output.groupby(IDX).sum()
    #
    # losses_unit = demand.index.unique("unit").item()
    # losses = insert_index_level(losses, losses_unit, "unit", pos=2).mul(-1)
    #
    # var[f"{SECONDARY}|Losses|{bc_in}|{technology}"] = losses[losses > 0]
    # losses_neg = losses[losses < 0]
    # if not losses_neg.empty:
    #     # negative losses are expected for X-to-heat technologies that
    #     # utilize enthalpy heat and for heat pumps
    #     assert technology in ("CHP", "Boiler", "Ground Heat Pump", "Air Heat Pump"), (
    #         f"Unknown technology with efficiencies > 1: {technology}."
    #     )
    #     var[f"{SECONDARY}|Ambient Heat|{bc_in}|{technology}"] = rename_aggregate(
    #         losses_neg, "MWh_LHV", level="unit"
    #     ).mul(-1)
    #
    # # do not count ambient heat for demand + supply + losses = 0 to stay true
    # bc_out_alias = [BC_ALIAS.get(bc, bc) for bc in bc_out] + ["Losses"]
    # pattern = rf"{SECONDARY}\|({'|'.join(bc_out_alias)})\|{bc_in}\|{technology}"
    # total_vars = (
    #     pd.concat({k: v for k, v in var.items() if re.match(pattern, k)})
    #     .groupby(YEAR_LOC)
    #     .sum()
    # )
    # # total demand + supply + losses - ambient heat = 0
    # ambient_heat = var.get(
    #     f"{SECONDARY}|Ambient Heat|{bc_in}|{technology}", pd.Series()
    # )
    # if not ambient_heat.empty:
    #     assert (
    #         total_vars.add(demand.groupby(YEAR_LOC).sum()).sub(
    #             ambient_heat, fill_value=0
    #         )
    #         <= 1e-5
    #     ).all()
    # else:
    #     assert (total_vars.add(demand.groupby(YEAR_LOC).sum()) <= 1e-5).all()
    #
    # # optionally skipping global statistics update is useful during development
    # if debug:
    #     return var
    #
    # # remove from global statistic to prevent double counting
    # SUPPLY.drop(supply.index, inplace=True)
    # # adding demand to IAMC is redundant, because supply + losses == demand,
    # # and we the demand bus_carrier is known from the variable name.
    # DEMAND.drop(demand.index, inplace=True)

    # return var


def transform_load(var, carrier: str) -> dict:
    """

    Parameters
    ----------
    var
    carrier

    Returns
    -------
    :
    """

    if carrier.startswith(("shipping", "land transport")) or carrier.endswith(
        "aviation"
    ):
        sector = "Transport"
    elif carrier == "NH3":
        sector = "Non-energy usage"
    elif carrier.endswith("heat"):
        sector = "HH & Services"
    elif carrier.startswith("agriculture"):
        sector = "Agriculture"
    elif carrier == "electricity":
        sector = "Base Load"  # todo: sector load split
        # Base Load contains Transport, Industry, Households and service
    elif "industry" in carrier:
        sector = "Industry"
    else:
        raise ValueError(f"Unknown sector for Load carrier: {carrier}.")

    bc = {
        BC_ALIAS.get(bc, bc)
        for bc in filter_by(DEMAND, carrier=carrier, component="Load").index.unique(
            "bus_carrier"
        )
    }
    assert len(bc) == 1, (
        f"Mixed target bus carrier are not supported. "
        f"Found bus_carrier {bc} for carrier {carrier}."
    )
    bc = bc.pop()

    load_demand = _extract(DEMAND, carrier=carrier, component="Load")
    load_supply = _extract(SUPPLY, carrier=carrier, component="Load")
    # positive load values are possible if industry produces a surplus, for example
    load = load_demand.add(load_supply, fill_value=0)

    # # add carbon capture variants to correctly extract Links
    # carbon_capture_links = ("solid biomass for industry", "gas for industry")
    # if is_cc_variant := carrier in carbon_capture_links:
    #     carrier = [carrier, carrier + " CC"]

    supply = _extract(SUPPLY, carrier=carrier, component="Link")
    demand = _extract(DEMAND, carrier=carrier, component="Link")

    # if not load_supply.empty:
    #     logger.warning(
    #         f"Positive values in Load component detected for carrier: {carrier}.\n{load_supply.head()}\n"
    #         f"Please raise an issue in PyPSA and note down the issue number here."
    #         f"Positive Load values are added with positive signs to conserve balances."
    #     )
    #     load = load.add(load_supply, fill_value=0)

    # pd.testing.assert_series_equal(demand, load.groupby(YEAR_LOC).sum(), check_names=False)
    losses = supply.groupby(YEAR_LOC).sum().add(demand.groupby(YEAR_LOC).sum())
    # losses_exist = losses.abs().lt(1e-5).all()
    # if losses_exist and is_cc_variant:
    #     # CC Links have bus0 efficiencies < 1, i.e. they have losses
    #     pass
    # else:
    assert losses.abs().lt(1e-5).all(), (
        f"Supply and demand are not equal. Please check for losses and efficiencies != 1 for carrier: {carrier}.\n{losses}"
    )

    var[f"{FINAL}|{sector}|{bc}"] = load

    return var


#
# def _process_multi_input_link(var: dict, supply: pd.Series, demand: pd.Series, bc_out: pd.Index, bc_in: pd.Index, technology: str, debug: bool = False) -> dict:
#     """"""
#     # supply = filter_by(SUPPLY, carrier=carrier, component="Link")
#     # demand = filter_by(DEMAND, carrier=carrier, component="Link")
#     #
#     # if supply.empty and demand.empty:
#     #     logger.warning(
#     #         f"No supply or demand found for {carrier}. Skipping transformation."
#     #     )
#     #     return var
#     #
#     # bc_in = demand.index.unique("bus_carrier")
#     # bc_out = supply.index.unique("bus_carrier")
#
#     for bus_carrier_demand in bc_in:
#         demand_bc = filter_by(demand, bus_carrier=bus_carrier_demand)
#         demand_share = demand_bc.sum() / demand.sum()
#         supply_bc = supply * demand_share
#         var = _process_single_input_link(var, supply_bc, demand_bc, bc_out, bus_carrier_demand, technology, debug=debug)
#
#     return var


def _get_port_efficiency(substring: str, port: str, component: str = "Link"):
    """
    Return the efficiency for a component and port.

    Parameters
    ----------
    substring
    port
    component

    Returns
    -------
    :

    Raises
    ------
    NotImplementedError
        If the efficiency changes between years.

    ValueError
        If the substring filter yields more than one row.
    """
    _refining_efficiencies = set()
    for n in networks.values():
        _refining_efficiencies.add(
            port_efficiency(n, component, port).filter(like=substring).unique().item()
        )

    if len(_refining_efficiencies) != 1:
        raise NotImplementedError("Myopic efficiencies not supported.")

    return _refining_efficiencies.pop()


def _extract(ds: pd.Series, **filter_kwargs) -> pd.Series:
    """Extract and group filter results."""
    results = filter_by(ds, **filter_kwargs)
    ds.drop(results.index, inplace=True)
    return results.groupby(IDX).sum()


def process_biomass_boilers(var) -> dict:
    """

    Parameters
    ----------
    var

    Returns
    -------
    :
    """
    carrier = ["rural biomass boiler", "urban decentral biomass boiler"]
    if len(filter_by(DEMAND, carrier=carrier).index.unique("bus_carrier")) <= 1:
        transform_link(var, carrier, technology="Boiler")
        return var

    logger.warning(
        "Solid biomass boilers have negative values at heat buses. The applied workaround "
        "calculates balances to circumnavigate the bug. Please raise an issue at PyPSA-EUR "
        "and note down the issue number here."
    )
    balances = filter_by(LINK_BALANCE, carrier=carrier)
    # SUPPLY.drop(carrier, inplace=True, errors="ignore")
    # DEMAND.drop(carrier, inplace=True, errors="ignore")
    #
    # _supply = insert_index_level(balances[balances >= 0], "MWh_th", "unit", pos=5).pipe(insert_index_level, "Link", "component", pos=1)
    #
    #
    # SUPPLY.append(insert_index_level(balances[balances >= 0], "MWh_th", "unit", pos=5).pipe(insert_index_level, "Link", "component", pos=1))
    #
    # SUPPLY = pd.concat([SUPPLY, insert_index_level(balances[balances >= 0], "MWh_th", "unit", pos=5).pipe(insert_index_level, "Link", "component", pos=1)])
    # DEMAND = pd.concat([DEMAND, insert_index_level(balances[balances >= 0], "MWh_LHV", "unit", pos=5).pipe(insert_index_level, "Link", "component", pos=1)])
    # transform_link(var, technology="Boiler", carrier=carrier)

    var[f"{SECONDARY}|Heat|Biomass|Boiler"] = (
        balances.clip(lower=0)
        .pipe(insert_index_level, "MWH_th", "unit")
        .groupby(IDX)
        .sum()
    )
    _bal = insert_index_level(balances, "MWh_LHV", "unit").groupby(IDX).sum().mul(-1)
    losses = _bal[_bal.gt(0)]
    ambient_heat = _bal[_bal.le(0)].mul(-1)
    assert losses.gt(0).all()
    var[f"{SECONDARY}|Losses|Biomass|Boiler"] = losses
    if not ambient_heat.empty:
        assert (
            _bal.sub(losses, fill_value=0).sub(ambient_heat, fill_value=0) <= 1e-5
        ).all()
        var[f"{SECONDARY}|Ambient Heat|Biomass|Boiler"] = ambient_heat
    else:
        assert (_bal.sub(losses, fill_value=0) <= 1e-5).all()
    _extract(SUPPLY, carrier=carrier, component="Link")
    _extract(DEMAND, carrier=carrier, component="Link")

    return var


# def _get_traded_energy(n, var, bus_carrier, direction, subcat):
#     """
#     Calculate the trade statistics.
#
#     Parameters
#     ----------
#     n
#     var
#     bus_carrier
#     direction
#     subcat
#
#     Returns
#     -------
#     :
#     """
#     for scope in (TradeTypes.DOMESTIC, TradeTypes.FOREIGN):
#         trade = n.statistics.trade_energy(
#             bus_carrier=bus_carrier, direction=direction, scope=scope
#         )
#         trade = trade[trade.gt(0)].groupby("location").sum()
#         var[f"Primary Energy|{subcat}|{direction.title()} {scope.title()}"] = trade

#
# def _get_transformation_losses(supply: pd.Series, **filter_kwargs) -> pd.Series:
#     # single branch Links only
#     # unit_supply = supply.index.unique("unit").item()
#     # supply = supply.droplevel("unit")
#
#     demand = _extract(DEMAND, **filter_kwargs)
#     unit_demand = demand.index.unique("unit").item()
#     demand = demand.groupby(IDX).sum().droplevel("unit")
#
#     # unit_losses = unit_supply
#     # if unit_supply != unit_demand:
#     #     unit_losses = f"{unit_demand} to {unit_supply}"
#
#     losses = (
#         demand.add(supply.droplevel("unit"))
#         .pipe(insert_index_level, unit_demand, "unit")
#         .groupby(IDX)
#         .sum()
#     )
#
#     # implicitly True
#     assert demand.sum() + supply.sum() - losses.sum() <= 1e-5
#     # todo: myopic regional efficiencies
#     # eff = _get_port_efficiency(filter_kwargs["carrier"], port="1")
#     # assert abs(demand.sum() * eff + supply.sum()) <= 1e-4
#
#     return losses


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
    prefix = f"{PRIMARY}|Oil"
    bc = "oil"

    # assuming that all oil production is consumed locally.
    # Let's not filter_by components, to capture all but Stores.
    production = (
        filter_by(SUPPLY, bus_carrier=bc)
        .drop("EU", level="location")
        .groupby(IDX)
        .sum()
    )
    consumption = (
        filter_by(DEMAND, bus_carrier=bc)
        .drop("EU", level="location")
        .groupby(IDX)
        .sum()
    )
    assert "EU" not in consumption.index.unique("location")
    regional_oil_deficit = (
        consumption.add(production, fill_value=0).clip(upper=0).mul(-1)
    )
    regional_oil_surplus = consumption.add(production, fill_value=0).clip(lower=0)
    # assert primary oil and imports are equal to regional demands
    oil_refining_eff = _get_port_efficiency("oil refining", port="1")
    eu_oil = filter_by(SUPPLY, carrier=["oil refining", "import oil"], bus_carrier=bc)
    missing = (
        regional_oil_deficit.groupby("year").sum()
        - eu_oil.groupby("year").sum()
        - regional_oil_surplus.groupby("year").sum()
    )
    if not (missing <= 1e-5).abs().all():
        logger.warning(f"Missing oil amounts detected: {missing}")

    var[f"{prefix}|Fossil"] = regional_oil_deficit
    var[f"{prefix}|Refining Losses"] = regional_oil_deficit * (1 - oil_refining_eff)

    # remove EU imports and oil refining
    _extract(SUPPLY, carrier="import oil")
    _extract(SUPPLY, carrier="oil primary")
    _extract(SUPPLY, carrier="oil refining")
    _extract(DEMAND, carrier="oil refining")

    # unsustainable bioliquids have regional bus generators, but are already
    # accounted for in regional_oil_deficit. The "unsustainable bioliquids"
    # Link supplies to the oil bus.
    _extract(SUPPLY, carrier="unsustainable bioliquids", component="Generator")

    # var["Primary Energy|Oil"] = var["Primary Energy|Oil|Fossil Oil"].add(
    #     var["Primary Energy|Oil|Unsustainable Bioliquids"], fill_value=0
    # )
    var[prefix] = _sum_variables_by_prefix(var, prefix)

    return var


def primary_gas(var) -> dict:
    bc = "gas"
    prefix = f"{PRIMARY}|{BC_ALIAS.get(bc, bc)}"
    var[f"{prefix}|Import Foreign"] = _extract(IMPORT_FOREIGN, bus_carrier=bc)
    var[f"{prefix}|Import Domestic"] = _extract(IMPORT_DOMESTIC, bus_carrier=bc)
    var[f"{prefix}|Production"] = _extract(
        SUPPLY, carrier="gas", bus_carrier=bc, component="Generator"
    )
    var[f"{prefix}|Import Global"] = _extract(
        SUPPLY,
        carrier="import gas",
        bus_carrier=bc,
        component="Generator",
    )

    # todo: move to primary - this link is only needed to track CO2
    var[f"{prefix}|Biogas|w/o CC"] = _extract(
        SUPPLY, carrier="biogas to gas", bus_carrier=bc
    )
    var[f"{prefix}|Biogas|w CC"] = _extract(
        SUPPLY, carrier="biogas to gas CC", bus_carrier=bc
    )
    _extract(
        DEMAND, carrier=["biogas to gas", "biogas to gas CC"], bus_carrier="biogas"
    )

    var[prefix] = _sum_variables_by_prefix(var, prefix)

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
    bc = "municipal solid waste"
    prefix = f"{PRIMARY}|{BC_ALIAS.get(bc, bc)}"
    var[f"{prefix}|Import Foreign"] = _extract(IMPORT_FOREIGN, bus_carrier=bc)
    var[f"{prefix}|Import Domestic"] = _extract(IMPORT_DOMESTIC, bus_carrier=bc)
    var[f"{prefix}|Solid"] = _extract(SUPPLY, bus_carrier=bc, component="Generator")

    # municipal solid waste is only used to transform "municipal solid waste" to
    # "non-sequestered HVC" and to track CO2. Same as Biogas, include in primary
    _extract(
        SUPPLY,
        carrier="municipal solid waste",
        bus_carrier="non-sequestered HVC",
        component="Link",
    )
    _extract(DEMAND, carrier="municipal solid waste", bus_carrier=bc, component="Link")

    # HVC is a side product of naphtha for industry. The oil demand of
    # the link equals the naphtha output. There are no losses.
    var[f"{prefix}|HVC from naphtha processing"] = _extract(
        SUPPLY,
        carrier="naphtha for industry",
        bus_carrier="non-sequestered HVC",
        component="Link",
    )

    var[prefix] = _sum_variables_by_prefix(var, prefix)

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
    prefix = f"{PRIMARY}|Coal"
    var[f"{prefix}|Hard"] = (
        filter_by(DEMAND, bus_carrier="coal", component="Link").groupby(IDX).sum()
    ).mul(-1)
    var[f"{prefix}|Lignite"] = (
        filter_by(DEMAND, bus_carrier="lignite", component="Link").groupby(IDX).sum()
    ).mul(-1)
    var[prefix] = _sum_variables_by_prefix(var, prefix)

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
    bc = "H2"
    prefix = f"{PRIMARY}|{BC_ALIAS.get(bc, bc)}"
    var[f"{prefix}|Import Foreign"] = _extract(IMPORT_FOREIGN, bus_carrier=bc)
    var[f"{prefix}|Import Domestic"] = _extract(IMPORT_DOMESTIC, bus_carrier=bc)
    var[f"{prefix}|Import Global"] = _extract(
        SUPPLY, carrier="import H2", bus_carrier=bc, component="Generator"
    )

    var[prefix] = _sum_variables_by_prefix(var, prefix)

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
    bc = "solid biomass"
    prefix = f"{PRIMARY}|{BC_ALIAS.get(bc, bc)}"
    var[f"{prefix}|Import Foreign"] = _extract(IMPORT_FOREIGN, bus_carrier=bc)
    var[f"{prefix}|Import Domestic"] = _extract(IMPORT_DOMESTIC, bus_carrier=bc)

    var[f"{prefix}|Solid"] = _extract(SUPPLY, bus_carrier=bc, component="Generator")
    var[f"{prefix}|Biogas"] = _extract(
        SUPPLY, bus_carrier="biogas", component="Generator"
    )
    var[prefix] = _sum_variables_by_prefix(var, prefix)

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
    prefix = f"{PRIMARY}|Hydro"
    # var["Primary Energy|Hydro|PHS"] = _extract(
    #     SUPPLY, carrier="PHS", component="StorageUnit"
    # )
    var[f"{prefix}|Reservoir"] = _extract(
        SUPPLY, carrier="hydro", component="StorageUnit"
    )
    var[f"{prefix}|Run-of-River"] = _extract(
        SUPPLY, carrier="ror", component="Generator"
    )
    var[prefix] = _sum_variables_by_prefix(var, prefix)

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
    prefix = f"{PRIMARY}|Solar"
    var[f"{prefix}|Utility"] = _extract(SUPPLY, carrier="solar", component="Generator")
    var[f"{prefix}|HSAT"] = _extract(
        SUPPLY, carrier="solar-hsat", component="Generator"
    )
    var[f"{prefix}|Rooftop"] = _extract(
        SUPPLY, carrier="solar rooftop", component="Generator"
    )
    var[prefix] = _sum_variables_by_prefix(var, prefix)

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
    # var[f"{PRIMARY}|Oil|Unsustainable Bioliquids"] = _extract(
    #     SUPPLY, carrier="unsustainable bioliquids", component="Link"
    # )
    # # remove liquids demand
    # _extract(DEMAND, carrier="unsustainable bioliquids", component="Link")

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
    # use localized uranium demands from nuclear power plants
    var[f"{PRIMARY}|Nuclear|Uranium"] = (
        filter_by(DEMAND, carrier="uranium", component="Link")
        .groupby(IDX)
        .sum()
        .mul(-1)
    )
    _extract(SUPPLY, bus_carrier="uranium", component="Generator")

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
    var[f"{PRIMARY}|Ammonium|Import"] = _extract(SUPPLY, carrier="import NH3")

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
    var[f"{PRIMARY}|Methanol|Import"] = _extract(SUPPLY, carrier="import methanol")

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
    prefix = f"{PRIMARY}|Wind"
    var[f"{prefix}|Onshore"] = _extract(SUPPLY, carrier="onwind", component="Generator")
    var[f"{prefix}|Offshore"] = _extract(
        SUPPLY, carrier=["offwind-ac", "offwind-dc"], component="Generator"
    )
    var[prefix] = _sum_variables_by_prefix(var, prefix)

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
    # heat_bus_carrier = ["rural heat", "urban decentral heat", "urban central heat"]
    prefix = f"{PRIMARY}|{BC_ALIAS.get('rural heat', 'Heat')}"

    # done in secondary when processing Links
    # enthalpy_heat = (
    #     LINK_BALANCE.pipe(filter_for_carrier_connected_to, heat_bus_carrier)
    #     .pipe(calculate_input_share, heat_bus_carrier)
    #     .pipe(filter_by, bus_carrier=["ambient heat", "latent heat"])
    #     .pipe(insert_index_level, "MWh_th", "unit")
    # )
    #
    # var[f"{prefix}|Latent"] = (
    #     filter_by(enthalpy_heat, bus_carrier="latent heat").groupby(IDX).sum()
    # )
    # var[f"{prefix}|Ambient"] = (
    #     filter_by(enthalpy_heat, bus_carrier="ambient heat").groupby(IDX).sum()
    # )

    solar_thermal_carr = [
        c for c in SUPPLY.index.unique("carrier") if "solar thermal" in c
    ]
    var[f"{prefix}|Solar"] = _extract(
        SUPPLY, carrier=solar_thermal_carr, component="Generator"
    )
    var[f"{prefix}|Geothermal"] = _extract(SUPPLY, carrier="geothermal heat")
    var[prefix] = _sum_variables_by_prefix(var, prefix)

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
    var = WriteOnceDict()
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


#
# def secondary_electricity_supply(var: dict) -> dict:
#     """
#
#     Parameters
#     ----------
#     var
#
#     Returns
#     -------
#     :
#     """
#     prefix = "Secondary Energy|Electricity"
#     bc = ["AC", "low voltage"]
#
#     transform_link(var, technology="CHP", carrier="urban central gas CHP")
#     transform_link(var, technology="CHP", carrier="urban central oil CHP")
#     transform_link(var, technology="CHP", carrier="urban central coal CHP")
#     transform_link(var, technology="CHP", carrier="urban central lignite CHP")
#     transform_link(
#         var,
#         technology="CHP",
#         carrier=["urban central H2 CHP", "urban central H2 retrofit CHP"],
#     )
#     transform_link(var, technology="CHP", carrier="urban central solid biomass CHP")
#     transform_link(var, technology="CHP w/o CC", carrier="waste CHP")
#     transform_link(var, technology="CHP w CC", carrier="waste CHP CC")
#
#     transform_link(var, technology="Powerplant", carrier=["CCGT", "OCGT"])
#     transform_link(var, technology="Powerplant", carrier="coal")
#     transform_link(var, technology="Powerplant", carrier="lignite")
#     transform_link(var, technology="Powerplant", carrier="solid biomass")
#     transform_link(var, technology="Powerplant", carrier="nuclear")
#
#     transform_link(var, technology="", carrier="SMR")
#
#     # var[f"{prefix}|Oil|CHP"] = _extract(
#     #     SUPPLY, carrier="urban central oil CHP", bus_carrier=bc
#     # )
#     # var[f"{prefix}|Coal|CHP Hard Coal"] = _extract(
#     #     SUPPLY, carrier="urban central coal CHP", bus_carrier=bc
#     # )
#     # var[f"{prefix}|Coal|Hard Coal"] = _extract(SUPPLY, carrier="coal", bus_carrier=bc)
#     # var[f"{prefix}|Coal|Lignite"] = _extract(SUPPLY, carrier="lignite", bus_carrier=bc)
#     # var[f"{prefix}|Coal|CHP Lignite"] = _extract(
#     #     SUPPLY, carrier="urban central lignite CHP", bus_carrier=bc
#     # )
#     # var[f"{prefix}|H2|CHP"] = _extract(
#     #     SUPPLY,
#     #     carrier=["urban central H2 CHP", "urban central H2 retrofit CHP"],
#     #     bus_carrier=bc,
#     # )
#     # var[f"{prefix}|Solid Biomass"] = _extract(
#     #     SUPPLY,
#     #     carrier=["solid biomass", "urban central solid biomass CHP"],
#     #     bus_carrier=bc,
#     # )
#
#     # var[f"{prefix}|Nuclear"] = _extract(SUPPLY, carrier="nuclear", bus_carrier=bc)
#
#     # var[f"{prefix}|Waste|CHP w/o CC"] = _extract(
#     #     SUPPLY, carrier="waste CHP", bus_carrier=bc
#     # )
#     # var[f"{prefix}|Waste|CHP w CC"] = _extract(
#     #     SUPPLY, carrier="waste CHP CC", bus_carrier=bc
#     # )
#
#     var[prefix] = _sum_variables_by_prefix(var, prefix)
#
#     # subcategory aggregation must happen after prefix aggregations, or the
#     # sum will be distorted
#     for subcat in ("Gas", "Coal", "Waste"):
#         var[f"{prefix}|{subcat}"] = _sum_variables_by_prefix(var, f"{prefix}|{subcat}")
#
#     # distribution grid losses are no supply, but we deal with it now remove all
#     # electricity from the global supply statistic
#     var["Secondary Energy|Losses|Electricity|Distribution Grid"] = (
#         _extract(  # todo: label and move to losses
#             SUPPLY, carrier="electricity distribution grid"
#         )
#         + _extract(DEMAND, carrier="electricity distribution grid")
#     )
#
#     assert filter_by(SUPPLY, bus_carrier=bc, component="Link").empty
#
#     return var


# def secondary_gas_supply(var: dict) -> dict:
#     """
#
#     Parameters
#     ----------
#     var
#
#     Returns
#     -------
#     :
#     """
#     # Secondary Energy|<output bus_carrier>|from <input bus_carrier>|<subcategory>
#     prefix = "Secondary Energy|Gas"
#     bc = "gas"
#     var[f"{prefix}|Biogas|w/o CC"] = _extract(
#         SUPPLY, carrier="biogas to gas", bus_carrier=bc
#     )
#     var[f"{prefix}|Biogas|w CC"] = _extract(
#         SUPPLY, carrier="biogas to gas CC", bus_carrier=bc
#     )
#     # _extract(DEMAND, carrier=["biogas to gas", "biogas to gas CC"], bus_carrier="biogas")
#
#     var[f"{prefix}|Solid Biomass|w/o CC"] = _extract(
#         SUPPLY, carrier="BioSNG", bus_carrier=bc
#     )
#     var[f"{prefix}|Solid Biomass|w CC"] = _extract(
#         SUPPLY, carrier="BioSNG CC", bus_carrier=bc
#     )
#     # _extract(DEMAND, carrier=["BioSNG", "BioSNG CC"], bus_carrier="solid biomass")
#
#     var[f"{prefix}|Sabatier"] = _extract(SUPPLY, carrier="Sabatier", bus_carrier=bc)
#
#     var[prefix] = _sum_variables_by_prefix(var, prefix)
#
#     assert filter_by(SUPPLY, bus_carrier=bc, component="Link").empty
#
#     # move to Carbon function
#     # var["Secondary Carbon|Biogas|atmosphere"] = _extract(DEMAND, carrier="biogas to gas", bus_carrier="co2")  # CO2 credit
#     # var["Secondary Carbon|Biogas|stored"] = _extract(SUPPLY, carrier="biogas to gas CC", bus_carrier="co2")  # carbon capture
#
#     return var
#
#
# def secondary_hydrogen_supply(var: dict) -> dict:
#     """
#
#     Parameters
#     ----------
#     var
#
#     Returns
#     -------
#     :
#     """
#     prefix = "Secondary Energy|Hydrogen"
#     bc = "H2"
#     var[f"{prefix}|Electricity|Electrolysis"] = _extract(
#         SUPPLY, carrier="H2 Electrolysis", bus_carrier=bc
#     )
#     var[f"{prefix}|Gas|SMR w/o CC"] = _extract(SUPPLY, carrier="SMR", bus_carrier=bc)
#     var[f"{prefix}|Gas|SMR w CC"] = _extract(SUPPLY, carrier="SMR CC", bus_carrier=bc)
#     var[f"{prefix}|Methanol|Steam Reforming w/o CC"] = _extract(
#         SUPPLY, carrier="Methanol steam reforming", bus_carrier=bc
#     )
#     var[f"{prefix}|Methanol|Steam Reforming w CC"] = _extract(
#         SUPPLY, carrier="Methanol steam reforming CC", bus_carrier=bc
#     )
#
#     var[prefix] = _sum_variables_by_prefix(var, prefix)
#
#     assert filter_by(SUPPLY, bus_carrier=bc, component="Link").empty
#
#     return var
#
#
# def secondary_methanol_supply(var: dict) -> dict:
#     """
#
#     Parameters
#     ----------
#     var
#
#     Returns
#     -------
#     :
#     """
#     prefix = "Secondary Energy|Methanol"
#     bc = "methanol"
#
#     # need to distinguish between methanol and heat output
#     methanolisation_inputs_for_methanol = (
#         LINK_BALANCE.pipe(filter_for_carrier_connected_to, bc)
#         .pipe(calculate_input_share, bc)
#         .pipe(filter_by, carrier="methanolisation")
#     )
#     var[f"{prefix}|H2"] = (
#         filter_by(methanolisation_inputs_for_methanol, bus_carrier="H2")
#         .pipe(insert_index_level, "MWh_LHV", "unit")
#         .groupby(IDX)
#         .sum()
#     )
#     var[f"{prefix}|Electricity"] = (
#         filter_by(methanolisation_inputs_for_methanol, bus_carrier="AC")
#         .pipe(insert_index_level, "MWh_el", "unit")
#         .groupby(IDX)
#         .sum()
#     )
#     _extract(SUPPLY, carrier="methanolisation", bus_carrier=bc)
#
#     var[prefix] = _sum_variables_by_prefix(var, prefix)
#
#     assert filter_by(SUPPLY, bus_carrier=bc, component="Link").empty
#     # biomass-to-methanol?
#
#     return var
#
#
# def secondary_oil_supply(var: dict) -> dict:
#     """
#
#     Parameters
#     ----------
#     var
#
#     Returns
#     -------
#     :
#     """
#     prefix = "Secondary Energy|Oil"
#     bc = "oil"
#
#     var[f"{prefix}|Solid Biomass|Biomass2Liquids w/o CC"] = _extract(
#         SUPPLY, carrier="biomass to liquid", bus_carrier=bc
#     )
#     var[f"{prefix}|Solid Biomass|Biomass2Liquids w CC"] = _extract(
#         SUPPLY, carrier="biomass to liquid CC", bus_carrier=bc
#     )
#
#     # electrobiofuels has 2 inputs: solid biomass and H2 and one output
#     electrobiofuels_inputs_for_oil = (
#         LINK_BALANCE.pipe(filter_for_carrier_connected_to, bc)
#         .pipe(calculate_input_share, bc)
#         .pipe(filter_by, carrier="electrobiofuels")
#     )
#     var[f"{prefix}|H2|Electrobiofuels"] = (
#         filter_by(electrobiofuels_inputs_for_oil, bus_carrier="H2")
#         .pipe(insert_index_level, "MWh_LHV", "unit")
#         .groupby(IDX)
#         .sum()
#     )
#     var[f"{prefix}|Electricity|Electrobiofuels"] = (
#         filter_by(electrobiofuels_inputs_for_oil, bus_carrier="solid biomass")
#         .pipe(insert_index_level, "MWh_LHV", "unit")
#         .groupby(IDX)
#         .sum()
#     )
#     _extract(SUPPLY, carrier="electrobiofuels", bus_carrier=bc)
#
#     var[f"{prefix}|H2|Fischer-Tropsch"] = _extract(
#         SUPPLY, carrier="Fischer-Tropsch", bus_carrier=bc
#     )
#
#     var[f"{prefix}|Unsustainable Bioliquids"] = _extract(
#         SUPPLY, carrier="unsustainable bioliquids", bus_carrier=bc
#     )
#
#     assert filter_by(SUPPLY, bus_carrier=bc, component="Link").empty
#
#     var[prefix] = _sum_variables_by_prefix(var, prefix)
#     for group in ("H2", "Electricity", "Solid Biomass"):
#         var[f"{prefix}|{group}"] = _sum_variables_by_prefix(var, f"{prefix}|{group}")
#
#     return var
#
#
# def secondary_ammonia_supply(var: dict) -> dict:
#     """
#
#     Parameters
#     ----------
#     var
#
#     Returns
#     -------
#     :
#     """
#     prefix = "Secondary Energy|NH3"
#     bc = "NH3"
#
#     haber_bosch_input_for_nh3 = (
#         LINK_BALANCE.pipe(filter_for_carrier_connected_to, bc)
#         .pipe(calculate_input_share, bc)
#         .pipe(filter_by, carrier="Haber-Bosch")
#         .pipe(insert_index_level, "n/a", "unit")
#     )
#     var[f"{prefix}|Electricity|Haber-Bosch"] = (
#         filter_by(haber_bosch_input_for_nh3, carrier="Haber-Bosch", bus_carrier="AC")
#         .groupby(IDX)
#         .sum()
#     )
#     var[f"{prefix}|H2|Haber-Bosch"] = (
#         filter_by(haber_bosch_input_for_nh3, carrier="Haber-Bosch", bus_carrier="H2")
#         .groupby(IDX)
#         .sum()
#     )
#     _extract(SUPPLY, carrier="Haber-Bosch", bus_carrier=bc)
#
#     return var
#
#
# def secondary_heat_supply(var: dict) -> dict:
#     """
#
#     Parameters
#     ----------
#     var
#
#     Returns
#     -------
#     :
#     """
#     prefix = "Secondary Energy|Heat"
#     bc = BusCarrier.heat_buses()
#
#     var[f"{prefix}|Solid Biomass|Boiler"] = _extract(
#         SUPPLY,
#         carrier=["rural biomass boiler", "urban decentral biomass boiler"],
#         bus_carrier=bc,
#     )
#     var[f"{prefix}|Solid Biomass|CHP"] = _extract(
#         SUPPLY, carrier="urban central solid biomass CHP", bus_carrier=bc
#     )
#
#     var[f"{prefix}|Electricity|Ground Heat Pump"] = _extract(
#         SUPPLY, carrier="rural ground heat pump", bus_carrier=bc
#     )
#     var[f"{prefix}|Electricity|Air Heat Pump"] = _extract(
#         SUPPLY,
#         carrier=[
#             "urban decentral air heat pump",
#             "rural air heat pump",
#             "urban central air heat pump",
#         ],
#         bus_carrier=bc,
#     )
#     var[f"{prefix}|Gas|Boiler"] = _extract(
#         SUPPLY,
#         carrier=[
#             "rural gas boiler",
#             "urban central gas boiler",
#             "urban decentral gas boiler",
#         ],
#         bus_carrier=bc,
#     )
#
#     var[f"{prefix}|Electricity|Resistive Heater"] = _extract(
#         SUPPLY,
#         carrier=[
#             "rural resistive heater",
#             "urban decentral resistive heater",
#             "urban central resistive heater",
#         ],
#         bus_carrier=bc,
#     )
#
#     var[f"{prefix}|Oil|Boiler"] = _extract(
#         SUPPLY,
#         carrier=["rural oil boiler", "urban decentral oil boiler"],
#         bus_carrier=bc,
#     )
#
#     var[f"{prefix}|H2|CHP"] = _extract(
#         SUPPLY,
#         carrier=["urban central H2 CHP", "urban central H2 retrofit CHP"],
#         bus_carrier=bc,
#     )
#
#     var[f"{prefix}|Gas|CHP"] = _extract(
#         SUPPLY, carrier="urban central gas CHP", bus_carrier=bc
#     )
#
#     var[f"{prefix}|Oil|CHP"] = _extract(
#         SUPPLY,
#         carrier="urban central oil CHP",
#         bus_carrier=bc,
#     )
#
#     var[f"{prefix}|Coal|CHP"] = _extract(
#         SUPPLY,
#         carrier=["urban central coal CHP", "urban central lignite CHP"],
#         bus_carrier=bc,
#     )
#
#     var[f"{prefix}|Waste|CHP w/o CC"] = _extract(
#         SUPPLY, carrier="waste CHP", bus_carrier=bc
#     )
#     var[f"{prefix}|Waste|CHP w CC"] = _extract(
#         SUPPLY, carrier="waste CHP CC", bus_carrier=bc
#     )
#
#     # Excess heat from technologies
#     var[f"{prefix}|H2|Sabatier"] = _extract(SUPPLY, carrier="Sabatier", bus_carrier=bc)
#     var[f"{prefix}|H2|Fischer-Tropsch"] = _extract(
#         SUPPLY, carrier="Fischer-Tropsch", bus_carrier=bc
#     )
#     var[f"{prefix}|Electricity|Electrolysis"] = _extract(
#         SUPPLY, carrier="H2 Electrolysis", bus_carrier=bc
#     )
#
#     # methanolisation and Haber-Bosch have multiple inputs
#     methanolisation_input_for_heat = (
#         LINK_BALANCE.drop(["co2", "co2 stored"], level="bus_carrier")
#         .pipe(filter_for_carrier_connected_to, bc)
#         .pipe(calculate_input_share, bc)
#         .pipe(filter_by, carrier=["methanolisation", "Haber-Bosch"])
#         .pipe(insert_index_level, "MWh_LHV", "unit")
#     )
#     var[f"{prefix}|Electricity|Methanolisation"] = (
#         filter_by(
#             methanolisation_input_for_heat, carrier="methanolisation", bus_carrier="AC"
#         )
#         .groupby(IDX)
#         .sum()
#     )
#     var[f"{prefix}|H2|Methanolisation"] = (
#         filter_by(
#             methanolisation_input_for_heat, carrier="methanolisation", bus_carrier="H2"
#         )
#         .groupby(IDX)
#         .sum()
#     )
#     _extract(SUPPLY, carrier="methanolisation", bus_carrier=bc)
#
#     var[f"{prefix}|Electricity|Haber-Bosch"] = (
#         filter_by(
#             methanolisation_input_for_heat, carrier="Haber-Bosch", bus_carrier="AC"
#         )
#         .groupby(IDX)
#         .sum()
#     )
#     var[f"{prefix}|H2|Haber-Bosch"] = (
#         filter_by(
#             methanolisation_input_for_heat, carrier="Haber-Bosch", bus_carrier="H2"
#         )
#         .groupby(IDX)
#         .sum()
#     )
#     _extract(SUPPLY, carrier="Haber-Bosch", bus_carrier=bc)
#
#     # deal with water pit charger losses now to clear all heat buses
#     var["Storage Losses|Heat|Water Pits"] = _extract(  # todo: label and move to storage
#         SUPPLY,
#         carrier="urban central water pits discharger",
#         bus_carrier="urban central heat",
#     ) + _extract(
#         DEMAND,
#         carrier="urban central water pits charger",
#         bus_carrier="urban central heat",
#     )
#     # drop the supply/demand at the other bus side of (dis)charger links
#     _extract(
#         SUPPLY,
#         carrier="urban central water pits charger",
#         bus_carrier="urban central water pits",
#     )
#     _extract(
#         DEMAND,
#         carrier="urban central water pits discharger",
#         bus_carrier="urban central water pits",
#     )
#
#     assert filter_by(SUPPLY, bus_carrier=bc, component="Link").empty
#
#     return var
#
#
# def secondary_waste_supply(var: dict) -> dict:
#     prefix = "Secondary Energy|Waste"
#     bc = "non-sequestered HVC"
#
#     var[f"{prefix}|Waste|Naptha"] = _extract(
#         SUPPLY, carrier="naphtha for industry", bus_carrier=bc
#     )
#     var[f"{prefix}|Waste|Municipal solid waste"] = _extract(
#         SUPPLY, carrier="municipal solid waste", bus_carrier=bc
#     )
#
#     assert filter_by(SUPPLY, bus_carrier=bc, component="Link").empty
#
#     var[prefix] = _sum_variables_by_prefix(var, prefix)
#
#     return var
#
#
# def secondary_solid_biomass_supply(var: dict) -> dict:
#     var["Secondary Energy|Solid Biomass|Boiler Error"] = _extract(
#         SUPPLY,
#         carrier=["rural biomass boiler", "urban decentral biomass boiler"],
#         bus_carrier="solid biomass",
#     )
#
#     if not var["Secondary Energy|Solid Biomass|Boiler Error"].empty:
#         logger.warning(
#             f"Solid biomass boilers supply to solid biomass bus. Total amount of energy supplied = {var['Secondary Energy|Solid Biomass|Boiler Error'].sum():.2f} MWh_LHV"
#         )
#
#     assert filter_by(SUPPLY, bus_carrier="solid biomass", component="Link").empty
#
#     return var


# def secondary_electricity_losses(var: dict) -> dict:
#     # needs secondary energy var
#     # needs a robust way to identify labels and carrier
#     # needs a way to test, because subtraction of supply and demand is error-prone  --> get_efficiency grouper for Links
#     # simply use LINK_BALANCE ? with carrier filter
#
#     # regex to find secondary energy supply by output group
#     # re.match("Secondary Energy\|[a-z, A-Z]*\|Electricity\|Electrolysis")
#
#     prefix = "Secondary Energy|Losses|Electricity"  # stay under Secondary Energy!
#     bc = ["AC", "low voltage"]
#
#     losses = (
#         filter_by(LINK_BALANCE, carrier="H2 Electrolysis")
#         .drop(["co2", "co2 stored"], level="bus_carrier", errors="ignore")
#         .groupby(["year", "location"])
#         .sum()
#     )
#     demand = (
#         filter_by(DEMAND, carrier="H2 Electrolysis").groupby(["year", "location"]).sum()
#     )
#     supply = (
#         var["Secondary Energy|Hydrogen|Electricity|Electrolysis"]
#         .droplevel("unit")
#         .add(
#             var["Secondary Energy|Heat|Electricity|Electrolysis"].droplevel("unit"),
#             fill_value=0,
#         )
#     )
#
#     assert (losses.sub(demand + supply) < 0.000001).all()
#
#     var["Transformation Losses|Electricity|Electrolysis"] = _get_transformation_losses(
#         var["Secondary Energy|Hydrogen|Electricity|Electrolysis"]
#         + var["Secondary Energy|Heat|Electricity|Electrolysis"],
#         carrier="H2 Electrolysis",
#         bus_carrier=bc,
#     )
#
#     # heat pumps Ground Heat Pump: no losses but ambient heat
#     var[f"{prefix}|Ground Heat Pump"] = _get_transformation_losses(
#         var["Secondary Energy|Heat|Electricity|Ground Heat Pump"],
#         carrier="rural ground heat pump",
#         bus_carrier=bc,
#     )
#
#     # heat pumps Air Heat Pump
#     var[f"{prefix}|Air Heat Pump"] = _get_transformation_losses(
#         var["Secondary Energy|Heat|Electricity|Air Heat Pump"],
#         carrier=[""],
#         bus_carrier=bc,
#     )
#
#     # resistive heater
#
#     # Secondary Energy|Heat|Solid Biomass|CHP
#     # Secondary Energy|AC|Solid Biomass|CHP
#     # Secondary Energy|*|Solid Biomass|CHP    -> supply regex
#
#     # Transformation Losses|Solid Biomass|CHP  --> losses label
#
#     assert filter_by(DEMAND, bus_carrier=bc, component="Link").empty
#
#     return var


def collect_storage_imbalances() -> pd.Series:
    """Extract all storage imbalances due to losses."""
    var = WriteOnceDict()

    for carrier in filter_by(SUPPLY, component="Store").index.unique("carrier"):
        supply = filter_by(SUPPLY, component="Store", carrier=carrier)
        demand = filter_by(DEMAND, component="Store", carrier=carrier)
        balance = supply.add(demand, fill_value=0)

        if balance.sum() != 0:
            logger.warning(
                f"Store imbalances detected for carrier {carrier} with total imbalance of {balance.groupby('year').sum()}."
            )
            bc = balance.index.unique("bus_carrier").item()
            if carrier == "urban central water pits":
                tech = "Water Pits"
            elif carrier == "urban central water tanks":
                tech = "Water Tank"
            elif carrier == "coal":
                tech = "Coal"
            else:
                raise ValueError(f"Unknown carrier: {carrier}")
            var[f"{SECONDARY}|Losses|{BC_ALIAS[bc]}|{tech}"] = balance.groupby(
                IDX
            ).sum()
        else:
            logger.debug(f"No Store imbalances detected for carrier: {carrier}.")

        SUPPLY.drop(supply.index, inplace=True)
        DEMAND.drop(demand.index, inplace=True)

    return combine_variables(var)


def collect_losses_energy() -> pd.Series:
    # DSM Links with losses that connect to buses with stores
    var = WriteOnceDict()
    prefix = f"{SECONDARY}|Losses"

    var[f"{prefix}|AC|Distribution Grid"] = _extract(
        SUPPLY, carrier="electricity distribution grid"
    ) + _extract(DEMAND, carrier="electricity distribution grid")

    var[f"{prefix}|AC|BEV charger"] = (
        _extract(SUPPLY, carrier="BEV charger", component="Link")
        .add(_extract(DEMAND, carrier="BEV charger", component="Link"))
        .mul(-1)
    )
    var[f"{prefix}|Heat|Water Pits"] = _extract(
        SUPPLY,
        carrier="urban central water pits discharger",
        bus_carrier="urban central heat",
        component="Link",
    ).add(
        _extract(
            DEMAND,
            carrier="urban central water pits charger",
            bus_carrier="urban central heat",
            component="Link",
        )
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

    # DAC has no outputs but CO2, which is ignored in energy flows
    var[f"{prefix}|AC|DAC"] = _extract(DEMAND, carrier="DAC", bus_carrier="AC").mul(-1)
    var[f"{prefix}|Heat|DAC"] = _extract(
        DEMAND,
        carrier="DAC",
        bus_carrier=["rural heat", "urban decentral heat", "urban central heat"],
    ).mul(-1)
    var[f"{prefix}|Waste|HVC to air"] = _extract(
        DEMAND,
        carrier="HVC to air",
        component="Link",
        bus_carrier="non-sequestered HVC",
    ).mul(-1)

    return combine_variables(var)


def collect_secondary_energy() -> pd.Series:
    """Extract all secondary energy variables from the networks."""
    var = WriteOnceDict()

    # secondary_electricity_supply(var)
    transform_link(var, technology="CHP", carrier="urban central gas CHP")
    transform_link(var, technology="CHP", carrier="urban central oil CHP")
    transform_link(var, technology="CHP", carrier="urban central coal CHP")
    transform_link(var, technology="CHP", carrier="urban central lignite CHP")
    transform_link(
        var,
        technology="CHP",
        carrier=["urban central H2 CHP", "urban central H2 retrofit CHP"],
    )
    transform_link(var, technology="CHP", carrier="urban central solid biomass CHP")
    transform_link(var, technology="CHP w/o CC", carrier="waste CHP")
    transform_link(var, technology="CHP w CC", carrier="waste CHP CC")

    transform_link(var, technology="Powerplant", carrier=["CCGT", "OCGT"])
    transform_link(var, technology="Powerplant", carrier="coal")
    transform_link(var, technology="Powerplant", carrier="lignite")
    transform_link(var, technology="Powerplant", carrier="solid biomass")
    transform_link(var, technology="Powerplant", carrier="nuclear")

    transform_link(var, technology="BioSNG w/o CC", carrier="BioSNG")
    transform_link(var, technology="BioSNG w CC", carrier="BioSNG CC")
    transform_link(var, technology="Sabatier", carrier="Sabatier")

    transform_link(var, technology="Electrolysis", carrier="H2 Electrolysis")
    transform_link(var, technology="SMR w/o CC", carrier="SMR")
    transform_link(var, technology="SMR w CC", carrier="SMR CC")
    transform_link(
        var, technology="Steam Reforming w/o CC", carrier="Methanol steam reforming"
    )
    transform_link(
        var, technology="Steam Reforming w CC", carrier="methanol steam reforming CC"
    )
    transform_link(var, technology="SMR w CC", carrier="SMR CC")

    transform_link(
        var, technology="Biomass2Liquids w/o CC", carrier="biomass to liquid"
    )
    transform_link(
        var, technology="Biomass2Liquids w CC", carrier="biomass to liquid CC"
    )
    transform_link(var, technology="Fischer-Tropsch", carrier="Fischer-Tropsch")
    transform_link(
        var, technology="Unsustainable Bioliquids", carrier="unsustainable bioliquids"
    )

    transform_link(
        var,
        technology="Boiler",
        carrier=["rural oil boiler", "urban decentral oil boiler"],
    )
    transform_link(
        var,
        technology="Boiler",
        carrier=[
            "rural gas boiler",
            "urban central gas boiler",
            "urban decentral gas boiler",
        ],
    )
    transform_link(
        var,
        technology="Resistive Heater",
        carrier=[
            "rural resistive heater",
            "urban decentral resistive heater",
            "urban central resistive heater",
        ],
    )
    transform_link(var, technology="Ground Heat Pump", carrier="rural ground heat pump")
    transform_link(
        var,
        technology="Air Heat Pump",
        carrier=[
            "urban decentral air heat pump",
            "rural air heat pump",
            "urban central air heat pump",
        ],
    )

    # # MHC is a side product of naphtha for industry. The oil demand of
    # # the link equals the naphtha output. There are no losses.
    # # todo: move to primary
    # var[f"{SECONDARY}|Waste|Oil|Naphtha Refining"] = _extract(
    #     SUPPLY,
    #     carrier="naphtha for industry",
    #     bus_carrier="non-sequestered HVC",
    #     component="Link",
    # )

    # solid biomass is produced by some boilers, which is wrong of course
    # but needs to be addressed nevertheless to correct balances
    process_biomass_boilers(var)

    # multi input links
    transform_link(var, technology="Methanolisation", carrier="methanolisation")
    transform_link(var, technology="Electrobiofuels", carrier="electrobiofuels")
    transform_link(var, technology="Haber-Bosch", carrier="Haber-Bosch")

    # Links that connect to buses with single loads. They are skipped in
    # IAMC variables, because their buses are only needed because of
    # separated loads at bus1.
    demand_carrier = [
        "agriculture machinery oil",
        "coal for industry",
        "gas for industry",
        "gas for industry CC",
        "industry methanol",
        "kerosene for aviation",
        "land transport oil",
        "naphtha for industry",
        "shipping methanol",
        "shipping oil",
        "solid biomass for industry",
        "solid biomass for industry CC",
        # "urban central water pits charger",
        # "urban central water pits discharger",
    ]
    remaining_supply = filter_by(SUPPLY, component="Link").drop(
        demand_carrier, level="carrier", errors="ignore"
    )
    assert remaining_supply.empty, f"{remaining_supply.index.unique('carrier')}"

    remaining_demand = filter_by(DEMAND, component="Link").drop(
        demand_carrier, level="carrier", errors="ignore"
    )
    assert remaining_demand.empty, f"{remaining_demand.index.unique('carrier')}"

    return combine_variables(var)


def collect_final_energy() -> pd.Series:
    """Extract all final energy variables from the networks."""
    var = WriteOnceDict()

    load_carrier = filter_by(DEMAND, component="Load").index.unique("carrier")
    for carrier in load_carrier:
        transform_load(var, carrier)

    assert filter_by(DEMAND, component="Load").empty
    assert filter_by(SUPPLY, component="Load").empty

    # CC Links have bus0 efficiencies < 1, i.e. they have losses
    for carrier in ("gas for industry CC", "solid biomass for industry CC"):
        bc = carrier.split(" for industry")[0]
        var[f"{SECONDARY}|Losses|{BC_ALIAS[bc]}|CC"] = (
            _extract(SUPPLY, component="Link", carrier=carrier)
            .add(_extract(DEMAND, component="Link", carrier=carrier))
            .mul(-1)
        )

    assert SUPPLY.empty, f"Supply is not empty: {SUPPLY}"
    assert DEMAND.empty, f"Demand is not empty: {SUPPLY}"

    return combine_variables(var)


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
    SUPPLY = collect_myopic_statistics(
        networks, "supply", groupby=groupby, **kwargs
    ).drop("t_co2", level="unit", errors="ignore")
    DEMAND = (
        collect_myopic_statistics(networks, "withdrawal", groupby=groupby, **kwargs)
        .drop("t_co2", level="unit", errors="ignore")
        .mul(-1)
    )
    IMPORT_FOREIGN = collect_myopic_statistics(
        networks, "trade_energy", scope=TradeTypes.FOREIGN, direction="import", **kwargs
    ).drop("t_co2", level="unit", errors="ignore")
    EXPORT_FOREIGN = collect_myopic_statistics(
        networks, "trade_energy", scope=TradeTypes.FOREIGN, direction="export", **kwargs
    ).drop("t_co2", level="unit", errors="ignore")
    IMPORT_DOMESTIC = collect_myopic_statistics(
        networks,
        "trade_energy",
        scope=TradeTypes.DOMESTIC,
        direction="import",
        **kwargs,
    ).drop("t_co2", level="unit", errors="ignore")
    EXPORT_DOMESTIC = collect_myopic_statistics(
        networks,
        "trade_energy",
        scope=TradeTypes.DOMESTIC,
        direction="export",
        **kwargs,
    ).drop("t_co2", level="unit", errors="ignore")

    # necessary for Links with multiple inputs
    LINK_BALANCE = collect_myopic_statistics(
        networks, comps="Link", statistic="energy_balance"
    ).drop(["co2", "co2 stored", "process emissions"], level=DM.BUS_CARRIER)

    # all transmission is already in trade_energy.
    transmission_carrier = [t[1] for t in get_transmission_techs(networks)]
    SUPPLY.drop(transmission_carrier, level="carrier", errors="ignore", inplace=True)
    DEMAND.drop(transmission_carrier, level="carrier", errors="ignore", inplace=True)

    # collect transformed energy system variables. Note, that the order of collection is relevant.
    #
    to_concat = [
        collect_primary_energy(),
        collect_storage_imbalances(),
        collect_losses_energy(),
        collect_secondary_energy(),
        collect_final_energy(),
        collect_system_cost(),
    ]
    df = pd.concat(to_concat)
    # energy_primary = collect_primary_energy()
    # energy_storage = collect_storage_imbalances()
    # energy_losses = collect_losses_energy()
    # energy_secondary = collect_secondary_energy()
    # energy_final = collect_final_energy()
    # system_cost = collect_system_cost()

    # df = pd.concat([energy_primary, energy_storage, energy_losses, energy_secondary, energy_final, system_cost])

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
