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
    filter_by,
    get_transmission_techs,
    insert_index_level,
    rename_aggregate,
)
from scripts._helpers import configure_logging, mock_snakemake

logger = logging.getLogger(__file__)
logger.setLevel(logging.DEBUG)

YEAR_LOC = [DM.YEAR, DM.LOCATION]
IDX = YEAR_LOC + ["unit"]
PRIMARY = "Primary Energy"
SECONDARY = "Secondary Energy"
FINAL = "Final Energy"

BC_ALIAS = {
    "AC": "AC",
    "H2": "H2",
    "NH3": "NH3",
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


class SeriesCollector(dict):
    """Prevent overwriting existing keys."""

    def __setitem__(self, key: str, value: pd.Series):
        if key in self:
            idx = self[key].index.names
            assert value.index.names == idx, (
                f"Denying to join existing index levels: {idx} with {value.index.names}"
            )
            old = self.pop(key)
            super().__setitem__(key, old.add(value, fill_value=0))
            logger.debug(f"Merged key {key} with existing Series.")
        else:
            super().__setitem__(key, value)


def _process_single_input_link(
    var: dict,
    supply: pd.Series,
    demand: pd.Series,
    bc_out: pd.Index,
    bc_in: str,
    technology: str,
) -> dict:
    demand_unit = demand.index.unique("unit").item()
    losses = supply.groupby(YEAR_LOC).sum() + demand.groupby(YEAR_LOC).sum()
    losses = insert_index_level(losses, demand_unit, "unit", pos=2).mul(-1)

    var[f"{SECONDARY}|Losses|{BC_ALIAS[bc_in]}|{technology}"] = losses[losses > 0]
    surplus = losses[losses < 0]

    for bc in bc_out:
        label = f"{SECONDARY}|{BC_ALIAS[bc]}|{BC_ALIAS[bc_in]}|{technology}"
        # supply_bc = filter_by(supply, bus_carrier=bc).groupby(IDX).sum()
        # if surplus.empty:
        var[label] = filter_by(supply, bus_carrier=bc).groupby(IDX).sum()
        # else:  # split ambient heat. 100% of demand is Heat|AC, surplus is Ambient Heat|AC
        #     assert technology in ("CHP", "Boiler", "Ground Heat Pump", "Air Heat Pump"), (
        #         f"Unknown technology with efficiencies > 1: {technology}."
        #     )
        #     supply_unit = supply_bc.index.unique("unit").item()
        #     ambient_heat = rename_aggregate(surplus, supply_unit, level="unit").mul(-1)
        #     var[label] = supply_bc.sub(ambient_heat, fill_value=0)
        #     var[f"{SECONDARY}|Ambient Heat|{BC_ALIAS[bc_in]}|{technology}"] = ambient_heat

    if not surplus.empty:
        assert technology in ("CHP", "Boiler", "Ground Heat Pump", "Air Heat Pump"), (
            f"Unknown technology with efficiencies > 1: {technology}."
        )
        heat_ambient = rename_aggregate(surplus, "MWh_th", level="unit").mul(-1)
        label_heat = f"{SECONDARY}|Heat|{BC_ALIAS[bc_in]}|{technology}"
        heat_total = var.pop(label_heat)
        var[label_heat] = heat_total.sub(heat_ambient, fill_value=0)
        var[f"{SECONDARY}|Ambient Heat|{BC_ALIAS[bc_in]}|{technology}"] = heat_ambient
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
    #     # need to subtract ambient heat amounts from primary heat output
    #     lbl = f"{SECONDARY}|Heat|{bc_in}|{technology}"
    #     full_output = var.pop(lbl)
    #     var[lbl] = full_output - var[f"{SECONDARY}|Ambient Heat|{bc_in}|{technology}"]

    # do not count ambient heat for demand + supply + losses = 0 to stay true
    # bc_out_alias = [BC_ALIAS[bc] for bc in bc_out] + ["Losses"]
    # pattern = rf"{SECONDARY}\|({'|'.join(bc_out_alias)})\|{bc_in}\|{technology}"
    # total_vars = (
    #     pd.concat({k: v for k, v in var.items() if re.match(pattern, k)})
    #     .groupby(YEAR_LOC)
    #     .sum()
    # )
    # total demand + supply + losses - ambient heat = 0
    # ambient_heat = var.get(label_ambient_heat, pd.Series())
    # if not losses_neg.empty:
    #     assert (
    #             total_vars.add(demand.groupby(YEAR_LOC).sum()).add(
    #                 losses_neg.droplevel("unit"), fill_value=0
    #             )
    #             <= 1e-5
    #     ).all(), (
    #         f"Imbalances detected for Link with bus carrier supply: {bc_in} and bus carrier out: {list(bc_out)} for technology {technology}."
    #     )
    # else:
    #     assert (total_vars.add(demand.groupby(YEAR_LOC).sum()) <= 1e-5).all()

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

    assert "|" not in technology, (
        f"Pipe operator '|' not allows in technology '{technology}' because it breaks the regex match further below."
    )

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
            # scaling takes into account that Link inputs and outputs are not equally large
            # scaling = abs(supply.sum() / demand.sum())
            supply_bc = supply * demand_share  # * scaling
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
    # positive load values are possible if industry produces a surplus
    load = load_demand.add(load_supply, fill_value=0)

    supply = _extract(SUPPLY, carrier=carrier, component="Link")
    demand = _extract(DEMAND, carrier=carrier, component="Link")

    losses = supply.groupby(YEAR_LOC).sum().add(demand.groupby(YEAR_LOC).sum())
    assert losses.abs().lt(1e-5).all(), (
        f"Supply and demand are not equal. Please check for losses and efficiencies != 1 for carrier: {carrier}.\n{losses}"
    )

    var[f"{FINAL}|{bc}|{sector}"] = load.mul(-1)

    return var


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
    balances = (
        collect_myopic_statistics(networks, comps="Link", statistic="energy_balance")
        .pipe(filter_by, carrier=carrier)
        .drop(["co2", "co2 stored", "process emissions"], level=DM.BUS_CARRIER)
    )
    # balances = filter_by(LINK_BALANCE, carrier=carrier)
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
    surplus = _bal[_bal.le(0)].mul(-1)
    var[f"{SECONDARY}|Losses|Biomass|Boiler"] = losses
    if not surplus.empty:
        # assert (
        #     _bal.sub(losses, fill_value=0).sub(surplus, fill_value=0) <= 1e-5
        # ).all()  # implicitly True
        var[f"{SECONDARY}|Ambient Heat|Biomass|Boiler"] = surplus
    # else:
    #     assert (_bal.sub(losses, fill_value=0) <= 1e-5).all()  # implicitly True
    _extract(SUPPLY, carrier=carrier, component="Link")
    _extract(DEMAND, carrier=carrier, component="Link")

    return var


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
    mapper = {"": "MWh_LHV"}  # Municipal solid waste is missing the unit
    var[f"{prefix}|Import Foreign"] = _extract(IMPORT_FOREIGN, bus_carrier=bc).rename(
        mapper, axis="index", level="unit"
    )
    var[f"{prefix}|Import Domestic"] = _extract(IMPORT_DOMESTIC, bus_carrier=bc).rename(
        mapper, axis="index", level="unit"
    )
    var[f"{prefix}|Solid"] = _extract(
        SUPPLY, bus_carrier=bc, component="Generator"
    ).rename(mapper, axis="index", level="unit")

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
    var[f"{prefix}|HVC from naphtha"] = _extract(
        SUPPLY,
        carrier="naphtha for industry",
        bus_carrier="non-sequestered HVC",
        component="Link",
    )

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

    return var


def primary_electricity(var: dict) -> dict:
    """
    Calculate the electricity generated per region.

    Parameters
    ----------
    var
        The collection to append new variables in.

    Returns
    -------
    :
        The updated variables' collection.
    """
    prefix = f"{PRIMARY}|AC"
    # var["Primary Energy|Hydro|PHS"] = _extract(
    #     SUPPLY, carrier="PHS", component="StorageUnit"
    # )
    var[f"{prefix}|Reservoir"] = _extract(
        SUPPLY, carrier="hydro", component="StorageUnit"
    )
    var[f"{prefix}|Run-of-River"] = _extract(
        SUPPLY, carrier="ror", component="Generator"
    )
    var[f"{prefix}|Wind Onshore"] = _extract(
        SUPPLY, carrier="onwind", component="Generator"
    )
    var[f"{prefix}|Wind Offshore"] = _extract(
        SUPPLY, carrier=["offwind-ac", "offwind-dc"], component="Generator"
    )
    var[f"{prefix}|Solar Utility"] = _extract(
        SUPPLY, carrier="solar", component="Generator"
    )
    var[f"{prefix}|Solar HSAT"] = _extract(
        SUPPLY, carrier="solar-hsat", component="Generator"
    )
    var[f"{prefix}|Solar Rooftop"] = _extract(
        SUPPLY, carrier="solar rooftop", component="Generator"
    )

    var[f"{prefix}|Import Domestic"] = _extract(IMPORT_DOMESTIC, bus_carrier="AC")
    var[f"{prefix}|Import Foreign"] = _extract(IMPORT_FOREIGN, bus_carrier="AC")

    return var


def primary_uranium(var: dict) -> dict:
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
    var[f"{PRIMARY}|Uranium"] = (
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
    var[f"{PRIMARY}|NH3|Import"] = _extract(SUPPLY, carrier="import NH3")

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
    prefix = f"{PRIMARY}|{BC_ALIAS.get('rural heat', 'Heat')}"

    carrier = [c for c in SUPPLY.index.unique("carrier") if "solar thermal" in c]
    var[f"{prefix}|Solar thermal"] = _extract(
        SUPPLY, carrier=carrier, component="Generator"
    )
    var[f"{prefix}|Geothermal"] = _extract(SUPPLY, carrier="geothermal heat")

    return var


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
    var = SeriesCollector()
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
    primary_electricity(var)
    primary_uranium(var)
    primary_ammonia(var)
    primary_heat(var)
    primary_methanol(var)

    assert filter_by(SUPPLY, component="Generator").empty
    assert IMPORT_DOMESTIC.empty, f"Import domestic is not empty: {IMPORT_DOMESTIC}"
    assert IMPORT_FOREIGN.empty, f"Import foreign is not empty: {IMPORT_FOREIGN}"

    return combine_variables(var)


def collect_storage_imbalances() -> pd.Series:
    """Extract all storage imbalances due to losses."""
    var = SeriesCollector()
    comps = ["Store", "StorageUnit"]

    imbalanced_techs = {
        "urban central water pits": "Water Pits",  # Storage losses
        "urban central water tanks": "Water Tank",  # Storage losses
        "coal": "Coal",  # FixMe: small unexplained imbalance accepted for now
        "PHS": "PHS",  # Pump efficiency
    }

    for carrier in filter_by(SUPPLY, component=comps).index.unique("carrier"):
        supply = filter_by(SUPPLY, component=comps, carrier=carrier)
        demand = filter_by(DEMAND, component=comps, carrier=carrier)
        balance = supply.add(demand, fill_value=0)

        if balance.sum() != 0:
            logger.warning(
                f"Store imbalances detected for carrier {carrier} with "
                f"total imbalance of {balance.groupby('year').sum()}."
            )
            # should raise KeyError if bc is not registered
            bc = balance.index.unique("bus_carrier").item()
            tech = imbalanced_techs[carrier]
            var[f"{SECONDARY}|Losses|{BC_ALIAS[bc]}|{tech}"] = balance.groupby(
                IDX
            ).sum()
        else:
            logger.debug(f"No Store imbalances detected for carrier: {carrier}.")

        SUPPLY.drop(supply.index, inplace=True)
        DEMAND.drop(demand.index, inplace=True)

    # collect water pit losses here, to combine them with storage
    # labels and prevent duplicates in the IAMC data frame.
    var[f"{SECONDARY}|Losses|Heat|Water Pits"] = (
        _extract(
            SUPPLY,
            carrier="urban central water pits discharger",
            bus_carrier="urban central heat",
            component="Link",
        )
        .add(
            _extract(
                DEMAND,
                carrier="urban central water pits charger",
                bus_carrier="urban central heat",
                component="Link",
            )
        )
        .mul(-1)
    )

    return combine_variables(var)


def collect_losses_energy() -> pd.Series:
    # DSM Links with losses that connect to buses with stores
    var = SeriesCollector()
    prefix = f"{SECONDARY}|Losses"

    var[f"{prefix}|AC|Distribution Grid"] = (
        _extract(SUPPLY, carrier="electricity distribution grid")
        + _extract(DEMAND, carrier="electricity distribution grid")
    ).mul(-1)

    var[f"{prefix}|AC|BEV charger"] = (
        _extract(SUPPLY, carrier="BEV charger", component="Link")
        .add(_extract(DEMAND, carrier="BEV charger", component="Link"))
        .mul(-1)
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

    # gas and hydrogen compressing cost energy
    var[f"{prefix}|AC|H2 Compressing"] = _extract(
        DEMAND,
        carrier=["H2 pipeline", "H2 pipeline (Kernnetz)", "H2 pipeline retrofitted"],
        component="Link",
        bus_carrier="AC",
    ).mul(-1)
    var[f"{prefix}|AC|Gas Compressing"] = _extract(
        DEMAND,
        carrier=["gas pipeline", "gas pipeline new"],
        component="Link",
        bus_carrier="AC",
    ).mul(-1)

    return combine_variables(var)


def collect_secondary_energy() -> pd.Series:
    """Extract all secondary energy variables from the networks."""
    var = SeriesCollector()

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

    # solid biomass is produced by some boilers, which is wrong
    # but needs to be addressed nevertheless to correct balances
    process_biomass_boilers(var)

    # multi input links
    transform_link(var, technology="Methanolisation", carrier="methanolisation")
    transform_link(var, technology="Electrobiofuels", carrier="electrobiofuels")
    transform_link(var, technology="Haber-Bosch", carrier="Haber-Bosch")
    # link_balance = collect_myopic_statistics(
    #     networks,
    #     "energy_balance",
    #     comps="Link",
    # )
    # from evals.utils import filter_for_carrier_connected_to, calculate_input_share
    # inputs_for_nh3 = (
    #     link_balance.filter(like="Haber-Bosch")
    #     .drop(["co2", "co2 stored"], level="bus_carrier")
    #     .pipe(filter_for_carrier_connected_to, "H2")
    #     .pipe(calculate_input_share, "NH3")
    # )
    # var[f"{SECONDARY}|NH3|H2|Haber-Bosch"] = inputs_for_nh3.xs("H2", level="bus_carrier").groupby(YEAR_LOC).sum()
    # var[f"{SECONDARY}|NH3|AC|Haber-Bosch"] = inputs_for_nh3.xs("AC", level="bus_carrier").groupby(YEAR_LOC).sum()
    #
    # var[f"{SECONDARY}|Losses|H2|Haber-Bosch"] = inputs_for_nh3.xs("AC", level="bus_carrier").groupby(YEAR_LOC).sum()
    # var[f"{SECONDARY}|Losses|AC|Haber-Bosch"] = inputs_for_nh3.xs("AC", level="bus_carrier").groupby(YEAR_LOC).sum()

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
    var = SeriesCollector()

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

    var[f"{SECONDARY}|Losses|Heat|Vent"] = _extract(
        DEMAND,
        component="Generator",
        carrier=[
            "urban central heat vent",
            "rural heat vent",
            "urban decentral heat vent",
        ],
    ).mul(-1)

    assert SUPPLY.empty, f"Supply is not empty: {SUPPLY}"
    assert DEMAND.empty, f"Demand is not empty: {DEMAND}"

    for bus_carrier in EXPORT_DOMESTIC.index.unique("bus_carrier"):
        prefix = f"{FINAL}|{BC_ALIAS[bus_carrier]}"
        var[f"{prefix}|Export Domestic"] = _extract(
            EXPORT_DOMESTIC, bus_carrier=bus_carrier
        )
        var[f"{prefix}|Export Foreign"] = _extract(
            EXPORT_FOREIGN, bus_carrier=bus_carrier
        )

    assert EXPORT_DOMESTIC.empty, f"Export domestic is not empty: {EXPORT_DOMESTIC}"
    assert EXPORT_FOREIGN.empty, f"Export foreign is not empty: {EXPORT_FOREIGN}"

    return combine_variables(var)


if __name__ == "__main__":
    if "snakemake" not in globals():
        snakemake = mock_snakemake(
            "export_iamc_variables",
            run="KN2045_Mix",
        )
    configure_logging(snakemake)

    networks = read_networks(snakemake.input.networks)

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
    EXPORT_FOREIGN = (
        collect_myopic_statistics(
            networks,
            "trade_energy",
            scope=TradeTypes.FOREIGN,
            direction="export",
            **kwargs,
        )
        .drop("t_co2", level="unit", errors="ignore")
        .mul(-1)
    )
    IMPORT_DOMESTIC = collect_myopic_statistics(
        networks,
        "trade_energy",
        scope=TradeTypes.DOMESTIC,
        direction="import",
        **kwargs,
    ).drop("t_co2", level="unit", errors="ignore")
    EXPORT_DOMESTIC = (
        collect_myopic_statistics(
            networks,
            "trade_energy",
            scope=TradeTypes.DOMESTIC,
            direction="export",
            **kwargs,
        )
        .drop("t_co2", level="unit", errors="ignore")
        .mul(-1)
    )

    # all transmission is already in trade_energy. The bus_carrier must be
    # considered in the filter to prevent dropping compression costs f
    transmission_bus_carrier = {
        "AC": "AC",
        "CO2 pipeline": "co2",
        "DC": "AC",
        "H2 pipeline": "H2",
        "H2 pipeline (Kernnetz)": "H2",
        "H2 pipeline retrofitted": "H2",
        "gas pipeline": "gas",
        "gas pipeline new": "gas",
        "municipal solid waste transport": "municipal solid waste",
        "solid biomass transport": "solid biomass",
    }
    # transmission_carrier = [t[1] for t in get_transmission_techs(networks)]
    for component, carrier in get_transmission_techs(networks):
        bus_carrier = transmission_bus_carrier[carrier]
        SUPPLY.drop(
            filter_by(
                SUPPLY, component=component, carrier=carrier, bus_carrier=bus_carrier
            ).index,
            inplace=True,
        )
        DEMAND.drop(
            filter_by(
                DEMAND, component=component, carrier=carrier, bus_carrier=bus_carrier
            ).index,
            inplace=True,
        )
    # SUPPLY.drop(transmission_carrier, level="carrier", errors="ignore", inplace=True)
    # DEMAND.drop(transmission_carrier, level="carrier", errors="ignore", inplace=True)

    # collect transformed energy system variables. Note, that the order of
    # collection is relevant for assertions statements.
    to_concat = [
        collect_primary_energy(),
        collect_storage_imbalances(),
        collect_losses_energy(),
        collect_secondary_energy(),
        collect_final_energy(),
        collect_system_cost(),
    ]
    df = pd.concat(to_concat)

    df = insert_index_level(df, "PyPSA-AT", "model")
    df = insert_index_level(df, snakemake.wildcards.run, "scenario")
    df.index = df.index.rename({"location": "region"})  # comply with IAMC data model

    if df.index.has_duplicates:
        print(df[df.index.duplicated()])

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
