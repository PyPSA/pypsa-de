"""
Export variables in IAMC data model for all regions.

IAMC variable naming convention:
Category|Subcategory|Specification

naming conventions:
* Primary Energy|<bus_carrier>|Subcategory
* Secondary Energy|<output_bus_carrier>|<input_bus_carrier>|Subcategory
  For example, Secondary Energy|H2|AC|Electrolysis should be read as
  Hydrogen supply from Electricity (AC + low voltage) using Electrolysis
* Final Energy|<bus_carrier>|Subcategory

Notes
-----
https://docs.ece.iiasa.ac.at/iamc.html
https://pyam-iamc.readthedocs.io/en/stable/
"""

import logging
import re

import pandas as pd
from pyam import IamDataFrame

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
TRANS_IN = "Transformation Input"
TRANS_OUT = "Transformation Output"
TRANS_BYPASS = "Transformation Bypass"

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
    "unsustainable bioliquids": "Oil",
    # load buses
    "naphtha for industry": "Oil",
    "agriculture machinery oil": "Oil",
    "coal for industry": "Coal",
    "gas for industry": "Gas",
    "kerosene for aviation": "Oil",
    "land transport oil": "Oil",
    "industry methanol": "Methanol",
    "shipping methanol": "Methanol",
    "shipping oil": "Oil",
    "solid biomass for industry": "Biomass",
    "EV battery": "AC",
    # Store
    "urban central water pits": "Heat",
    "urban central water tanks": "Heat",
    "urban decentral water tanks": "Heat",
    "rural water tanks": "Heat"
}


class SeriesCollector(dict):
    """Prevent overwriting existing keys."""

    def __setitem__(self, key: str, value: pd.Series):
        if key in self:
            idx = self[key].index.names
            assert value.index.names == idx, (
                f"Denying to join mismatching index: {idx} with {value.index.names}"
            )
            old = self.pop(key)
            super().__setitem__(key, old.add(value, fill_value=0))
            logger.debug(f"Merged key {key} with existing Series.")
        else:
            super().__setitem__(key, value)


def _process_single_input_link(
    supply: pd.Series,
    demand: pd.Series,
    bc_out: pd.Index,
    bc_in: str,
    technology: str,
):
    demand_unit = demand.index.unique("unit").item()
    losses = supply.groupby(YEAR_LOC).sum() + demand.groupby(YEAR_LOC).sum()
    losses = insert_index_level(losses, demand_unit, "unit", pos=2).mul(-1)

    var[f"{SECONDARY}|Demand|{BC_ALIAS[bc_in]}|{technology}"] = (
        demand.groupby(IDX).sum().mul(-1)
    )
    var[f"{SECONDARY}|Losses|{BC_ALIAS[bc_in]}|{technology}"] = losses[losses > 0]
    surplus = losses[losses < 0]

    for bc in bc_out:
        label = f"{SECONDARY}|{BC_ALIAS[bc]}|{BC_ALIAS[bc_in]}|{technology}"
        var[label] = filter_by(supply, bus_carrier=bc).groupby(IDX).sum()

    if not surplus.empty:
        assert technology in ("CHP", "CHP CC", "Boiler", "Ground Heat Pump", "Air Heat Pump"), (
            f"Unexpected technology with efficiencies > 1: {technology}."
        )
        heat_ambient = rename_aggregate(surplus, "MWh_th", level="unit").mul(-1)
        label_heat = f"{SECONDARY}|Heat|{BC_ALIAS[bc_in]}|{technology}"
        heat_total = var.pop(label_heat)
        var[label_heat] = heat_total.sub(heat_ambient, fill_value=0)
        var[f"{SECONDARY}|Ambient Heat|{BC_ALIAS[bc_in]}|{technology}"] = heat_ambient


def aggregate_variables(label: str, pattern: str):
    """Aggregate a subset of variables."""
    assert label not in var, (
        f"Adding to existing keys causes data duplication. key={label}"
    )
    to_sum = {k: v for k, v in var.items() if re.match(pattern, k)}
    if to_sum:
        # variables may have different units. Overwriting units to
        # yield one row per year and location to use in sankey diagrams
        year_sum = pd.concat(to_sum).groupby(YEAR_LOC).sum()
        var[label] = insert_index_level(year_sum, "MWh", "unit", pos=2)
    else:
        logger.debug(
            f"No matches for label {label} and pattern {pattern}. Skipping aggregaion."
        )


def merge_variables(collection: SeriesCollector | dict) -> pd.Series | pd.DataFrame:
    """
    Combine variables into a single data series.

    Parameters
    ----------
    collection

    Returns
    -------
    :
    """
    to_concat = {k: v for k, v in collection.items() if not v.empty}

    if len(to_concat) == 0:
        return pd.Series()

    ds = pd.concat(to_concat)
    ds.index = ds.index.rename({None: "Variable"})

    return ds


def transform_link(carrier: str | list, technology: str) -> None:
    """
    Transform a Link component into supply and transformation losses.

    The Link demand is equal to supply + losses and not included in
    the output to avoid redundant data. Losses have positive signs
    and the demand bus unit.

    Parameters
    ----------
    carrier
    technology

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
        return

    bc_in = demand.index.unique("bus_carrier")
    bc_out = supply.index.unique("bus_carrier")

    if len(bc_in) > 1:
        for bus_carrier_demand in bc_in:
            demand_bc = filter_by(demand, bus_carrier=bus_carrier_demand)
            demand_share = demand_bc.sum() / demand.sum()
            # scaling takes into account that Link inputs and outputs are not equally large
            # scaling = abs(supply.sum() / demand.sum())
            supply_bc = supply * demand_share  # * scaling
            _process_single_input_link(
                supply_bc,
                demand_bc,
                bc_out,
                bus_carrier_demand,
                technology,
            )
    else:
        _process_single_input_link(supply, demand, bc_out, bc_in.item(), technology)

    # remove from global statistic to prevent double counting
    SUPPLY.drop(supply.index, inplace=True)
    # adding demand to IAMC is redundant, because supply + losses == demand,
    # and we the demand bus_carrier is known from the variable name.
    DEMAND.drop(demand.index, inplace=True)


def transform_load(carrier: str):
    """

    Parameters
    ----------
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
    elif carrier.startswith("agriculture"):  # must come before heat
        sector = "Agriculture"
    elif carrier.endswith("heat"):
        sector = "HH & Services"
    elif carrier == "electricity":
        sector = "Base Load"  # todo: sector load split
        # Base Load contains a mix Transport, Industry, Households and service
    elif "industry" in carrier:
        sector = "Industry"
    else:
        raise ValueError(f"Unknown sector for Load carrier: {carrier}.")

    bc = {
        BC_ALIAS[bc]
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


def _extract(ds: pd.Series, **filter_kwargs) -> pd.Series:
    """Extract and group filter results."""
    results = filter_by(ds, **filter_kwargs)
    ds.drop(results.index, inplace=True)
    return results.groupby(IDX).sum()


def drop_transmission_technologies():
    # all transmission is already in trade_energy. The bus_carrier must be
    # considered in the filter to prevent dropping compression costs
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


def collect_regional_nh3_loads():
    # use regional surplus as regional Load
    bc = "NH3"

    nh3_regional_supply = merge_variables(
        {
            k: v
            for k, v in var.items()
            if re.match(rf"^Secondary Energy\|{BC_ALIAS[bc]}", k)
        }
    )
    nh3_eu_demand = DEMAND.filter(like=BC_ALIAS[bc])
    nh3_eu_import = merge_variables(
        {
            k: v
            for k, v in var.items()
            if re.match(rf"^Primary Energy\|{BC_ALIAS[bc]}$", k)
        }
    )
    imbalances_iamc = (
        nh3_regional_supply.groupby("year")
        .sum()
        .add(nh3_eu_import.groupby("year").sum(), fill_value=0)
        .add(nh3_eu_demand.groupby("year").sum(), fill_value=0)
    )
    imbalances_bus = (
        collect_myopic_statistics(
            networks,
            "energy_balance",
            groupby=["location", "carrier"],
            bus_carrier=bc,
            aggregate_components=None,
        )
        .groupby("year")
        .sum()
    )

    # check that imbalances are equal to imbalances in the network
    pd.testing.assert_series_equal(imbalances_iamc, imbalances_bus, check_names=False)
    if not imbalances_bus.empty:
        logger.warning(f"Imbalances detected for bus carrier {bc}:\n{imbalances_bus}.")
    var[f"{FINAL}|{BC_ALIAS[bc]}|Agriculture"] = nh3_regional_supply.groupby(IDX).sum()
    _extract(DEMAND, component="Load", carrier=bc, bus_carrier=bc)


def process_biomass_boilers() -> None:
    """
    Special processing biomass boilers that have solid biomass supply.
    """
    carrier = ["rural biomass boiler", "urban decentral biomass boiler"]
    if len(filter_by(DEMAND, carrier=carrier).index.unique("bus_carrier")) <= 1:
        transform_link(carrier, technology="Boiler")
        return

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
        var[f"{SECONDARY}|Ambient Heat|Biomass|Boiler"] = surplus

    _extract(SUPPLY, carrier=carrier, component="Link")
    _extract(DEMAND, carrier=carrier, component="Link")


def primary_oil():
    """
    Calculate the amounts of oil entering a region.

    Returns
    -------
    :

    Notes
    -----
    proper tests must make sure that no oil amounts
    in the network are skipped.
    """
    prefix = f"{PRIMARY}|Oil"
    bc = "oil"

    # assuming that all oil production is consumed locally.
    # Let's not filter_by components, to capture anything but Stores.
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
    regional_deficit = consumption.add(production, fill_value=0).clip(upper=0).mul(-1)
    regional_surplus = consumption.add(production, fill_value=0).clip(lower=0)

    var[f"{prefix}|Import"] = regional_deficit
    var[f"{FINAL}|Oil|Export"] = regional_surplus

    aggregate_variables(prefix, pattern=rf"^{prefix.replace('|', r'\|')}")

    # remove EU imports and oil refining
    _extract(SUPPLY, carrier="import oil")
    _extract(SUPPLY, carrier="oil primary")
    _extract(SUPPLY, carrier="oil refining")
    _extract(DEMAND, carrier="oil refining")

    # unsustainable bioliquids have regional bus generators.
    # The "unsustainable bioliquids" Link forwards all generated energy
    # to the oil bus.
    _extract(SUPPLY, carrier="unsustainable bioliquids", component="Generator")


def primary_gas():
    """
    Calculate the amount of gas entering a region.

    Returns
    -------
    :
    """
    bc = "gas"
    prefix = f"{PRIMARY}|{BC_ALIAS[bc]}"
    var[f"{prefix}|Import Foreign"] = _extract(IMPORT_FOREIGN, bus_carrier=bc)
    var[f"{prefix}|Import Domestic"] = _extract(IMPORT_DOMESTIC, bus_carrier=bc)

    var[f"{prefix}|Global Import LNG"] = _extract(SUPPLY, bus_carrier=bc, component="Generator", carrier="lng gas")
    var[f"{prefix}|Global Import Pipeline"] = _extract(SUPPLY, bus_carrier=bc, component="Generator", carrier="pipeline gas")
    var[f"{prefix}|Domestic Production"] = _extract(SUPPLY, bus_carrier=bc, component="Generator", carrier="production gas")
    var[f"{prefix}|Green Global Import"] = _extract(SUPPLY, carrier="import gas", bus_carrier=bc, component="Link")

    var[f"{prefix}|Biogas"] = _extract(SUPPLY, carrier="biogas to gas", bus_carrier=bc)
    var[f"{prefix}|Biogas CC"] = _extract(
        SUPPLY, carrier="biogas to gas CC", bus_carrier=bc
    )

    aggregate_variables(prefix, pattern=rf"^{prefix.replace('|', r'\|')}")

    # drop biogas withdrawal from biogas processing
    _extract(
        DEMAND, carrier=["biogas to gas", "biogas to gas CC"], bus_carrier="biogas"
    )


def primary_waste():
    """
    Collect the amounts of municipal solid wate and HVC entering a region.

    Returns
    -------
    :
    """
    bc = "municipal solid waste"
    prefix = f"{PRIMARY}|{BC_ALIAS[bc]}"
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

    # HVC is a side product of naphtha for industry. The oil demand of
    # the link equals the naphtha output. There are no losses.
    var[f"{prefix}|HVC from naphtha"] = _extract(
        SUPPLY,
        carrier="naphtha for industry",
        bus_carrier="non-sequestered HVC",
        component="Link",
    )

    aggregate_variables(prefix, pattern=rf"^{prefix.replace('|', r'\|')}")

    # municipal solid waste is only used to transform "municipal solid waste" to
    # "non-sequestered HVC" and to track CO2. Same as Biogas processing.
    _extract(
        SUPPLY,
        carrier="municipal solid waste",
        bus_carrier="non-sequestered HVC",
        component="Link",
    )
    _extract(DEMAND, carrier="municipal solid waste", bus_carrier=bc, component="Link")


def primary_coal():
    """
    Calculate the amounts of coal consumed in a region.

    Coal is not produced by any Link, therefore it's safe to assume
    all Link withdrawal is imported fossil coal or lignite.

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

    aggregate_variables(prefix, pattern=rf"^{prefix.replace('|', r'\|')}")

    # remove EU coal generators from the to-do list
    coal_generators = filter_by(  # todo: use _extract() shorthand
        SUPPLY, bus_carrier=["coal", "lignite"], component="Generator"
    )
    SUPPLY.drop(coal_generators.index, inplace=True)


def primary_hydrogen():
    """
    Calculate the amounts of hydrogen imported into a region.

    There are global import of Hydrogen, found in `Generator`
    components, and various types of H2 pipelines, that bring
    Hydrogen into regions.

    Returns
    -------
    :
        The updated variables' collection.
    """
    bc = "H2"
    prefix = f"{PRIMARY}|{BC_ALIAS[bc]}"
    var[f"{prefix}|Import Foreign"] = _extract(IMPORT_FOREIGN, bus_carrier=bc)
    var[f"{prefix}|Import Domestic"] = _extract(IMPORT_DOMESTIC, bus_carrier=bc)
    var[f"{prefix}|Green Import Global"] = _extract(
        SUPPLY, carrier="import H2", bus_carrier=bc, component="Generator"
    )

    aggregate_variables(prefix, pattern=rf"^{prefix.replace('|', r'\|')}")


def primary_biomass():
    """
    Calculate the amounts of biomass generated in a region.

    Returns
    -------
    :
        The updated variables' collection.
    """
    bc = "solid biomass"
    prefix = f"{PRIMARY}|{BC_ALIAS[bc]}"
    var[f"{prefix}|Import Foreign"] = _extract(IMPORT_FOREIGN, bus_carrier=bc)
    var[f"{prefix}|Import Domestic"] = _extract(IMPORT_DOMESTIC, bus_carrier=bc)

    var[f"{prefix}|Solid"] = _extract(SUPPLY, bus_carrier=bc, component="Generator")

    aggregate_variables(prefix, pattern=rf"^{prefix.replace('|', r'\|')}")

    # biogas is a separate bus carrier group
    var[f"{PRIMARY}|Biogas"] = (
        _extract(  # todo: Biogas is simplified, either include in BC_ALIAS or drop
            SUPPLY, bus_carrier="biogas", component="Generator"
        )
    )


def primary_electricity():
    """
    Calculate the electricity generated per region.

    Returns
    -------
    :
    """
    prefix = f"{PRIMARY}|AC"

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

    aggregate_variables(prefix, pattern=rf"^{prefix.replace('|', r'\|')}")


def primary_uranium():
    """
    Calculate the uranium demand for nuclear power plants per region.

    Returns
    -------
    :
        The updated variables' collection.
    """
    bc = "uranium"
    prefix = f"{PRIMARY}|{BC_ALIAS[bc]}"
    var[f"{PRIMARY}|Uranium|Import"] = (
        filter_by(DEMAND, carrier="nuclear", bus_carrier="uranium", component="Link")
        .groupby(IDX)
        .sum()
        .mul(-1)
    )

    # todo: assert var[f"{PRIMARY}|Uranium"].sum() == EU Generator
    _extract(SUPPLY, bus_carrier="uranium", component="Generator", location="EU")

    aggregate_variables(prefix, pattern=rf"^{prefix.replace('|', r'\|')}")


def primary_ammonia():
    """
    Calculate the ammonium imported per region.

    Returns
    -------
    :
    """
    # there is no regional ammonium demand
    bc = "NH3"
    prefix = f"{PRIMARY}|{BC_ALIAS[bc]}"
    # todo: needed?
    var[f"{prefix}|Green Global Import"] = _extract(SUPPLY, carrier="import NH3")

    aggregate_variables(prefix, pattern=rf"^{prefix.replace('|', r'\|')}")


def primary_methanol():
    """
    Calculate methanol imported per region.

    Returns
    -------
    :
    """
    bc = "methanol"
    prefix = f"{PRIMARY}|{BC_ALIAS[bc]}"
    regional_demand = (
        filter_by(DEMAND, bus_carrier=bc, component="Link").groupby(IDX).sum()
    )
    regional_production = (
        filter_by(SUPPLY, bus_carrier=bc, component="Link")
        .drop("import methanol", level="carrier", errors="ignore")
        .groupby(IDX)
        .sum()
    )

    deficit = regional_demand.add(regional_production, fill_value=0)
    var[f"{prefix}|Green Global Import"] = deficit.mul(-1)

    aggregate_variables(prefix, pattern=rf"^{prefix.replace('|', r'\|')}")

    _extract(SUPPLY, carrier="import methanol", location="EU")


def primary_heat():
    """
    Calculate heat generation and enthalpy of evaporation.

    Returns
    -------
    :
    """
    prefix = f"{PRIMARY}|{BC_ALIAS.get('rural heat', 'Heat')}"

    carrier = [c for c in SUPPLY.index.unique("carrier") if "solar thermal" in c]
    var[f"{prefix}|Solar thermal"] = _extract(
        SUPPLY, carrier=carrier, component="Generator"
    )
    var[f"{prefix}|Geothermal"] = _extract(SUPPLY, carrier="geothermal heat")

    aggregate_variables(prefix, pattern=rf"^{prefix.replace('|', r'\|')}")


def collect_system_cost() -> pd.Series:
    """
    Extract total energy system cost per region.

    Returns
    -------
    :
    """
    # Nodal OPEX and nodal CAPEX in billion EUR2020
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

    var["System Costs|OPEX"] = myopic_opex.groupby(IDX).sum()
    var["System Costs|CAPEX"] = myopic_capex.groupby(IDX).sum()
    var["System Costs"] = var["System Costs|CAPEX"] + var["System Costs|OPEX"]

    return merge_variables(var)


def collect_primary_energy():
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
    """
    primary_gas()
    primary_oil()
    primary_hydrogen()
    primary_waste()
    primary_coal()
    primary_biomass()
    primary_electricity()
    primary_uranium()
    primary_ammonia()
    primary_heat()
    primary_methanol()

    assert filter_by(SUPPLY, component="Generator").empty, (
        f"Generators are not empty: {filter_by(SUPPLY, component='Generator')}"
    )
    assert IMPORT_DOMESTIC.empty, f"Import domestic is not empty: {IMPORT_DOMESTIC}"
    assert IMPORT_FOREIGN.empty, f"Import foreign is not empty: {IMPORT_FOREIGN}"


def collect_storage_imbalances():
    """Extract all storage imbalances due to losses."""
    comps = ["Store", "StorageUnit"]

    imbalanced_techs = {
        # Storage losses:
        "urban central water pits": "Water Pits",
        "urban central water tanks": "Water Tank",
        "urban decentral water tanks": "Water Tank",
        "rural water tanks": "Water Tank",
        "coal": "Coal",  # FixMe: small unexplained imbalance accepted for now
        "PHS": "PHS",  # Pump efficiency
        "non-sequestered HVC": "Waste",
    }

    for carrier in filter_by(SUPPLY, component=comps).index.unique("carrier"):
        supply = filter_by(SUPPLY, component=comps, carrier=carrier)
        demand = filter_by(DEMAND, component=comps, carrier=carrier)
        balance = supply.add(demand, fill_value=0).mul(-1)

        if balance.sum() != 0:
            logger.warning(
                f"Store imbalances detected for carrier {carrier} with "
                f"total imbalance of {balance.groupby('year').sum()}."
            )
            bc = balance.index.unique("bus_carrier").item()
            label = f"{SECONDARY}|Losses|{BC_ALIAS[bc]}|{imbalanced_techs[carrier]}"
            var[label] = balance.groupby(IDX).sum()
        else:
            logger.debug(f"No Store imbalances detected for carrier: {carrier}.")

        SUPPLY.drop(supply.index, inplace=True)
        DEMAND.drop(demand.index, inplace=True)




def collect_storage_charger_discharger_pairs():
    # Assuming, that Links used to supply to storages have efficiencies of 1.0
    # i.e. they do not have losses themselves and the supply/demand balance
    # from Store components contain all standing losses.
    # drop the supply/demand at the other bus side of (dis)charger links
    storage_systems = ("rural water tanks", "urban central water tanks", "urban decentral water tanks", "urban central water pits", "battery", "home battery")

    for storage_system in storage_systems:
        charger_losses = _extract(SUPPLY, carrier=f"{storage_system} charger").add(_extract(
            DEMAND,
            carrier=f"{storage_system} charger"
        ))
        assert charger_losses.abs().le(1.5).all(), f"Charger Losses detected for carrier: {storage_system}"
        discharger_losses = _extract(
            SUPPLY,
            carrier=f"{storage_system} discharger"
        ).add(_extract(
            DEMAND,
            carrier=f"{storage_system} discharger"
        ))
        assert discharger_losses.abs().le(1.5).all(), f"Storage system imbalances detected for carrier: {storage_system}"


def collect_losses_energy():
    prefix = f"{SECONDARY}|Losses"

    var[f"{prefix}|AC|Distribution Grid"] = (
        _extract(SUPPLY, carrier="electricity distribution grid")
        .add(_extract(DEMAND, carrier="electricity distribution grid")
    ).mul(-1))

    var[f"{prefix}|AC|BEV charger"] = (
        _extract(SUPPLY, carrier="BEV charger", component="Link")
        .add(_extract(DEMAND, carrier="BEV charger", component="Link"))
        .mul(-1)
    )
    var[f"{prefix}|AC|V2G"] = (
        _extract(SUPPLY, carrier="V2G", component="Link")
        .add(_extract(DEMAND, carrier="V2G", component="Link"))
        .mul(-1)
    )
    var[f"{prefix}|AC|Large Batteries"] = (
        _extract(SUPPLY, carrier="home battery discharger", component="Link")
        .add(_extract(DEMAND, carrier="battery charger", component="Link"))
        .mul(-1)
    )
    var[f"{prefix}|AC|Home Batteries"] = (
        _extract(SUPPLY, carrier="battery discharger", component="Link")
        .add(_extract(DEMAND, carrier="battery charger", component="Link"))
        .mul(-1)
    )

    # DAC has no outputs but CO2, which is ignored in energy flows
    var[f"{SECONDARY}|Demand|AC|DAC"] = _extract(
        DEMAND, carrier="DAC", bus_carrier="AC"
    ).mul(-1)
    var[f"{SECONDARY}|Demand|Heat|DAC"] = _extract(
        DEMAND,
        carrier="DAC",
        bus_carrier=["rural heat", "urban decentral heat", "urban central heat"],
    ).mul(-1)
    var[f"{SECONDARY}|Demand|Waste|HVC to air"] = _extract(
        DEMAND,
        carrier="HVC to air",
        component="Link",
        bus_carrier="non-sequestered HVC",
    ).mul(-1)

    # gas and hydrogen compressing cost energy
    var[f"{SECONDARY}|Demand|AC|H2 Compressing"] = _extract(
        DEMAND,
        carrier=["H2 pipeline", "H2 pipeline (Kernnetz)", "H2 pipeline retrofitted"],
        component="Link",
        bus_carrier="AC",
    ).mul(-1)
    var[f"{SECONDARY}|Demand|AC|Gas Compressing"] = _extract(
        DEMAND,
        carrier=["gas pipeline", "gas pipeline new"],
        component="Link",
        bus_carrier="AC",
    ).mul(-1)


def collect_secondary_energy():
    """Extract all secondary energy variables from the networks."""

    transform_link(technology="CHP", carrier="urban central gas CHP")
    transform_link(technology="CHP", carrier="urban central oil CHP")
    transform_link(technology="CHP", carrier="urban central coal CHP")
    transform_link(technology="CHP", carrier="urban central lignite CHP")
    transform_link(
        technology="CHP",
        carrier=["urban central H2 CHP", "urban central H2 retrofit CHP"],
    )
    transform_link(technology="CHP", carrier="urban central solid biomass CHP")
    transform_link(technology="CHP", carrier="waste CHP")

    transform_link(technology="CHP CC", carrier="waste CHP CC")
    transform_link(technology="CHP CC", carrier="urban central gas CHP CC")
    transform_link(technology="CHP CC", carrier="urban central solid biomass CHP CC")

    transform_link(technology="Powerplant", carrier=["CCGT", "OCGT"])
    transform_link(technology="Powerplant", carrier="H2 OCGT")
    transform_link(technology="Powerplant", carrier="H2 Fuel Cell")
    transform_link(technology="Powerplant", carrier=["OCGT methanol", "CCGT methanol", "CCGT methanol CC", "allam methanol"])
    transform_link(technology="Powerplant", carrier="coal")
    transform_link(technology="Powerplant", carrier="oil")
    transform_link(technology="Powerplant", carrier="lignite")
    transform_link(technology="Powerplant", carrier="solid biomass")
    transform_link(technology="Powerplant", carrier="nuclear")

    transform_link(technology="BioSNG", carrier="BioSNG")
    transform_link(technology="BioSNG CC", carrier="BioSNG CC")
    transform_link(technology="Sabatier", carrier="Sabatier")

    transform_link(technology="Electrolysis", carrier="H2 Electrolysis")
    transform_link(technology="SMR", carrier="SMR")
    transform_link(technology="SMR CC", carrier="SMR CC")
    transform_link(technology="Steam Reforming", carrier="Methanol steam reforming")
    transform_link(
        technology="Steam Reforming CC", carrier="Methanol steam reforming CC"
    )

    transform_link(technology="Ammonia2Hydrogen", carrier="ammonia cracker")
    transform_link(technology="Biomass2Hydrogen", carrier="solid biomass to hydrogen")
    transform_link(technology="Biomass2Liquids", carrier="biomass to liquid")
    transform_link(technology="Biomass2Liquids CC", carrier="biomass to liquid CC")
    transform_link(technology="Biomass2Methanol", carrier="biomass-to-methanol")
    transform_link(technology="Biomass2Methanol CC", carrier="biomass-to-methanol CC")
    transform_link(technology="Fischer-Tropsch", carrier="Fischer-Tropsch")
    transform_link(technology="Methanol2Oil", carrier="methanol-to-kerosene")
    transform_link(
        technology="Unsustainable Bioliquids", carrier="unsustainable bioliquids"
    )

    transform_link(
        technology="Boiler",
        carrier=["rural oil boiler", "urban decentral oil boiler"],
    )
    transform_link(
        technology="Boiler",
        carrier=[
            "rural gas boiler",
            "urban central gas boiler",
            "urban decentral gas boiler",
        ],
    )
    transform_link(
        technology="Resistive Heater",
        carrier=[
            "rural resistive heater",
            "urban decentral resistive heater",
            "urban central resistive heater",
        ],
    )
    transform_link(technology="Ground Heat Pump", carrier="rural ground heat pump")
    transform_link(
        technology="Air Heat Pump",
        carrier=[
            "urban decentral air heat pump",
            "rural air heat pump",
            "urban central air heat pump",
            "urban central ptes heat pump",
        ],
    )

    # solid biomass is produced by some boilers, which is wrong
    # but needs to be addressed nevertheless to correct balances
    process_biomass_boilers()

    # multi input links
    transform_link(technology="Methanolisation", carrier="methanolisation")
    transform_link(technology="Electrobiofuels", carrier="electrobiofuels")
    transform_link(technology="Haber-Bosch", carrier="Haber-Bosch")

    # Links that connect to buses with single loads. They are skipped in
    # IAMC variables, because their buses are only needed to track different
    # kinds of Loads and carbon.
    demand_carrier = [
        "agriculture machinery oil",
        "coal for industry",
        "gas for industry",
        "gas for industry CC",
        "industry methanol",
        "land transport oil",
        "naphtha for industry",
        "solid biomass for industry",
        "solid biomass for industry CC",
        "shipping methanol",
        "shipping oil",
        "kerosene for aviation",
    ]
    remaining_supply = filter_by(SUPPLY, component="Link").drop(
        demand_carrier, level="carrier", errors="ignore"
    )
    assert remaining_supply.empty, f"{remaining_supply.index.unique('carrier')}"

    remaining_demand = filter_by(DEMAND, component="Link").drop(
        demand_carrier, level="carrier", errors="ignore"
    )
    assert remaining_demand.empty, f"{remaining_demand.index.unique('carrier')}"


def collect_final_energy():
    """Extract all final energy variables from the networks."""

    load_carrier = filter_by(DEMAND, component="Load").index.unique("carrier")

    # NH3 has Loads on EU bus and we need regional demands
    if "NH3" in load_carrier:
        collect_regional_nh3_loads()
        load_carrier = load_carrier.drop("NH3")

    for carrier in load_carrier:
        transform_load(carrier)

    assert filter_by(DEMAND, component="Load").empty, (
        f"Missing demand from Loads detected: {filter_by(DEMAND, component='Load')}"
    )
    assert filter_by(SUPPLY, component="Load").empty, (
        f"Missing supply from Loads detected: {filter_by(SUPPLY, component='Load')}"
    )

    # CC Links have bus0 efficiencies < 1, i.e. they have losses
    for carrier in ("gas for industry CC", "solid biomass for industry CC"):
        bc = carrier.split(" for industry")[0]
        var[f"{SECONDARY}|Losses|{BC_ALIAS[bc]}|Industry CC"] = (
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


def calculate_sankey_totals():
    """Calculate energy total inputs and outputs for sankey diagrams."""
    # negative lookahead regex to exclude Ambient by default
    exclude = "(?!.*Ambient Heat)"

    for bc in sorted(set(BC_ALIAS.values())):
        aggregate_variables(f"{FINAL}|{bc}", rf"^{FINAL}\|{bc}")
        aggregate_variables(f"{TRANS_OUT}|{bc}", rf"^{SECONDARY}\|{bc}")

        if bc == "Heat":
            exclude = ""
        # find keys that start with 'Secondary Energy|', are not 'Ambient Heat',
        # continue with any alphanumerics, and continues with the bus_carrier
        aggregate_variables(
            f"{TRANS_IN}|{bc}", rf"^{SECONDARY}\|{exclude}[a-zA-Z0-9\s]*\|{bc}"
        )

        # bypass amounts connect primary with final energy
        transformation_out = var.get(f"{TRANS_OUT}|{bc}", pd.Series())
        final_demand = var.get(f"{FINAL}|{bc}", pd.Series())
        if not final_demand.empty and not transformation_out.empty:
            bypass = final_demand.sub(transformation_out, fill_value=0).clip(lower=0)
        elif not final_demand.empty and transformation_out.empty:
            bypass = final_demand  # all comes from primary
        else:  # all used in transformation
            bypass = pd.Series()

        if not bypass.empty:
            var[f"{TRANS_BYPASS}|{bc}"] = bypass


if __name__ == "__main__":
    if "snakemake" not in globals():
        snakemake = mock_snakemake(
            "export_iamc_variables",
            run="KN2045_Mix",
            prefix="test-sector-myopic-at10",
            config="config/test/config.at10.yaml"
        )
    configure_logging(snakemake)

    # networks = read_networks(snakemake.input.networks)
    networks = read_networks("/IdeaProjects/pypsa-at/results/test-sector-myopic-at10/KN2045_Mix")

    groupby = ["location", "carrier", "bus_carrier", "unit"]
    kwargs = {
        "aggregate_components": False,
        "drop_zeros": False,
        "drop_unit": False,
    }
    # calculate all statistics and process them to IAMC data model. The idea is to
    # calculate everything once and remove rows from the global statistic. This way
    # we make sure that nothing is counted twice or forgotten.
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

    # Remove transmission technologies from SUPPLY/DEMAND because
    # they are tracken in separate IMPORT/EXPORT statistics
    drop_transmission_technologies()

    # collect transformed energy system variables. Note, that the order of
    # collection is relevant for assertions statements.
    var = SeriesCollector()

    collect_primary_energy()
    collect_storage_imbalances()
    collect_storage_charger_discharger_pairs()
    collect_losses_energy()
    collect_secondary_energy()
    collect_final_energy()
    collect_system_cost()

    calculate_sankey_totals()

    df = merge_variables(var)

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
            "Quality Assessment": "demo",
            "Release for publication": "no",
        }
    ).to_frame("value")
    meta.index.name = "key"
    with pd.ExcelWriter(snakemake.output.exported_variables) as writer:
        iamc.to_excel(writer, sheet_name="data", index=False)
        meta.to_excel(writer, sheet_name="meta", index=False)
