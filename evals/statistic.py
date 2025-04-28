"""Collect statistics for ESM evaluations."""  # noqa: A005

import logging
from functools import partial
from inspect import getmembers
from itertools import product
from pathlib import Path

import pandas as pd
import pypsa
from pandas import DataFrame
from pypsa.statistics import (
    StatisticsAccessor,
    aggregate_timeseries,
    get_weightings,
    groupers,
)

from evals.constants import (
    TRANSMISSION_CARRIER,
    UNITS,
    BusCarrier,
    Carrier,
    DataModel,
    Group,
    Regex,
)
from evals.fileio import get_resources_directory, read_csv_files
from evals.utils import (
    add_grid_lines,
    align_edge_directions,
    filter_by,
    get_trade_type,
    insert_index_level,
    split_location_carrier,
    trade_mask,
)

logger = logging.getLogger(__file__)


def get_location(
    n: pypsa.Network,
    c: str,
    port: str = "",
    location_port: str = "",
    avoid_eu_locations: bool = True,
) -> pd.Series:
    """
    Return the grouper series for the location of a component.

    The additional location port argument will swap the bus
    location to the specified bus port locations. The default
    location is the location from buses at the "port" argument.
    But be careful, the location override will happen for all
    ports of the component.

    Note, that the bus_carrier will still be the bus_carrier
    from the "port" argument, i.e. only the location is swapped.

    Parameters
    ----------
    n
        The network to evaluate.
    c
        The component name, e.g. 'Load', 'Generator', 'Link', etc.
    port
        Limit results to this branch port.
    location_port
        Use the specified port bus for the location, defaults to
        using the location of the 'port' bus.
    avoid_eu_locations
        Look into the port 0 and port 1 location in branch components
        and prefer locations that are not 'EU'.

    Returns
    -------
    :
        A list of series to group statistics by.
    """
    if avoid_eu_locations and c in n.branch_components:
        bus0 = n.static(c)["bus0"].map(n.static("Bus").location).rename("loc0")
        bus1 = n.static(c)["bus1"].map(n.static("Bus").location).rename("loc1")
        buses = pd.concat([bus0, bus1], axis=1)

        def _select_location(row) -> str:
            if row.loc0 != "EU" or pd.isna(row.loc1):
                return row.loc0
            return row.loc1

        return buses.apply(_select_location, axis=1).rename("location")

        # selection order: country code > EU > NaN

    # todo: probably obsolete?
    if location_port and c in n.branch_components:
        buses = n.static(c)[f"bus{location_port}"]
        return buses.map(n.static("Bus").location).rename(DataModel.LOCATION)

    return n.static(c)[f"bus{port}"].map(n.buses.location).rename("location")


def get_location_from_name_at_port(
    n: pypsa.Network, c: str, location_port: str = ""
) -> pd.Series:
    """
    Return the location from the component name.

    Parameters
    ----------
    n
        The network to evaluate.
    c
        The component name, e.g. 'Load', 'Generator', 'Link', etc.
    location_port
        Limit results to this branch port.

    Returns
    -------
    :

    """
    group = f"({Regex.region.pattern})"
    return (
        n.static(c)[f"bus{location_port}"]
        .str.extract(group, expand=False)
        .str.strip()  # some white spaces still go through regex
        .rename(f"bus{location_port}")
    )


def collect_myopic_statistics(
    networks: dict,
    statistic: str,
    aggregate_components: str | None = "sum",
    drop_zero_rows: bool = True,
    **kwargs: object,
) -> pd.DataFrame | pd.Series:
    """
    Build a myopic statistic from loaded networks.

    This method calls ESMStatisticsAccessor methods. It calls the
    statistics method for every year and optionally aggregates
    components, e.g. Links and Lines often should become summed up.

    Parameters
    ----------
    networks
        The loaded networks in a dictionary with the year as keys.
    statistic
        The name of the metric to build.
    aggregate_components
        The aggregation function to combine components by.
    drop_zero_rows
        Whether to drop rows from the returned statistic that have
        only zeros as values.
    **kwargs
        Any key word argument accepted by the statistics function.

    Returns
    -------
    :
        The built statistic with the year as the outermost index level.

    Raises
    ------
    ValueError
        In case a non-existent statistics function was requested.
    """
    kwargs = kwargs or {}

    pypsa_statistics = [m[0] for m in getmembers(pypsa.statistics.StatisticsAccessor)]

    if statistic in pypsa_statistics:  # register a default to reduce verbosity
        kwargs.setdefault("groupby", ["location", "carrier", "bus_carrier", "unit"])

    year_statistics = []
    for year, n in networks.items():
        func = getattr(n.statistics, statistic)
        assert func, (
            f"Statistic '{statistic}' not found. "
            f"Available statistics are: "
            f"'{[m[0] for m in getmembers(n.statistics)]}'."
        )
        year_statistic = func(**kwargs)
        year_statistic = insert_index_level(year_statistic, year, DataModel.YEAR)
        year_statistics.append(year_statistic)

    statistic = pd.concat(year_statistics, axis=0, sort=True)
    if DataModel.LOCATION in statistic.index.names:
        if "EU" in statistic.index.unique(DataModel.LOCATION):
            logger.debug(
                f"EU node found in statistic:\n"
                f"{filter_by(statistic, location='EU')}"
                f"\n\nPlease check if this is intentional!"
            )

    if aggregate_components and "component" in statistic.index.names:
        _names = statistic.index.droplevel("component").names
        statistic = statistic.groupby(_names).agg(aggregate_components)

    if kwargs.get("aggregate_time") is False:
        statistic.columns.name = DataModel.SNAPSHOTS

    if drop_zero_rows and isinstance(statistic, pd.Series):
        statistic = statistic.loc[statistic != 0]
    elif drop_zero_rows and isinstance(statistic, pd.DataFrame):
        statistic = statistic.loc[(statistic != 0).any(axis=1)]
    else:
        raise TypeError(f"Unknown statistic type '{type(statistic)}'")

    # assign the correct unit the statistic if possible
    if "unit" in statistic.index.names:
        if not statistic.empty:
            try:
                statistic.attrs["unit"] = statistic.index.unique("unit").item()
            except ValueError:
                logger.warning(
                    f"Mixed units detected in statistic: {statistic.index.unique('unit')}."
                )
        statistic = statistic.droplevel("unit")

    return statistic.sort_index()


class ESMStatistics(StatisticsAccessor):
    """
    Provides additional statistics for ESM evaluations.

    Extends the StatisticsAccessor with additional metrics.

    Note, that the __call__ method of the base class is not
    updated. Metrics registered with this class need to
    be called explicitly and are not included in the output
    of n.statistics().

    The actual patching is done directly after reading in the
    network files in read_networks(). This means, that
    io.read_networks() must be used to load networks, or the
    statistics will not be available under n.statistics().

    Parameters
    ----------
    n
        The loaded postnetwork.

    result_path
        The output path including the subdirectory, i.e. the path
        where the evaluation results are stored.
    """

    def __init__(self, n: pypsa.Network, result_path: Path) -> None:
        super().__init__(n)
        self.result_path = result_path
        pypsa.options.params.statistics.nice_names = False
        pypsa.options.params.statistics.drop_zero = True
        groupers.add_grouper("location", get_location)
        groupers.add_grouper(
            "bus0", partial(get_location_from_name_at_port, location_port="0")
        )
        groupers.add_grouper(
            "bus1", partial(get_location_from_name_at_port, location_port="1")
        )

    def ac_load_split(self) -> pd.DataFrame:
        """
        Split energy amounts for electricity Loads.

        The following AC loads can be distinguished:
          - industry,
          - rail, and
          - households and services.

        Industry and rail data are read from CSV files.
        HH & services data is the remainder of total
        electricity minus rail and industry parts.

        Returns
        -------
        :
            The data series with split AC loads.

        Notes
        -----
        Currently broken: Energy demands only exist for historical years and
        industry and rail demands probably are not substracted from the correct
        series.
        """
        year = self._n.meta["wildcards"]["planning_horizons"]
        clusters = self._n.meta["wildcards"]["clusters"]
        run = self._n.meta["run"]
        res = get_resources_directory(self._n)(run["prefix"])

        indu = read_csv_files(
            res,
            glob=f"industrial_energy_demand_base_s_{clusters}_*.csv",
            sub_directory=run["name"][0],
        )
        indu = indu.loc[year, "current electricity"] * UNITS["TW"]  # to MWH

        rail = (
            read_csv_files(
                res,
                glob="pop_weighted_energy_totals_s_adm.csv",
                sub_directory=run["name"][0],
            )["electricity rail"]
            * UNITS["TW"]
        )  # fixme: data for 2019 only

        # read_csv_files(
        #     res(run["prefix"]),
        #     glob="pop_weighted_energy_totals_s_adm.csv",
        #     sub_directory=run["name"][0],
        # ).filter(regex="electricity residential|electricity services").sum(axis=1)

        p = (
            self.energy_balance(
                comps="Load",
                groupby=["location", "carrier", "bus_carrier"],
                bus_carrier="low voltage",
            )
            .droplevel(DataModel.BUS_CARRIER)
            .unstack()
        )

        # load p is negative, because it is demand (withdrawal), but csv
        # data (industry, transport) has positive values only. Must
        # reverse the sign for industry and rail demands.
        homes_and_trade = Carrier.domestic_homes_and_trade
        p[Carrier.industry] = indu.mul(-1)
        p[Carrier.electricity_rail] = rail.mul(-1)
        p[homes_and_trade] = p["electricity"] + indu + rail

        if any(p[homes_and_trade] > 0):
            logger.warning(
                msg=f"Positive values found for {homes_and_trade} "
                f"demand. This happens if the combined electricity demand "
                f"from Industry and Rail nodal energy files is larger than "
                f"the electricity Loads in the network.\n"
                f"{p[p[homes_and_trade] > 0][homes_and_trade]}\n\n"
                f"All values larger than zero will be set to zero. "
                f"(Note that this is different to the Toolbox implementation "
                f"where signs are flipped).\n"
            )
            # fixme: just a note. There is a bug in the old Toolbox that
            #  counts the aforementioned amounts as demand (although
            #  the amounts should be clipped.)
            p[homes_and_trade] = p[homes_and_trade].clip(upper=0)

        # rename to avoid mixing up electricity with
        p = p.rename({"electricity": "industry + hh & services load"}, axis=1)

        df = insert_index_level(p.stack(), "low voltage", DataModel.BUS_CARRIER)
        df = df.reorder_levels(DataModel.IDX_NAMES)

        df.attrs["name"] = "Electricity split "
        df.attrs["unit"] = "MWh"

        return df

    def bev_v2g(self, drop_v2g_withdrawal: bool = True) -> DataFrame:
        """
        Calculate BEV and V2G energy amounts.

        Parameters
        ----------
        drop_v2g_withdrawal
            Whether to exclude vehicle to grid technologies from the
            results. This option is included since the Toolbox
            implementation drops them too.

        Returns
        -------
        :
            A DataFrame containing the calculated BEV and V2G energy
            amounts.
        """
        c = Carrier
        names_supply = {
            c.bev_charger: c.bev_charger_supply,
            c.v2g: c.v2g_supply,
        }
        names_withdrawal = {
            c.bev: c.bev_passenger_withdrawal,
            c.bev_charger: c.bev_charger_draw,
            c.v2g: c.v2g_withdrawal,
        }
        carrier = [Carrier.bev, Carrier.bev_charger, Carrier.v2g]
        supply = self.supply(
            comps="Link",
            groupby=["location", "carrier", "bus_carrier"],
            bus_carrier=[BusCarrier.AC, BusCarrier.LI_ION],
        )
        supply = filter_by(supply, carrier=carrier)

        withdrawal = self.withdrawal(
            comps="Link",
            groupby=["location", "carrier", "bus_carrier"],
            bus_carrier=[BusCarrier.AC, BusCarrier.LI_ION],
        )
        withdrawal = filter_by(withdrawal, carrier=carrier)
        withdrawal = withdrawal.mul(-1)  # to keep withdrawal negative

        # rename carrier to avoid name clashes for supply/withdrawal
        supply = supply.rename(names_supply, level=DataModel.CARRIER)
        withdrawal = withdrawal.rename(names_withdrawal, level=DataModel.CARRIER)

        # join along index, sum duplicates and pivot carriers to columns
        p = (
            pd.concat([withdrawal, supply])
            .groupby([DataModel.LOCATION, DataModel.CARRIER])
            .sum()
            .unstack()
        )

        ratio = (p[c.bev_charger_draw] / p[c.bev_charger_supply]).abs()

        p[c.bev_charger_losses] = p[c.bev_charger_draw] + p[c.bev_charger_supply]
        p[c.bev_demand] = ratio * p[c.bev_passenger_withdrawal]
        p[c.bev_losses] = p[c.bev_demand] - p[c.bev_passenger_withdrawal]
        p[c.v2g_demand] = ratio * p[c.v2g_withdrawal] if c.v2g_withdrawal in p else 0
        p[c.v2g_losses] = p[c.v2g_demand] + p[c.v2g_supply] if c.v2g_supply in p else 0

        ser = insert_index_level(p.stack(), BusCarrier.AC, DataModel.BUS_CARRIER, pos=2)
        ser.attrs["name"] = "BEV&V2G"
        ser.attrs["unit"] = "MWh"

        if drop_v2g_withdrawal:
            ser = ser.drop(c.v2g_withdrawal, level=DataModel.CARRIER, errors="ignore")

        return ser

    def phs_split(
        self, aggregate_time: str = "sum", drop_hydro_cols: bool = True
    ) -> pd.DataFrame:
        """
        Split energy amounts for StorageUnits.

        Parameters
        ----------
        aggregate_time
            The aggregation function used to aggregate time steps.

        drop_hydro_cols
            Whether, or not to drop 'hydro' carriers from the result.
            This is required to stay consistent with the old Toolbox
            implementation.

        Returns
        -------
        :
            A DataFrame containing the split energy amounts for
            PHS and hydro.
        """
        n = self._n

        idx = n.static("StorageUnit").index
        phs = pd.DataFrame(index=idx)
        for time_series in ("p_dispatch", "p_store", "spill", "inflow"):
            p = n.pnl("StorageUnit")[time_series].reindex(columns=idx, fill_value=0)
            weights = get_weightings(n, "StorageUnit")
            phs[time_series] = aggregate_timeseries(p, weights, agg=aggregate_time)

        efficiency = phs["p_store"] * n.static("StorageUnit")["efficiency_dispatch"]
        part_inflow = phs["inflow"] / (phs["inflow"] + efficiency)

        phs["Dispatched Power from Inflow"] = phs["p_dispatch"] * part_inflow
        phs["Dispatched Power from Stored"] = phs["p_dispatch"] * (1 - part_inflow)
        phs["Spill from Inflow"] = phs["spill"] * part_inflow
        phs["Spill from Stored"] = phs["spill"] * (1 - part_inflow)

        # use evaluation output carrier names
        mapper = {
            "p_dispatch": "Dispatched Power",
            "p_store": "Stored Power",
            "inflow": "Inflow",
            "spill": "Spill",
        }
        phs = phs.rename(mapper, axis=1)

        ser = phs.stack()
        ser.index = ser.index.swaplevel(0, 1)
        ser.index = split_location_carrier(ser.index, names=DataModel.IDX_NAMES)

        # merge 'carrier' with 'bus_carrier' level and keep original
        # bus_carrier. Needed to stay consistent with the old Toolbox
        # naming conventions.
        ser.index = pd.MultiIndex.from_tuples(
            [(r[1], f"{r[2]} {r[0]}", r[2]) for r in ser.index],
            names=DataModel.IDX_NAMES,
        )

        ser = ser.rename(
            index={"PHS": BusCarrier.AC, "hydro": BusCarrier.AC},
            level=DataModel.BUS_CARRIER,
        )

        ser.attrs["name"] = "PHS&Hydro"
        ser.attrs["unit"] = "MWh"

        if drop_hydro_cols:
            cols = [
                "hydro Dispatched Power from Inflow",
                "hydro Dispatched Power from Stored",
                "hydro Spill from Inflow",
                "hydro Spill from Stored",
            ]
            ser = ser.drop(cols, level=DataModel.CARRIER)

        return ser.sort_index()

    def phs_hydro_operation(self) -> pd.DataFrame:
        """
        Calculate Hydro- and Pumped Hydro Storage unit statistics.

        Returns
        -------
        :
            Cumulated or constant time series for storage units.
        """
        n = self._n
        ts_efficiency_name_agg = [
            ("p_dispatch", "efficiency_dispatch", Group.turbine_cum, "cumsum"),
            ("p_store", "efficiency_store", Group.pumping_cum, "cumsum"),
            ("spill", None, Group.spill_cum, "cumsum"),
            ("inflow", None, Group.inflow_cum, "cumsum"),
            ("state_of_charge", None, Group.soc, None),
        ]

        weights = get_weightings(n, "StorageUnit")

        su = n.static("StorageUnit").query("carrier in ['PHS', 'hydro']")

        results = []
        for time_series, efficiency, index_name, agg in ts_efficiency_name_agg:
            df = n.pnl("StorageUnit")[time_series].filter(su.index, axis=1)
            if agg:
                df = df.mul(weights, axis=0).agg(agg)
            if efficiency == "efficiency_dispatch":
                df = df / su[efficiency]
            elif efficiency == "efficiency_store":
                df = df * su[efficiency]
            # The actual bus carrier is "AC" for both, PHS and hydro.
            # Since only PHS and hydro are considered, we can use the
            # bus_carrier level.
            result = insert_index_level(df, index_name, DataModel.BUS_CARRIER, axis=1)
            results.append(result.T)

        # broadcast storage volume to time series (not quite the
        # same as utils.scalar_to_time_series, because it's a series)
        volume = su["p_nom_opt"] * su["max_hours"]
        volume_ts = pd.concat([volume] * len(n.snapshots), axis=1)
        volume_ts.columns = n.snapshots
        volume_ts = insert_index_level(volume_ts, Group.soc_max, DataModel.BUS_CARRIER)
        results.append(volume_ts)

        statistic = pd.concat(results)
        statistic.index = split_location_carrier(
            statistic.index,
            names=[DataModel.BUS_CARRIER, DataModel.LOCATION, DataModel.CARRIER],
        )
        statistic = statistic.reorder_levels(DataModel.IDX_NAMES)

        statistic.columns.names = [DataModel.SNAPSHOTS]
        statistic.attrs["name"] = "StorageUnit Operation"
        statistic.attrs["unit"] = "MWh"

        return statistic

    def trade_energy(
        self,
        scope: str | tuple,
        direction: str = "saldo",
        bus_carrier: str = None,
        aggregate_time: str = "sum",
    ) -> pd.DataFrame:
        """
        Calculate energy amounts exchanged between locations.

        Returns positive values for 'import' (supply) and negative
        values for 'export' (withdrawal).

        Parameters
        ----------
        scope
            The scope of energy exchange. Must be one of "foreign",
            "domestic", or "local".

        direction
            The direction of the trade. Can be one of "saldo", "export",
            or "import".

        bus_carrier
            The bus carrier for which to calculate the energy exchange.
            Defaults to using all bus carrier.

        aggregate_time
            The method of aggregating the energy exchange over time.
            Can be one of "sum", "mean", "max", "min".

        Returns
        -------
        :
            A DataFrame containing the calculated energy exchange
            between locations.
        """
        n = self._n
        results_comp = []

        buses = n.static("Bus").reset_index()
        if bus_carrier:
            _bc = [bus_carrier] if isinstance(bus_carrier, str) else bus_carrier
            buses = buses.query("carrier in @_bc")

        for port, c in product((0, 1), ("Link", "Line")):
            mask = trade_mask(n.static(c), scope).to_numpy()
            comp = n.static(c)[mask].reset_index()

            p = buses.merge(
                comp,
                left_on="Bus",
                right_on=f"bus{port}",
                suffixes=("_bus", ""),
            ).merge(n.pnl(c).get(f"p{port}").T, on=c)

            _location = (
                DataModel.LOCATION + "_bus"
                if "location" in comp
                else DataModel.LOCATION
            )
            p = p.set_index([_location, DataModel.CARRIER, "carrier_bus"])
            p.index.names = DataModel.IDX_NAMES
            # branch components have reversed sign
            p = p.filter(n.snapshots, axis=1).mul(-1.0)
            if direction == "export":
                p = p.clip(upper=0)  # keep negative values (withdrawal)
            elif direction == "import":
                p = p.clip(lower=0)  # keep positive values (supply)
            elif direction != "saldo":
                raise ValueError(f"Direction '{direction}' not supported.")

            results_comp.append(insert_index_level(p, c, "component"))

        result = pd.concat(results_comp).groupby(DataModel.IDX_NAMES).sum()

        if aggregate_time:
            # assuming Link and Line have the same weights
            weights = get_weightings(n, "Link")
            result = result.multiply(weights, axis=1)
            result = result.agg(aggregate_time, axis=1)

        name = " & ".join(scope) if isinstance(scope, tuple) else scope
        result.attrs["name"] = f"{name} {direction}"
        result.attrs["unit"] = "MWh"

        return result.sort_index()

    def trade_capacity(
        self,
        scope: str,
        bus_carrier: str = "",
    ) -> pd.DataFrame:
        """
        Calculate exchange capacity between locations.

        Parameters
        ----------
        scope
            The scope of energy exchange. Must be one of
            constants.TRADE_TYPES.
        bus_carrier
            The bus carrier for which to calculate the energy exchange.
            Defaults to using all bus carrier.

        Returns
        -------
        :
            Energy exchange capacity between locations.
        """
        n = self._n

        capacity = self.optimal_capacity(
            comps=n.branch_components,
            bus_carrier=bus_carrier,
            groupby=["bus0", "bus1", "carrier", "bus_carrier"],
            nice_names=False,
        ).to_frame()
        trade_type = capacity.apply(
            lambda row: get_trade_type(row.name[1], row.name[2]), axis=1
        )

        trade_capacity = capacity[trade_type == scope]

        # duplicate capacities to list them for source and destination
        # locations. For example, the trade capacity for AT -> DE gas
        # pipeline will be shown in location AT and in location DE.
        df_list = []
        for bus in ("bus0", "bus1"):
            df = trade_capacity.droplevel(bus)
            df.index.names = [DataModel.COMPONENT] + DataModel.IDX_NAMES
            df_list.append(df)

        trade_capacity = pd.concat(df_list).drop_duplicates()

        return trade_capacity.squeeze()

    def ambient_heat(self) -> pd.Series | pd.DataFrame:
        """Calculate ambient heat energy amounts used by heat pumps."""
        energy_balance = self._n.statistics.energy_balance(
            comps="Link",
            bus_carrier=BusCarrier.heat_buses() + ["low voltage"],
            groupby=DataModel.IDX_NAMES,
        )
        heat_pump = energy_balance.filter(like="heat pump", axis=0)

        def _heat_minus_ac(ser: pd.Series) -> pd.Series:
            """
            Return the sum of AC withdrawal and heat supply.

            This function assumes, that 1 MWh (electricity) is 1 MWh (thermal).

            Parameters
            ----------
            ser
                The input Series with location, carrier and
                bus_carrier MultiIndex levels. The carrier is
                expected to be one of the heat pump carriers.
                The bus_carrier per carrier are AC and one of
                the heat bus_carrier.

            Returns
            -------
            :
                The sum of AC withdrawal and heat supply per
                carrier.
            """
            hp = ser.unstack(DataModel.BUS_CARRIER)
            assert hp.shape[1] == 2, f"Unexpected number of bus_carrier: {hp.columns}."
            assert "low voltage" in hp.columns, (
                f"AC missing in bus_carrier: {hp.columns}."
            )
            return hp.T.sum()

        ambient_heat = heat_pump.groupby(DataModel.CARRIER, group_keys=False).apply(
            _heat_minus_ac
        )

        new_index_items = []
        for loc, carr in ambient_heat.index:
            if carr.startswith("rural"):
                new_index_items.append((loc, carr, BusCarrier.HEAT_RURAL))
            elif carr.startswith("urban decentral"):
                new_index_items.append((loc, carr, BusCarrier.HEAT_URBAN_DECENTRAL))
            elif carr.startswith("urban central"):
                new_index_items.append((loc, carr, BusCarrier.HEAT_URBAN_CENTRAL))
            else:
                raise ValueError(f"Carrier {carr} not recognized.")
        ambient_heat.index = pd.MultiIndex.from_tuples(
            new_index_items, names=DataModel.IDX_NAMES
        )

        # def _add_bus_carrier_from_carrier_name(idx):
        #     """"""
        #     loc, carr = idx
        #     if carr.startswith("rural"):
        #         bus_carr = BusCarrier.HEAT_RURAL
        #     elif carr.startswith("urban decentral"):
        #         bus_carr = BusCarrier.HEAT_URBAN_DECENTRAL
        #     elif carr.startswith("urban central"):
        #         bus_carr = BusCarrier.HEAT_URBAN_CENTRAL
        #     else:
        #         raise ValueError(f"Carrier {carr} not recognized.")
        #     # return (loc, carr, bus_carr)
        #     return bus_carr
        #
        # bus_carrier = ambient_heat.index.map(_add_bus_carrier_from_carrier_name)
        # bus_carrier.name = DataModel.BUS_CARRIER
        #
        # ambient_heat = insert_index_level(
        #     ambient_heat, "ambient heat", DataModel.BUS_CARRIER, pos=2
        # )

        ambient_heat.attrs["name"] = "Ambient Heat"
        ambient_heat.attrs["unit"] = "MWh"

        return ambient_heat

    def grid_capactiy(
        self,
        comps: list = None,
        bus_carrier: list = None,
        carrier: list = None,
        append_grid: bool = True,
        align_edges: bool = True,
    ) -> pd.DataFrame:
        """
        Return transmission grid capacities.

        Parameters
        ----------
        comps
            The network components to consider, defaults to all
            pypsa.Networks.branch_components.
        bus_carrier
            The bus carrier to consider.
        carrier
            The carrier to consider, defaults to all
            constants.TRANSMISSION_CARRIER.
        append_grid
            Whether to add the grid lines to the result.
        align_edges
            Whether to adjust edges between the same nodes but in
            reversed direction. For example, AC and DC grids have
            edges between IT0 0 and FR0 0 as IT->FR and FR->IT,
            respectively. If enabled, both will have the same bus0 and
            bus1.

        Returns
        -------
        :
            The optimal capacity for transmission technologies between
            nodes.

        Notes
        -----
        The "pypsa.statistics.transmission" statistic does not work here
        because it returns energy amounts whereas this statistic returns
        the optimal capacity.
        """
        n = self._n
        carrier = carrier or list(TRANSMISSION_CARRIER)
        capacities = n.statistics.optimal_capacity(
            comps=comps or n.branch_components,
            bus_carrier=bus_carrier,
            groupby=["bus0", "bus1", "carrier", "bus_carrier"],
        )
        result = filter_by(capacities, carrier=carrier)

        result.attrs["name"] = "Capacity"
        result.attrs["unit"] = "MW"
        result.name = f"{result.attrs['name']} ({result.attrs['unit']})"

        if align_edges:
            result = align_edge_directions(result)

        if append_grid:
            result = add_grid_lines(n.static("Bus"), result)

        return result.sort_index()

    def grid_flow(
        self,
        comps: list = None,
        bus_carrier: list = None,
        carrier: list = None,
        aggregate_time: str = "sum",
        append_grid: bool = True,
    ) -> pd.DataFrame:
        """
        Return the transmission grid energy flow.

        Parameters
        ----------
        comps
            The network components to consider, defaults to all
            pypsa.Networks.branch_components.
        bus_carrier
            The bus carrier to consider.
        carrier
            The carrier to consider, defaults to all
            constants.TRANSMISSION_CARRIER.
        aggregate_time
            The aggregation function aggregate by.
        append_grid
            Whether to add the grid lines to the result.

        Returns
        -------
        :
            The amount of energy transfer for transmission technologies
            between nodes.
        """
        n = self._n
        carrier = carrier or list(TRANSMISSION_CARRIER)
        comps = comps or n.branch_components

        energy_transmission = n.statistics.transmission(
            comps=comps,
            groupby=["bus0", "bus1", "carrier", "bus_carrier"],
            bus_carrier=bus_carrier,
            aggregate_time=False,
        )
        energy_transmission = filter_by(energy_transmission, carrier=carrier)

        # split directions:
        # positive values are from bus0 to bus1
        bus0_to_bus1 = energy_transmission.clip(lower=0)

        # negative values are from bus1 to bus0
        idx_names = list(energy_transmission.index.names)
        bus1_to_bus0 = energy_transmission.clip(upper=0).mul(-1)
        # we reverse the node index levels to show positive values and
        # have a consistent way of interpreting the energy flow
        bus1_to_bus0 = bus1_to_bus0.swaplevel("bus0", "bus1")
        pos0, pos_1 = idx_names.index("bus0"), idx_names.index("bus1")
        idx_names[pos_1], idx_names[pos0] = idx_names[pos0], idx_names[pos_1]
        bus1_to_bus0.index.names = idx_names

        result = pd.concat([bus0_to_bus1, bus1_to_bus0])
        result = result.groupby(idx_names).sum()

        assert aggregate_time, "Time Series is not supported."
        unit = "MW"
        if aggregate_time in ("max", "min"):
            result = result.agg(aggregate_time, axis=1)
        elif aggregate_time:  # mean, median, etc.
            weights = get_weightings(n, comps)
            result = result.mul(weights, axis=1).agg(aggregate_time, axis=1)
            unit = "MWh"

        result.attrs["name"] = "Energy"
        result.attrs["unit"] = unit
        result.name = f"{result.attrs['name']} ({result.attrs['unit']})"

        if append_grid:
            result = add_grid_lines(n, result)

        return result.sort_index()
