# -*- coding: utf-8 -*-
"""Collect statistics for ESM evaluations."""  # noqa: A005

from inspect import getmembers
from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd
import pypsa
from pandas import DataFrame
from pypsa.statistics import (
    StatisticsAccessor,
    aggregate_timeseries,
    get_bus_and_carrier,
    get_operation,
    get_weightings,
    port_efficiency,
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
from evals.fileio import read_csv_files
from evals.metric import logger
from evals.utils import (
    filter_by,
    get_trade_type,
    insert_index_level,
    split_location_carrier,
    trade_mask,
)


def get_location(
    n: pypsa.Network, c: str, port: str = "", location_port: str = ""
) -> pd.Series:
    """Return the grouper series for the location of a component.

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

    Returns
    -------
    :
        A list of series to group statistics by.
    """
    if location_port and c in n.branch_components:
        bus_location = n.static(c)[f"bus{location_port}"]
        return bus_location.map(n.static("Bus").location).rename(DataModel.LOCATION)

    return n.static(c)[f"bus{port}"].map(n.buses.location).rename("location")


def get_location_and_carrier_and_bus_carrier(
    n: pypsa.Network,
    c: str,
    port: str = "",
    nice_names: bool = False,
    location_port: str = "",
) -> list[pd.Series]:
    """Get location, carrier, and bus carrier.

    Used in groupby statements to group statistics results.

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
    nice_names
        Whether to return the carrier alias.
    location_port
        Use the specified port bus for the location, defaults to
        using the location of the 'port' bus.

    Returns
    -------
    :
        A list of series to group statistics by.
    """
    bus, carrier = get_bus_and_carrier(n, c, port, nice_names=nice_names)

    if location_port and c in n.branch_components:
        bus_location = n.df(c)[f"bus{location_port}"]
        location = bus_location.map(n.df("Bus").location).rename(DataModel.LOCATION)
    else:
        location = bus.map(n.df("Bus").location).rename(DataModel.LOCATION)

    # Note, that we still use the original bus carrier from 'bus' here,
    # even if location_port is being used.
    bus_carrier = bus.map(n.df("Bus").carrier).rename(DataModel.BUS_CARRIER)

    return [location, carrier, bus_carrier]


def get_buses_and_carrier_and_bus_carrier(
    n: pypsa.Network, c: str, port: str = "", nice_names: bool = True
) -> list[pd.Series]:
    """Get src_bus, dst_bus, carrier, and bus carrier.

    Used in groupby statements to group statistics results.

    Parameters
    ----------
    n
        The network to evaluate.
    c
        The component name, e.g. 'Load', 'Generator', 'Link', etc.
    port
        Limit results to this port.
    nice_names
        Whether to return the carrier alias, defaults to True.

    Returns
    -------
    :
        A list of series to group statistics by.
    """
    pat = f"({Regex.region.pattern})"
    bus0 = n.df(c)["bus0"].str.extract(pat, expand=False)
    bus1 = n.df(c)["bus1"].str.extract(pat, expand=False)
    bus, carrier = get_bus_and_carrier(n, c, port, nice_names=nice_names)
    bus_carrier = bus.map(n.df("Bus").carrier).rename(DataModel.BUS_CARRIER)
    return [bus0, bus1, carrier, bus_carrier]


def collect_myopic_statistics(
    networks: dict,
    statistic: str,
    aggregate_components: str = "sum",
    carrier: list | tuple = None,
    drop_zero_rows: bool = True,
    **kwargs: object,
) -> pd.DataFrame | pd.Series:
    """Build a myopic statistic from loaded networks.

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
    carrier
        A list of carrier names used to filter the statistic. Only
        carrier in the input will be returned. Returns all by default.
    drop_zero_rows
        Whether to drop rows from the returned statistic that have
        only zeros as values.
    **kwargs
        Any key word argument accepted by the metric method.

    Returns
    -------
    :
        The built metric with the year as the outermost index level.

    Raises
    ------
    ValueError
        In case a non-existent metric was requested.
    """
    kwargs = kwargs or {}

    pypsa_statistics = [m[0] for m in getmembers(pypsa.statistics.StatisticsAccessor)]

    # groupby custom function by default to reduce visual noise
    if statistic in pypsa_statistics:
        # PyPSA 0.32 support registering grouper functions
        # https://github.com/PyPSA/PyPSA/pull/1078
        # kwargs.setdefault("groupby", get_location_and_carrier_and_bus_carrier)
        kwargs.setdefault("groupby", ["location", "carrier", "bus_carrier"])

    year_statistics = []
    for year, n in networks.items():
        statistic_func = getattr(n.statistics, statistic)
        if not statistic_func:
            raise ValueError(
                f"Statistic '{statistic}' not found. "
                f"Available statistics are: "
                f"'{[m[0] for m in getmembers(n.statistics)]}'."
            )
        year_statistic = statistic_func(**kwargs)
        year_statistic = insert_index_level(year_statistic, year, DataModel.YEAR)
        year_statistics.append(year_statistic)

    statistic = pd.concat(year_statistics, axis=0, sort=True)

    if aggregate_components and "component" in statistic.index.names:
        _names = list(statistic.index.names)
        _names.remove("component")
        statistic = statistic.groupby(_names).agg(aggregate_components)

    if carrier:
        statistic = filter_by(statistic, carrier=carrier)

    if kwargs.get("aggregate_time") is False:
        statistic.columns.name = DataModel.SNAPSHOTS

    # todo: verify zeros are being dropped out of the box
    # # drop entries with all zero rows. They only clutter results.
    # if drop_zero_rows:
    #     statistic = (
    #         statistic.loc[statistic != 0]  # Series
    #         if isinstance(statistic, pd.Series)
    #         else statistic.loc[(statistic != 0).any(axis=1)]  # DataFrame
    #     )

    return statistic.sort_index()


class ESMStatistics(StatisticsAccessor):
    """Provides additional statistics for ESM evaluations.

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

        # configure statistics:
        self.set_parameters(nice_names=False)
        # register grouper here once PyPSA >= 0.32

    def ac_load_split(self) -> pd.DataFrame:
        """Split energy amounts for electricity Loads.

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
        """
        year = self.n.meta["wildcards"]["planning_horizon"]

        indu = read_csv_files(
            self.result_path,
            glob="industrial_demand_elec*.csv",
            sub_directory="esm_run/interpolated_data",
        )
        indu = indu.loc[year, "current electricity"] * UNITS["TW"]  # to MWH

        rail = read_csv_files(
            self.result_path,
            glob="nodal_energy_totals_*.csv",
            sub_directory="esm_run/resources",
        )
        rail = rail.loc[year, "electricity rail"] * UNITS["TW"]  # to MWH

        p = self.energy_balance(
            comps="Load",
            groupby=get_location_and_carrier_and_bus_carrier,
            bus_carrier=BusCarrier.AC,
            nice_names=False,
        )
        p = p.droplevel(DataModel.BUS_CARRIER).unstack()

        # load p is negative, because it is demand (withdrawal), but csv
        # data (industry, transport) has positive values only. Must
        # reverse the sign for industry and rail demands.
        p[Carrier.industry] = indu.mul(-1)
        p[Carrier.electricity_rail] = rail.mul(-1)
        p[Carrier.domestic_homes_and_trade] = p["electricity"] + indu + rail

        if any(p[Carrier.domestic_homes_and_trade] > 0):
            logger.warning(
                msg=f"Positive values found for {Carrier.domestic_homes_and_trade} "
                f"demand. This happens if the combined electricity demand "
                f"from Industry and Rail nodal energy files is larger than "
                f"the electricity Loads in the network.\n"
                f"{p[p[Carrier.domestic_homes_and_trade] > 0]
                    [Carrier.domestic_homes_and_trade]}\n\n"
                f"All values larger than zero will be set to zero. "
                f"(Note that this is different to the Toolbox implementation "
                f"where signs are flipped).\n"
            )
            # fixme: just a note. There is a bug in the old Toolbox that
            #  counts the a aforementioned amounts as demand (although
            #  the amounts should be clipped.)
            p[Carrier.domestic_homes_and_trade] = p[
                Carrier.domestic_homes_and_trade
            ].clip(upper=0)

        p = p.rename({"electricity": "industry + hh & services load"}, axis=1)

        df = insert_index_level(p.stack(), BusCarrier.AC, DataModel.BUS_CARRIER)
        df = df.reorder_levels(DataModel.IDX_NAMES)

        df.attrs["name"] = "Energy "
        df.attrs["unit"] = "MWh"

        return df

    def bev_v2g(self, drop_v2g_withdrawal: bool = True) -> DataFrame:
        """Calculate BEV and V2G energy amounts.

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
            groupby=get_location_and_carrier_and_bus_carrier,
            bus_carrier=[BusCarrier.AC, BusCarrier.LI_ION],
        )
        supply = filter_by(supply, carrier=carrier)

        withdrawal = self.withdrawal(
            comps="Link",
            groupby=get_location_and_carrier_and_bus_carrier,
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
        """Split energy amounts for StorageUnits.

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
        n = self.n

        idx = n.df("StorageUnit").index
        phs = pd.DataFrame(index=idx)
        for time_series in ("p_dispatch", "p_store", "spill", "inflow"):
            p = n.pnl("StorageUnit")[time_series].reindex(columns=idx, fill_value=0)
            weights = get_weightings(n, "StorageUnit")
            phs[time_series] = aggregate_timeseries(p, weights, agg=aggregate_time)

        efficiency = phs["p_store"] * n.df("StorageUnit")["efficiency_dispatch"]
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
        """Calculate Hydro- and Pumped Hydro Storage unit statistics.

        Returns
        -------
        :
            Cumulated or constant time series for storage units.
        """
        n = self.n
        ts_efficiency_name_agg = [
            ("p_dispatch", "efficiency_dispatch", Group.turbine_cum, "cumsum"),
            ("p_store", "efficiency_store", Group.pumping_cum, "cumsum"),
            ("spill", None, Group.spill_cum, "cumsum"),
            ("inflow", None, Group.inflow_cum, "cumsum"),
            ("state_of_charge", None, Group.soc, None),
        ]

        weights = get_weightings(n, "StorageUnit")

        su = n.df("StorageUnit").query("carrier in ['PHS', 'hydro']")

        results = []
        for time_series, efficiency, index_name, agg in ts_efficiency_name_agg:
            df = n.pnl("StorageUnit")[time_series].filter(su.index, axis=1)
            if agg:
                df = df.mul(weights, axis=0).agg(agg)
            if efficiency == "efficiency_dispatch":
                df = df / su[efficiency]
            elif efficiency == "efficiency_store":
                df = df * su[efficiency]
            # The actual bus carrier in "AC" for both, PHS and hydro.
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
        """Calculate energy amounts exchanged between locations.

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
        n = self.n
        results_comp = []

        buses = n.df("Bus").reset_index()
        if bus_carrier:
            _bc = [bus_carrier] if isinstance(bus_carrier, str) else bus_carrier
            buses = buses.query("carrier in @_bc")

        for port, c in product((0, 1), ("Link", "Line")):
            mask = trade_mask(n.df(c), scope).to_numpy()
            comp = n.df(c)[mask].reset_index()

            p = buses.merge(
                comp,
                left_on="Bus",
                right_on=f"bus{port}",
                suffixes=("_bus", ""),
            ).merge(n.pnl(c).get(f"p{port}").T, on=c)

            p = p.set_index([DataModel.LOCATION, DataModel.CARRIER, "carrier_bus"])
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
        """Calculate exchange capacity between locations.

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
        n = self.n

        capacity = self.optimal_capacity(
            comps=n.branch_components,
            bus_carrier=bus_carrier,
            groupby=get_buses_and_carrier_and_bus_carrier,
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

    def energy_input(
        self,
        comps: str,
        aggregate_time: str = "sum",
        bus_carrier: list = None,
        include_losses: bool = False,
        location_port: str = "",
    ) -> pd.DataFrame:
        """
        Calculate the energy needed to produce another form of energy.

        This statistic calculates the amount of energy needed to
        produce energy of a certain kind, e.g. how much AC, gas, H2, or
        biomass is needed to produce heat? The function calculates the
        amount of energy for every branch port where the requested
        bus_carrier (heat in the example) is connected. It will apply
        the respective efficiency share to return the amount of energy
        in the form of the input branch carrier, i.e. the bus0
        bus_carrier, needed to produce the energy at the branch.

        The bus0 branch is also referred to the source branch, and the
        branch with the requested bus_carrier energies is referred to
        as the target branch.

        Beware, that this statistic may not be correct for CO2?

        Parameters
        ----------
        comps
            The PyPSA component (singular) to analyse.
            "comps" is kept in favor of "comp" for API consistency.

        aggregate_time
            The aggregation function to combine snapshots, defaults
            to "sum".

        bus_carrier
            The bus carrier to calculate energy demand for, i.e. the
            target branch bus carrier (heat in the example).

        include_losses
            Flag to return the full amount of Link energy input at
            branch port 1. This is useful to return the final
            energy demand for Links where the requested bus_carrier is
            connected to port 1.

        location_port
            Use the location from the specified bus port instead of the
            bus0 location.

        Returns
        -------
        :
            Energy amounts per carrier needed to produce energy for the
            input bus_carrier.
        """
        n = self.n

        def _calculate_efficiency_share() -> float:
            """Calculate the efficiency fraction for a branch port.

            Separate energy needed to produce energy at the target
            branch port from energy needed to produce energy at
            other branch ports of the same Link, e.g. Electrolysis
            takes AC (port 0) and produces H2 (port 1) and heat
            (port 2). We only want to know how much AC is needed
            to produce H2 or heat.

            If include_losses is True, the efficiency fraction is
            skipped for port "1" branches and the full amount of energy
            needed to produce energy of branch "1" will ultimately be
            returned.

            Returns
            -------
            The efficiency share relative to the bus0 energies for the
            target branch port.

            Notes
            -----
            Efficiencies of branch ports > 1 are relative to the first
            branch 1 efficiency.
            https://pypsa.readthedocs.io/en/latest/user-guide/components.html#multilink
            """
            if port == "0":
                return port_efficiency(n, comps, port=port)  # -1.0
            elif port == "1" and include_losses:
                return 1.0
            elif port == "1" and not include_losses:
                return port_efficiency(n, comps, port=port)

            eff_port_1 = port_efficiency(n, comps, port="1")
            eff_target = port_efficiency(n, comps, port=port)
            return eff_target / (eff_port_1 + eff_target)

        buses = n.df("Bus").reset_index()
        if bus_carrier:
            buses = buses.query("carrier in @bus_carrier")

        comp = n.df(comps).reset_index()

        ports = [col[3:] for col in n.df(comps).filter(like="bus")]

        port_results = []
        for port in ports:
            time_series = get_operation(n, comps).T
            efficiency_share = _calculate_efficiency_share()
            time_series = time_series.mul(efficiency_share, axis=0)

            bus_comp = buses.merge(
                comp, left_on="Bus", right_on=f"bus{port}", suffixes=("_bus", "")
            )
            p = bus_comp.merge(time_series, on=comps)

            # bus_carrier column needs to be corrected from target
            # branch carriers to the bus_carrier of energy
            # withdrawal, i.e. the bus0 (=input) bus_carrier.
            bus = "bus0" if comps in n.branch_components else "bus"
            p[DataModel.BUS_CARRIER] = p[bus].map(n.df("Bus")[DataModel.CARRIER])

            # support location switching from EU to country nodes
            if location_port and comps in n.branch_components:
                p[DataModel.LOCATION] = p[f"bus{location_port}"].map(
                    n.df("Bus")[DataModel.LOCATION]
                )

            carrier_col = "type" if comps == "Line" else DataModel.CARRIER
            p = p.set_index([DataModel.LOCATION, carrier_col, DataModel.BUS_CARRIER])
            p.index.names = DataModel.IDX_NAMES
            p = p.filter(n.snapshots, axis=1)

            port_results.append(insert_index_level(p, comps, DataModel.COMPONENT))

        result = pd.concat(port_results).groupby(DataModel.IDX_NAMES).sum()

        unit = "MW"
        if aggregate_time in ("max", "min"):
            result = result.agg(aggregate_time, axis=1)
        elif aggregate_time:  # mean, median, etc.
            weights = get_weightings(n, comps)
            result = result.mul(weights, axis=1).agg(aggregate_time, axis=1)
            unit = "MWh"

        result.attrs["name"] = "Input Energy"
        result.attrs["unit"] = unit

        return result.sort_index()

    def ambient_heat(self) -> pd.Series | pd.DataFrame:
        """Calculate ambient heat energy amounts used by heat pumps."""
        energy_balance = self.n.statistics.energy_balance(
            comps="Link",
            bus_carrier=[
                "residential rural heat",
                "services rural heat",
                "urban central heat",
                "AC",
            ],
            groupby=get_location_and_carrier_and_bus_carrier,
        )
        heat_pump = energy_balance.filter(like="heat pump", axis=0)

        def _heat_minus_ac(ser: pd.Series) -> pd.Series:
            """Return the sum of AC withdrawal and heat supply.

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
            assert "AC" in hp.columns, f"AC missing in bus_carrier: {hp.columns}."
            return hp.T.sum()

        ambient_heat = heat_pump.groupby(DataModel.CARRIER, group_keys=False).apply(
            _heat_minus_ac
        )

        ambient_heat = insert_index_level(
            ambient_heat, "ambient heat", DataModel.BUS_CARRIER, pos=2
        )

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
        """Return transmission grid capacities.

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
        n = self.n
        carrier = carrier or list(TRANSMISSION_CARRIER)
        capacities = n.statistics.optimal_capacity(
            comps=comps or n.branch_components,
            bus_carrier=bus_carrier,
            groupby=get_buses_and_carrier_and_bus_carrier,
        )
        result = filter_by(capacities, carrier=carrier)

        result.attrs["name"] = "Capacity"
        result.attrs["unit"] = "MW"
        result.name = f"{result.attrs['name']} ({result.attrs['unit']})"

        if align_edges:
            result = align_edge_directions(result)

        if append_grid:
            result = add_grid_lines(n.df("Bus"), result)

        return result.sort_index()

    def grid_flow(
        self,
        comps: list = None,
        bus_carrier: list = None,
        carrier: list = None,
        aggregate_time: str = "sum",
        append_grid: bool = True,
    ) -> pd.DataFrame:
        """Return the transmission grid energy flow.

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
        n = self.n
        carrier = carrier or list(TRANSMISSION_CARRIER)
        comps = comps or n.branch_components

        energy_transmission = n.statistics.transmission(
            comps=comps,
            groupby=get_buses_and_carrier_and_bus_carrier,
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
        result.name = f"{result.attrs['name']} " f"({result.attrs['unit']})"

        if append_grid:
            result = add_grid_lines(n, result)

        return result.sort_index()


def add_grid_lines(buses: pd.DataFrame, statistic: pd.Series) -> pd.DataFrame:
    """Add a column with gridlines to a statistic.

    Parameters
    ----------
    buses
        The Bus component data frame from a pypsa network.

    statistic
        A pandas object with a multiindex. There must be a "bus0" and
        a "bus1" multiindex level, that hold the node names.

    Returns
    -------
    :
        A data frame with an additional "line" column that holds x/y
        coordinate pairs between the respective bus0 and bus1 locations.
    """
    if isinstance(statistic, pd.Series):
        statistic = statistic.to_frame()

    bus0 = statistic.index.get_level_values("bus0")
    bus1 = statistic.index.get_level_values("bus1")
    ac_buses = filter_by(buses, carrier="AC")[["x", "y"]]

    def _get_bus_lines(_nodes: tuple[str]) -> np.ndarray:
        """Draw a line between buses using AC bus coordinates.

        Note, that only AC buses have coordinates assigned.

        Parameters
        ----------
        _nodes
            The start node name and the end node name in a tuple.

        Returns
        -------
        :
            A one dimensional array with lists of coordinate pairs,
            i.e. grid lines.
        """
        return ac_buses.loc[[*_nodes]][["y", "x"]].to_numpy()

    # generate lines [(x0, y0), (x1,y1)] between buses for every
    # row in grid and store it in a new column
    statistic["line"] = [*map(_get_bus_lines, zip(bus0, bus1, strict=True))]

    return statistic


def align_edge_directions(
    df: pd.DataFrame, lvl0: str = "bus0", lvl1: str = "bus1"
) -> pd.DataFrame:
    """Align the directionality of edges between two nodes.

    Parameters
    ----------
    df
        The input data frame with a multiindex.
    lvl0
        The first MultiIndex level name to swap values.
    lvl1
        The second MultiIndex level name to swap values.

    Returns
    -------
    :
        The input data frame with aligned edge directions between the
        nodes in lvl1 and lvl0.
    """
    seen = []

    def _reverse_values_if_seen(df_slice: pd.DataFrame) -> pd.DataFrame:
        """Reverse index levels if they have a duplicated permutation.

        Parameters
        ----------
        df_slice
            A slice of a data frame with the bus0 and bus1 index level.

        Returns
        -------
        :
            The slice with exchanged level values if the combination of
            lvl1 and lvl2 is not unique and the original slice
            otherwise.
        """
        buses = {df_slice.index.unique(lvl0)[0], df_slice.index.unique(lvl1)[0]}
        if buses in seen:
            reversed_slice = df_slice.swaplevel(lvl0, lvl1)
            # keep original names since we only want to swap values
            reversed_slice.index.names = df_slice.index.names
            return reversed_slice
        else:
            seen.append(buses)
            return df_slice

    return df.groupby([lvl0, lvl1], group_keys=False).apply(
        _reverse_values_if_seen,
    )
