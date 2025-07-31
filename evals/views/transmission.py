# SPDX-FileCopyrightText: 2023-2025 Austrian Gas Grid Management AG
#
# SPDX-License-Identifier: MIT
# For license information, see the LICENSE.txt file in the project root.
"""Module for trade evaluations."""

from pathlib import Path

from evals.constants import BusCarrier, Carrier, DataModel
from evals.plots.gridmap import GridMapConfig, TransmissionGridMap
from evals.statistic import collect_myopic_statistics
from evals.utils import filter_by

"""
Todo Notes
* click country label
  -> Table energy amounts import + export
* global import mit icon
* global import label click -> show table
* global import edge to bus node
  -> generator capacity is line width
* display extended capacity
* labels grid must not overlap
* solid waste transport
* biomass transport

"""


def view_grid_capacity(
    result_path: str | Path,
    networks: dict,
    config: dict,
    subdir: str | Path = "evaluation",
) -> None:  # numpydoc ignore=PR01
    """Export transmission grids to file using Folium."""

    # update bus coordinates to improve map readability
    for n in networks.values():
        # Lower Austria
        n.df("Bus").loc["AT12", "x"] -= 0.1  # Lon, left
        n.df("Bus").loc["AT12", "y"] += 0.3  # Lat, up
        # Burgenland
        n.df("Bus").loc["AT11", "x"] -= 0.1  # Lon, left
        n.df("Bus").loc["AT11", "y"] -= 0.26  # Lat, down
        # Salzburg
        n.df("Bus").loc["AT32", "x"] += 0.35  # Lon, right
        n.df("Bus").loc["AT32", "y"] -= 0.05  # Lat, down
        # Vienna
        n.df("Bus").loc["AT13", "x"] += 0.1  # Lon, right
        n.df("Bus").loc["AT13", "y"] += 0.1  # Lat, up

    grid_capactiy = collect_myopic_statistics(
        networks,
        statistic="grid_capacity",
        drop_zeros=False,
        comps=["Link", "Line"],
    )

    # cannot use utils.scale(), because of the additional "line" column
    col = "Capacity (MW)"
    grid_capactiy[col] = grid_capactiy[col] * 1e-3
    grid_capactiy = grid_capactiy.rename(columns={col: "Capacity (GW)"})
    grid_capactiy.attrs["unit"] = "GW"

    import_energy = collect_myopic_statistics(
        networks,
        statistic="supply",
        comps="Generator",
        bus_carrier=[BusCarrier.CH4, BusCarrier.H2, "biogas", "AC"],  # "gas primary",
    )
    import_energy *= 1e-6
    import_energy.attrs["name"] = "Import Energy"
    import_energy.attrs["unit"] = "TWh"
    metric_name = f"{import_energy.attrs['name']} ({import_energy.attrs['unit']})"
    import_energy.name = metric_name

    # the optimal capacity for pipelines is larger than the maximal
    # energy flow in the time series, because pipelines are oversized.
    # We use the maximal flow here since it is more interesting.
    import_capacity = collect_myopic_statistics(
        networks,
        statistic="supply",
        comps="Generator",
        bus_carrier=[BusCarrier.CH4, BusCarrier.H2, "biogas"],  # "gas primary"
        aggregate_time="max",
    ).drop("2015", level=DataModel.YEAR, errors="ignore")
    import_capacity *= 1e-3
    import_capacity.attrs["name"] = "Import Capacity"
    import_capacity.attrs["unit"] = "GW"
    metric_name = f"{import_capacity.attrs['name']} ({import_capacity.attrs['unit']})"
    import_capacity.name = metric_name

    config = GridMapConfig(show_year="2030")  # fixme: show_year is broken =(
    buses = networks[next(reversed(networks))].df("Bus")

    # every list item will become one HTML file with a map for the
    # specified carrier and bus_carrier
    # ToDo: Add CO2 once it works properly
    carriers_bus_carrier_groups = (
        ([Carrier.AC, Carrier.DC], BusCarrier.AC),
        ([Carrier.gas_pipepline, Carrier.gas_pipepline_new], BusCarrier.CH4),
        (
            [  # todo: use get_transmission_techs() instead of hardcoding
                Carrier.h2_pipeline,
                Carrier.h2_pipeline_retro,
                Carrier.h2_pipeline_kernnetz,
            ],
            BusCarrier.H2,
        ),
    )

    for carriers, bus_carrier in carriers_bus_carrier_groups:
        df_grid = filter_by(grid_capactiy, carrier=carriers)
        df_import_energy = filter_by(import_energy, bus_carrier=bus_carrier)
        df_import_energy = df_import_energy[df_import_energy > 0]
        df_import_capacity = filter_by(import_capacity, bus_carrier=bus_carrier)
        df_import_capacity = df_import_capacity[df_import_capacity > 0]
        grid_map = TransmissionGridMap(
            df_grid, df_import_energy, df_import_capacity, buses, config
        )
        grid_map.draw_grid_by_carrier_groups_myopic()
        grid_map.save(result_path, f"gridmap_{bus_carrier}", subdir)
