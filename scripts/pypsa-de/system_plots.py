# SPDX-FileCopyrightText: Contributors to PyPSA-DE <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
Per-scenario system plots: storage capacity metrics/maps, energy balances,
electricity capacity supply/demand comparison, and (new, beyond the source
script) capacity-by-region maps for the other capa_groups categories
(Wind + Solar / Backup / Storage Discharge / Demand-Side Flex, see utils.py).

One output folder per `run` scenario (see the `system_plots` rule), covering
that run's full myopic-pathway years. Ported ~1:1 from a sibling repo's
system_plots.py, which iterates `networks[year]` for a single scenario (as
opposed to plot_grid_scenario_comparison.py / that repo's
system_plots_scenario_comparison.py, which compare *across* scenarios).
"""

import logging
import os
import sys
from collections import defaultdict
from pathlib import Path

import cartopy
import cartopy.crs as ccrs
import geopandas as gpd
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pypsa
from matplotlib import colormaps, colors

sys.path.append(os.path.abspath(os.path.dirname(__file__)))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from utils import (
    c1_groups,
    c1_groups_name,
    capa_groups,
    df_to_png,
    plot_grid_network,
    tech_colors,
)

from scripts._helpers import configure_logging, mock_snakemake

logger = logging.getLogger(__name__)

plt.switch_backend("Agg")


def calculate_storage_capacity(n, scenario, year, region="DE", save_plot=True, plot_dir=None):
    """
    Calculate storage capacity metrics for a single network.

    Parameters
    ----------
    n : pypsa.Network
        PyPSA network object
    scenario : str
        Scenario name (e.g., 'LowFlex', 'MediumFlex')
    year : int
        Year of the scenario
    region : str, optional
        Region code to filter (default: "DE")
    save_plot : bool, optional
        Whether to save the result as PNG (default: True)
    plot_dir : str, optional
        Directory to save plots (required if save_plot=True)

    Returns
    -------
    pd.DataFrame
        DataFrame with capacity metrics by technology
    """

    techs = [
        "battery", "home battery", "EV battery",
        "rural water tanks", "urban central water pits",
        "urban central water tanks", "urban decentral water tanks",
        "H2 Store", "PHS", "iron-air battery"
    ]

    result = pd.DataFrame(
        index=techs,
        columns=["energy (GWh)", "discharge (GW)", "charge (GW)"]
    )

    kwargs = {
        "groupby": pypsa.statistics.groupers["bus", "carrier"],
        "nice_names": False,
    }

    # ============= Storage - Stores =============
    stores_capa = (
        n.stores[n.stores.bus.str.contains(region)]
        [["carrier", "e_nom_opt"]]
        .groupby("carrier")
        .sum() / 1e3
    )

    result["energy (GWh)"] = stores_capa.reindex(result.index)["e_nom_opt"].values

    # ============= Storage - Storage Units =============
    su_capa = (
        n.storage_units[n.storage_units.bus.str.contains(region)]
        .assign(energy_capacity=lambda x: x["p_nom_opt"] * x["max_hours"] / 1e3)
        .groupby("carrier")["energy_capacity"]
        .sum()
    )
    if "PHS" in su_capa.index:
        result.loc["PHS", "energy (GWh)"] = su_capa.loc["PHS"]

    # ============= Elec Components =============
    discharge_map = {
        "EV battery": "V2G",
        "iron-air battery": "iron-air battery discharger",
        "battery": "battery discharger",
        "home battery": "home battery discharger",
    }
    charge_map = {
        "EV battery": "BEV charger",
        "iron-air battery": "iron-air battery charger",
        "battery": "battery charger",
        "home battery": "home battery charger",
    }

    elec_capas = (
        n.statistics.optimal_capacity(bus_carrier=["AC", "low voltage"], **kwargs)
        .filter(like=region)
        .groupby("carrier")
        .sum()
        .div(1e3)
    )

    # Default: carrier name == result index (covers PHS, water tanks, H2 Store)
    discharge_labels = [discharge_map.get(t, t) for t in result.index]
    charge_labels = [charge_map.get(t, t) for t in result.index]

    result["discharge (GW)"] = elec_capas.reindex(discharge_labels).values
    result["charge (GW)"] = elec_capas.reindex(charge_labels).values

    # ============= Heat Components =============
    heat_capas = (
        n.statistics.optimal_capacity(
            bus_carrier=["urban central heat", "rural heat", "urban decentral heat"],
            **kwargs
        )
        .filter(like=region)
        .groupby("carrier")
        .sum()
        .div(1e3)
    )

    heat_techs = [
        "rural water tanks", "urban central water pits",
        "urban central water tanks", "urban decentral water tanks",
    ]

    # Only update heat_techs rows - avoids overwriting electricity values with NaNs
    heat_discharge_labels = [f"{t} discharger" for t in heat_techs]
    heat_charge_labels = [f"{t} charger" for t in heat_techs]

    result.loc[heat_techs, "discharge (GW)"] = heat_capas.reindex(heat_discharge_labels).values
    result.loc[heat_techs, "charge (GW)"] = heat_capas.reindex(heat_charge_labels).values

    # ============= H2 Components =============
    h2_capas = (
        n.statistics.optimal_capacity(bus_carrier=["H2"], **kwargs)
        .filter(like=region)
        .drop("Store", errors="ignore")
        .groupby("carrier")
        .sum()
        .div(1e3)
    )

    result.loc["H2 Store", ["discharge (GW)", "charge (GW)"]] = [
        h2_capas.clip(upper=0).sum(),
        h2_capas.clip(lower=0).sum(),
    ]

    # ============= Energy-to-Power Ratio =============
    result["energy-to-power (h)"] = (
        result["energy (GWh)"] / result["discharge (GW)"].where(result["discharge (GW)"] != 0, np.nan)
    )

    # ============= Maximum Usage Statistics =============
    kwargs_stats = {
        "groupby": pypsa.statistics.groupers["name", "bus", "carrier"],
        "nice_names": False,
    }
    buses = n.buses[(n.buses.index.str[:2] == region)].index

    supply = n.statistics.supply(aggregate_time=False, **kwargs_stats)
    demand = n.statistics.withdrawal(aggregate_time=False, **kwargs_stats)

    supply = (
        supply[supply.index.get_level_values("bus").isin(buses)]
        .groupby("carrier")
        .sum()
    )
    demand = (
        demand[demand.index.get_level_values("bus").isin(buses)]
        .groupby("carrier")
        .sum()
    )

    def safe_max(df, carrier):
        return df.loc[carrier].max() if carrier in df.index else np.nan

    discharge_map = {
        "PHS": (supply, "PHS"),
        "EV battery": (supply, "V2G"),
        "iron-air battery": (supply, "iron-air battery discharger"),
        "battery": (supply, "battery discharger"),
        "home battery": (supply, "home battery discharger"),
        **{t: (demand, f"{t} discharger") for t in heat_techs},
    }
    charge_map = {
        "PHS": (demand, "PHS"),
        "EV battery": (demand, "BEV charger"),
        "iron-air battery": (demand, "iron-air battery charger"),
        "battery": (demand, "battery charger"),
        "home battery": (demand, "home battery charger"),
        **{t: (supply, f"{t} charger") for t in heat_techs},
    }

    for key, (df, carrier) in discharge_map.items():
        result.loc[key, "max discharge (GW)"] = safe_max(df, carrier)
    for key, (df, carrier) in charge_map.items():
        result.loc[key, "max charge (GW)"] = safe_max(df, carrier)

    result.loc["H2 Store", "max discharge (GW)"] = supply.reindex(h2_capas[h2_capas < 0].index).sum().max()
    result.loc["H2 Store", "max charge (GW)"] = supply.reindex(h2_capas[h2_capas > 0].index).sum().max()

    result[["max charge (GW)", "max discharge (GW)"]] /= 1e3

    # ============= Save Plot =============
    if save_plot:
        if plot_dir is None:
            raise ValueError("plot_dir must be specified when save_plot=True")
        df = result.astype(float).round(2)
        df_to_png(df, f'{plot_dir}/storage_capacity_{scenario.lower()}_{year}.png')

    return result


def plot_storage_map(
    network,
    technologies,
    onshore_regions,
    stores_capa,
    output_path,
    scenario_name="",
    extent=None,
    figsize=None,
    dpi=300,
    cmap="inferno_r",
    wspace=0.15,
    hspace=0.2,
):
    """
    Plot energy storage capacity maps for different technologies.

    Parameters
    ----------
    network : pypsa.Network
        PyPSA network with bus location data
    technologies : list
        List of storage carrier names to plot
    onshore_regions : gpd.GeoDataFrame
        Geographic regions for background
    stores_capa : pd.DataFrame
        Storage capacities grouped by ["bus", "carrier"] with "e_nom_opt" column
    output_path : str
        Path to save the figure
    scenario_name : str, optional
        Scenario name for title
    extent : list, optional
        Map extent [lon_min, lon_max, lat_min, lat_max]
    figsize : tuple, optional
        Figure size
    dpi : int, optional
        Resolution
    cmap : str, optional
        Colormap name
    wspace : float, optional
        Horizontal spacing between subplots
    hspace : float, optional
        Vertical spacing between subplots
    """
    if extent is None:
        extent = [5.5, 15.5, 47, 56]

    aspect_ratio = (extent[1] - extent[0]) / (extent[3] - extent[2])
    display_projection = ccrs.EqualEarth()

    # Filter available technologies
    available_techs = [t for t in technologies if t in stores_capa.index.get_level_values("carrier").unique()]

    if not available_techs:
        logger.info(f"Skipping storage map for {scenario_name}: no matching technologies")
        return

    # Calculate layout
    n_plots = len(available_techs)
    n_cols = min(3, n_plots)
    n_rows = (n_plots + n_cols - 1) // n_cols

    if figsize is None:
        figsize = (3.5 * n_cols, 5 * n_rows)

    # Create subplots
    fig, axes = plt.subplots(
        n_rows, n_cols,
        subplot_kw={"projection": display_projection},
        figsize=figsize
    )

    axes = np.atleast_1d(axes).flatten()

    for i, tech in enumerate(available_techs):
        ax = axes[i]

        # Follow the working pattern
        df = onshore_regions.copy()
        df = df[df.index.str.contains("DE")]

        capas = stores_capa.xs(tech, level="carrier", drop_level=True)["e_nom_opt"]
        capas.index = network.buses.location.loc[capas.index].values
        capas = capas.groupby(level=0).sum()

        df[tech] = capas.reindex(df.index).fillna(0) / 1e3  # GWh
        total_capacity = df[tech].sum()

        vmin, vmax = df[tech].min(), df[tech].max()

        # Background
        ax.add_feature(cartopy.feature.BORDERS, edgecolor="black", linewidth=0.5)
        ax.coastlines(edgecolor="black", linewidth=0.5)
        ax.set_facecolor("white")
        ax.add_feature(cartopy.feature.OCEAN, color="azure")
        ax.set_title(f"{tech}\nTotal: {total_capacity:.1f} GWh", fontsize=10, pad=15)

        # Plot
        df_plot_crs = df.to_crs(display_projection.proj4_init)
        df_plot_crs.plot(
            column=tech,
            ax=ax,
            linewidth=0.05,
            edgecolor="grey",
            legend=False,
            vmin=vmin,
            vmax=vmax if vmax > 0 else 1,
            cmap=cmap,
        )

        # Annotate GWh value at each region centroid
        centroids = df_plot_crs.geometry.centroid
        for name, value in df[tech].items():
            if value > 0.05:
                pt = centroids.loc[name]
                ax.text(
                    pt.x, pt.y, f"{value:.1f}",
                    transform=display_projection, ha="center", va="center",
                    fontsize=6, color="black",
                    path_effects=[pe.withStroke(linewidth=1.5, foreground="white")],
                )

        ax.set_extent(extent, ccrs.PlateCarree())
        ax.set_aspect(aspect_ratio)

        # Colorbar
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax if vmax > 0 else 1))
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax, shrink=0.6, pad=0.02, orientation="horizontal")
        cbar.set_label("Energy Capacity (GWh)", fontsize=9)

    # Hide unused subplots
    for idx in range(n_plots, len(axes)):
        axes[idx].set_visible(False)

    fig.suptitle(f"Energy Storage Capacity - {scenario_name}", fontsize=14, y=0.98)
    plt.tight_layout()
    plt.subplots_adjust(wspace=wspace, hspace=hspace)
    plt.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def plot_group_capacity_maps(
    network,
    group_name,
    tech_dict,
    onshore_regions,
    output_path,
    scenario_name="",
    region="DE",
    bus_carrier=("AC", "low voltage"),
    extent=None,
    figsize=None,
    dpi=300,
    cmap="viridis",
    wspace=0.15,
    hspace=0.2,
):
    """Map optimal capacity (GW) by region for each tech label in a capa_groups
    category (see utils.py). For "Storage Discharge" this is discharge power
    (its carriers - "battery discharger", "PHS", "V2G", etc. - are already
    discharge-side components), for the other categories it's installed
    capacity - same n.statistics.optimal_capacity call either way.
    """
    if extent is None:
        extent = [5.5, 15.5, 47, 56]

    aspect_ratio = (extent[1] - extent[0]) / (extent[3] - extent[2])
    display_projection = ccrs.EqualEarth()

    caps = (
        network.statistics.optimal_capacity(
            bus_carrier=list(bus_carrier),
            groupby=pypsa.statistics.groupers["bus", "carrier"],
            nice_names=False,
        )
        .filter(like=region)
        .groupby(["bus", "carrier"])
        .sum()
        .div(1e3)  # GW
    )

    per_tech = {}
    for tech_label, carriers in tech_dict.items():
        sel = caps[caps.index.get_level_values("carrier").isin(carriers)]
        if sel.empty:
            continue
        by_bus = sel.groupby("bus").sum()
        by_bus.index = network.buses.location.loc[by_bus.index].values
        by_bus = by_bus.groupby(level=0).sum()
        if by_bus.abs().sum() > 0.05:  # skip near-zero techs (GW)
            per_tech[tech_label] = by_bus

    if not per_tech:
        logger.info(f"Skipping capacity map for '{group_name}' ({scenario_name}): no data")
        return

    n_plots = len(per_tech)
    n_cols = min(3, n_plots)
    n_rows = (n_plots + n_cols - 1) // n_cols
    if figsize is None:
        figsize = (3.5 * n_cols, 5 * n_rows)

    fig, axes = plt.subplots(
        n_rows, n_cols, subplot_kw={"projection": display_projection}, figsize=figsize
    )
    axes = np.atleast_1d(axes).flatten()

    for i, tech_label in enumerate(per_tech):
        ax = axes[i]

        df = onshore_regions.copy()
        df = df[df.index.str.contains(region)]
        df[tech_label] = per_tech[tech_label].reindex(df.index).fillna(0)
        total = df[tech_label].sum()

        vmin, vmax = df[tech_label].min(), df[tech_label].max()

        ax.add_feature(cartopy.feature.BORDERS, edgecolor="black", linewidth=0.5)
        ax.coastlines(edgecolor="black", linewidth=0.5)
        ax.set_facecolor("white")
        ax.add_feature(cartopy.feature.OCEAN, color="azure")
        ax.set_title(f"{tech_label}\nTotal: {total:.1f} GW", fontsize=10, pad=15)

        df_plot_crs = df.to_crs(display_projection.proj4_init)
        df_plot_crs.plot(
            column=tech_label, ax=ax, linewidth=0.05, edgecolor="grey",
            legend=False, vmin=vmin, vmax=vmax if vmax > 0 else 1, cmap=cmap,
        )

        centroids = df_plot_crs.geometry.centroid
        for name, value in df[tech_label].items():
            if value > 0.05:
                pt = centroids.loc[name]
                ax.text(
                    pt.x, pt.y, f"{value:.1f}",
                    transform=display_projection, ha="center", va="center",
                    fontsize=6, color="black",
                    path_effects=[pe.withStroke(linewidth=1.5, foreground="white")],
                )

        ax.set_extent(extent, ccrs.PlateCarree())
        ax.set_aspect(aspect_ratio)

        sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax if vmax > 0 else 1))
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax, shrink=0.6, pad=0.02, orientation="horizontal")
        cbar.set_label("Capacity (GW)", fontsize=9)

    for idx in range(n_plots, len(axes)):
        axes[idx].set_visible(False)

    fig.suptitle(f"{group_name} Capacity - {scenario_name}", fontsize=14, y=0.98)
    plt.tight_layout()
    plt.subplots_adjust(wspace=wspace, hspace=hspace)
    plt.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def power_capacity_by_bus_carrier(network, region="DE", bus_carrier=("AC", "low voltage")):
    """(bus, carrier)-indexed optimal power capacity (GW) - the same data
    plot_group_capacity_maps plots, factored out so diff maps (see
    group_capacity_by_tech/plot_group_capacity_diff_maps below) can compute
    it for two networks and subtract."""
    return (
        network.statistics.optimal_capacity(
            bus_carrier=list(bus_carrier),
            groupby=pypsa.statistics.groupers["bus", "carrier"],
            nice_names=False,
        )
        .filter(like=region)
        .groupby(["bus", "carrier"])
        .sum()
        .div(1e3)  # GW
    )


def storage_energy_capacity_by_bus_carrier(network, region="DE"):
    """(bus, carrier)-indexed energy capacity (GWh) for Stores plus PHS
    storage_units (implied e_nom_opt = p_nom_opt * max_hours) - the same
    data plot_storage_map plots per-carrier; grouped here into
    utils.storage_energy_capa_groups' {tech_label: [carriers]} shape
    instead so it can go through the same group-map machinery as the
    power-capacity capa_groups categories.
    """
    stores_capa = (
        network.stores[network.stores.bus.str.contains(region)][["bus", "carrier", "e_nom_opt"]]
        .groupby(["bus", "carrier"])
        .sum()["e_nom_opt"]
        / 1e3  # GWh
    )
    su = network.storage_units[
        network.storage_units.bus.str.contains(region) & (network.storage_units.carrier == "PHS")
    ]
    su_capa = (
        su.assign(e_nom_opt=lambda x: x["p_nom_opt"] * x["max_hours"] / 1e3)
        .groupby(["bus", "carrier"])["e_nom_opt"]
        .sum()
    )
    return pd.concat([stores_capa, su_capa])


def group_capacity_by_tech(network, tech_dict, capacity_series):
    """Group a (bus, carrier)-indexed capacity Series (see
    power_capacity_by_bus_carrier/storage_energy_capacity_by_bus_carrier)
    into per-tech-label, per-region (bus location) totals, matching
    capa_groups'/storage_energy_capa_groups' {tech_label: [carriers]}
    structure. Shared data-prep step behind plot_group_capacity_diff_maps.

    Every tech_label in tech_dict gets an entry (possibly all-zero/empty),
    unlike plot_group_capacity_maps' own near-zero-total cutoff - diff maps
    need the full, consistent carrier set so a tech that's genuinely
    identical between the two networks being compared (e.g. offshore wind
    built out to its full potential either way) still gets its subplot
    rather than silently disappearing.
    """
    per_tech = {}
    for tech_label, carriers in tech_dict.items():
        sel = capacity_series[capacity_series.index.get_level_values("carrier").isin(carriers)]
        if sel.empty:
            per_tech[tech_label] = pd.Series(dtype=float)
            continue
        by_bus = sel.groupby("bus").sum()
        by_bus.index = network.buses.location.loc[by_bus.index].values
        per_tech[tech_label] = by_bus.groupby(level=0).sum()
    return per_tech


def plot_group_capacity_diff_maps(
    per_tech_a,
    per_tech_b,
    group_name,
    onshore_regions,
    output_path,
    label_a,
    label_b,
    region="DE",
    unit="GW",
    extent=None,
    figsize=None,
    dpi=300,
    cmap="RdBu_r",
    wspace=0.15,
    hspace=0.2,
):
    """Map the regional capacity difference (label_a minus label_b) for each
    tech label in a capacity-group category - diverging colormap centered
    on zero, same annotation style as plot_group_capacity_maps.

    per_tech_a/per_tech_b : dict[tech_label] -> pd.Series indexed by region
        (bus location), as returned by group_capacity_by_tech, for the two
        networks being compared.
    """
    if extent is None:
        extent = [5.5, 15.5, 47, 56]
    aspect_ratio = (extent[1] - extent[0]) / (extent[3] - extent[2])
    display_projection = ccrs.EqualEarth()

    # Same tech_label set for both (group_capacity_by_tech always returns
    # one entry per tech_dict key), kept in tech_dict's original order.
    tech_labels = list(per_tech_a) + [t for t in per_tech_b if t not in per_tech_a]
    diffs = {
        tech_label: per_tech_a.get(tech_label, pd.Series(dtype=float)).subtract(
            per_tech_b.get(tech_label, pd.Series(dtype=float)), fill_value=0
        )
        for tech_label in tech_labels
    }

    if not diffs:
        logger.info(f"Skipping capacity diff map for '{group_name}' ({label_a} vs {label_b}): empty group")
        return

    n_plots = len(diffs)
    n_cols = min(3, n_plots)
    n_rows = (n_plots + n_cols - 1) // n_cols
    if figsize is None:
        figsize = (3.5 * n_cols, 5 * n_rows)

    fig, axes = plt.subplots(
        n_rows, n_cols, subplot_kw={"projection": display_projection}, figsize=figsize
    )
    axes = np.atleast_1d(axes).flatten()

    for i, tech_label in enumerate(diffs):
        ax = axes[i]

        df = onshore_regions.copy()
        df = df[df.index.str.contains(region)]
        df[tech_label] = diffs[tech_label].reindex(df.index).fillna(0)
        total = df[tech_label].sum()

        vabs = max(df[tech_label].abs().max(), 0.05)
        vmin, vmax = -vabs, vabs

        ax.add_feature(cartopy.feature.BORDERS, edgecolor="black", linewidth=0.5)
        ax.coastlines(edgecolor="black", linewidth=0.5)
        ax.set_facecolor("white")
        ax.add_feature(cartopy.feature.OCEAN, color="azure")
        ax.set_title(f"{tech_label}\nTotal Δ: {total:+.1f} {unit}", fontsize=10, pad=15)

        df_plot_crs = df.to_crs(display_projection.proj4_init)
        df_plot_crs.plot(
            column=tech_label, ax=ax, linewidth=0.05, edgecolor="grey",
            legend=False, vmin=vmin, vmax=vmax, cmap=cmap,
        )

        centroids = df_plot_crs.geometry.centroid
        for name, value in df[tech_label].items():
            if abs(value) > 0.05:
                pt = centroids.loc[name]
                ax.text(
                    pt.x, pt.y, f"{value:+.1f}",
                    transform=display_projection, ha="center", va="center",
                    fontsize=6, color="black",
                    path_effects=[pe.withStroke(linewidth=1.5, foreground="white")],
                )

        ax.set_extent(extent, ccrs.PlateCarree())
        ax.set_aspect(aspect_ratio)

        sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax, shrink=0.6, pad=0.02, orientation="horizontal")
        cbar.set_label(f"Δ Capacity ({unit})", fontsize=9)

    for idx in range(n_plots, len(axes)):
        axes[idx].set_visible(False)

    fig.suptitle(f"{group_name} Capacity Diff: {label_a} − {label_b}", fontsize=14, y=0.98)
    plt.tight_layout()
    # Explicit top margin (on top of tight_layout's own spacing) - every
    # tech_label is now always plotted (see group_capacity_by_tech), so the
    # grid's last row is often partially filled; tight_layout alone doesn't
    # reserve room for the suptitle in that case and it overlaps the first
    # row's subplot titles.
    plt.subplots_adjust(top=0.88, wspace=wspace, hspace=hspace)
    plt.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def plot_grid_topology_maps(
    n_prenetwork, n_solved, onshore_regions, output_dir, scenario, year
):
    """AC/DC grid-topology maps (utils.plot_grid_network - see
    notebooks/grid_topology.ipynb), for both the "de_internal" and
    "de_and_interconnectors" scopes, for the network as handed to the
    solver ("prenetwork" - whatever grid CSV override is already applied,
    but before capacity/dispatch optimization) and after solving - 8 maps
    total. `n_prenetwork` may be None (skips that stage) - see
    _grid_prenetwork_path in stochastic_grid.smk for cases with no clean
    prenetwork counterpart.
    """
    stages = [("prenetwork", n_prenetwork), ("solved", n_solved)]
    scopes = ["de_internal", "de_and_interconnectors"]
    components = {"Line": "AC", "Link": "DC"}

    for stage_name, n in stages:
        if n is None:
            continue
        for scope in scopes:
            for component, comp_label in components.items():
                fig, ax = plot_grid_network(
                    n,
                    components=component,
                    scope=scope,
                    horizon=year,
                    annotate=True,
                    figsize=(8, 8),
                    regions=onshore_regions,
                    title=f"{scenario} {stage_name} ({year}) - {comp_label}, {scope}",
                    export_path=output_dir
                    / f"grid_topology_{stage_name}_{scope}_{comp_label}_{scenario}_{year}.png",
                )
                plt.close(fig)


def plot_balance(nb, title="title", tech_colors=None, save_path=None):
    import matplotlib.colors as mcolors

    resample = "D"
    nb = nb.resample(resample).mean()
    df = nb

    # Check if dataframe is empty or has no numeric data
    if df.empty:
        logger.info(f"WARNING: Skipping plot '{title}' - no data available")
        return

    # Check if all values are zero or NaN
    if df.abs().sum().sum() == 0 or df.isna().all().all():
        logger.info(f"WARNING: Skipping plot '{title}' - all values are zero or NaN")
        return

    # split into df with positive and negative values
    df_neg, df_pos = df.clip(upper=0), df.clip(lower=0)
    df_pos = df_pos[df_pos.sum().sort_values(ascending=False).index]
    df_neg = df_neg[df_neg.sum().sort_values().index]

    # Check if both are empty after splitting
    if df_pos.empty and df_neg.empty:
        logger.info(f"WARNING: Skipping plot '{title}' - no positive or negative values to plot")
        return

    # Check if we have any non-zero columns
    df_pos = df_pos.loc[:, (df_pos != 0).any(axis=0)]
    df_neg = df_neg.loc[:, (df_neg != 0).any(axis=0)]

    if df_pos.empty and df_neg.empty:
        logger.info(f"WARNING: Skipping plot '{title}' - all columns contain only zeros")
        return

    # Helper function to validate and fix colors
    def validate_color(col, color_value):
        if color_value is None:
            logger.info(f"WARNING: Missing color for carrier '{col}', using 'grey'")
            return 'grey'
        try:
            # Try to convert to rgba to validate
            mcolors.to_rgba(color_value)
            return color_value
        except (ValueError, TypeError) as e:
            logger.info(f"WARNING: Invalid color for carrier '{col}': {repr(color_value)} (type: {type(color_value).__name__})")
            logger.info(f"         Error: {e}")
            return 'grey'

    # get colors with validation
    c_neg = [
        validate_color(col, tech_colors.get(col, 'grey') if tech_colors else 'grey')
        for col in df_neg.columns
    ]
    c_pos = [
        validate_color(col, tech_colors.get(col, 'grey') if tech_colors else 'grey')
        for col in df_pos.columns
    ]

    fig, ax = plt.subplots(figsize=(14, 8))

    # plot positive values (only if not empty)
    if not df_pos.empty:
        ax = df_pos.plot.area(ax=ax, stacked=True, color=c_pos, linewidth=0.0)

    # rename negative values that are also present on positive side, so that they are not shown and plot negative values
    if not df_neg.empty:
        def f(c):
            return "out_" + c

        cols = [f(c) if (c in df_pos.columns) else c for c in df_neg.columns]
        cols_map = dict(zip(df_neg.columns, cols))
        ax = df_neg.rename(columns=cols_map).plot.area(
            ax=ax, stacked=True, color=c_neg, linewidth=0.0
        )

    # explicitly filter out duplicate labels
    handles, labels = ax.get_legend_handles_labels()

    if not handles:
        logger.info(f"WARNING: Skipping plot '{title}' - no data to display in legend")
        plt.close(fig)
        return

    filtered_handles_labels = [
        (h, l) for h, l in zip(handles, labels) if not l.startswith("out_")
    ]

    if not filtered_handles_labels:
        logger.info(f"WARNING: Skipping plot '{title}' - no valid legend entries")
        plt.close(fig)
        return

    handles, labels = zip(*filtered_handles_labels)

    # rescale the y-axis
    y_min = df_neg.sum(axis=1).min() if not df_neg.empty else 0
    y_max = df_pos.sum(axis=1).max() if not df_pos.empty else 0

    if y_min == 0 and y_max == 0:
        logger.info(f"WARNING: Skipping plot '{title}' - no variation in data")
        plt.close(fig)
        return

    ax.set_ylim([1.05 * y_min, 1.05 * y_max])
    ax.legend(
        handles,
        labels,
        ncol=1,
        loc="upper center",
        bbox_to_anchor=(1.13, 1.01),
    )
    ax.set_title(title)
    ax.grid(True)

    if save_path is not None:
        plt.savefig(save_path, bbox_inches="tight", dpi=300)
        logger.info(f"Saved plot: {save_path}")

    plt.close(fig)


def plot_electricity_capacity(
    networks,
    years,
    c1_groups,
    c1_groups_name,
    tech_colors,
    CAPACITY_DIR,
    *,
    de_only=True,
    scenario_name="scenario",
    ylim=None,
    kwargs=None,
):
    """
    Plot electricity capacity for one scenario.

    Parameters
    ----------
    networks : dict[int, pypsa.Network]
        Dictionary like networks[year]
    years : list[int]
        Years to plot
    c1_groups, c1_groups_name : list
        Technology aggregation definition
    tech_colors : dict
        Mapping tech -> color
    CAPACITY_DIR : str or Path
        Output directory
    de_only : bool, default True
        If True: only DE capacities (bus_carrier=["AC"], filter "DE")
        If False: whole system
    scenario_name : str
        Used in output filename
    ylim : float or None
        y-axis limit (auto if None)
    kwargs : dict or None
        Passed to statistics.optimal_capacity
    """

    if kwargs is None:
        kwargs = {}

    # -----------------------
    # Capacity extraction
    # -----------------------
    cap_all = pd.DataFrame()

    for year in years:
        n = networks[year]

        if de_only:
            cap = (
                n.statistics.optimal_capacity(bus_carrier=["AC"], **kwargs)
                .filter(like="DE")
                .groupby("carrier")
                .sum()
                .div(1e3)  # MW -> GW
                .to_frame(name=year)
            )
        else:
            cap = (
                n.statistics.optimal_capacity(bus_carrier=["AC"], nice_names=False)
                .div(1e3)  # MW -> GW
                .droplevel(0)
                .to_frame(name=year)
            )

        cap_all = cap_all.combine_first(cap) if not cap_all.empty else cap

    # Drop non-technologies
    cap_all = cap_all.drop(["load", "load-shedding"], errors="ignore")

    # -----------------------
    # Aggregate technologies
    # -----------------------
    c1_groups = c1_groups.copy()
    c1_groups_name = c1_groups_name.copy()

    if "solar" in c1_groups_name:
        i = c1_groups_name.index("solar")
        del c1_groups[i]
        del c1_groups_name[i]

    df_new = cap_all.copy()
    all_grouped_rows = []

    for group, name in zip(c1_groups, c1_groups_name):
        existing = [d for d in group if d in df_new.index]
        if existing:
            df_new.loc[name] = df_new.loc[existing].sum()
            all_grouped_rows.extend([d for d in existing if d != name])

    cap_all_agg = df_new.drop(index=all_grouped_rows)

    # -----------------------
    # Select top technologies
    # -----------------------
    techs = (
        cap_all_agg.abs()
        .sum(axis=1)
        .sort_values(ascending=False)
        .head(16)
        .index
    )

    df_plot = cap_all_agg.loc[techs].fillna(0)

    sort_year = years[-1]
    sorted_techs = (
        df_plot[sort_year].abs().sort_values(ascending=False).index.tolist()
    )

    # -----------------------
    # Plot
    # -----------------------
    n_years = len(years)

    fig, axes = plt.subplots(
        1,
        n_years,
        figsize=(6 * n_years, 5),
        sharex=False,
        sharey=True,
    )

    # Make axes always iterable
    if n_years == 1:
        axes = [axes]

    for i, year in enumerate(years):
        ax = axes[i]

        year_data = df_plot[year].loc[sorted_techs]

        supply = year_data[year_data >= 0]
        demand = -year_data[year_data < 0]

        pos_supply = np.arange(len(supply))
        pos_demand = np.arange(len(demand)) + len(supply) + 1

        if ylim is not None:
            ax.set_ylim(0, ylim)

        ax.bar(
            pos_supply,
            supply.values,
            color=[tech_colors.get(t, "#1f77b4") for t in supply.index],
            alpha=0.8,
        )
        ax.bar(
            pos_demand,
            demand.values,
            color=[tech_colors.get(t, "#d62728") for t in demand.index],
            alpha=0.8,
        )

        ax.axvline(len(supply) - 0.5, color="black", linestyle="--", linewidth=1.0)

        ax.text(0.02, 0.95, "Supply", transform=ax.transAxes,
                ha="left", va="top", fontsize=10, fontweight="bold")
        ax.text(0.98, 0.95, "Demand", transform=ax.transAxes,
                ha="right", va="top", fontsize=10, fontweight="bold")

        ymax = max(supply.max() if not supply.empty else 0,
                   demand.max() if not demand.empty else 0)

        for x, y in zip(pos_supply, supply.values):
            ax.text(x, y + 0.01 * ymax, f"{y:.0f}",
                    ha="center", va="bottom", fontsize=6)
        for x, y in zip(pos_demand, demand.values):
            ax.text(x, y + 0.01 * ymax, f"{y:.0f}",
                    ha="center", va="bottom", fontsize=6)

        tech_labels = list(supply.index) + list(demand.index)
        ax.set_xticks(list(pos_supply) + list(pos_demand))
        ax.set_xticklabels(tech_labels, rotation=45, ha="right", fontsize=8)

        ax.set_title(f"{year}", fontsize=16)
        ax.set_ylabel("GW", fontsize=12)
        ax.grid(True, axis="y", alpha=0.3)
        ax.axhline(0, color="black", linewidth=0.8)

    plt.tight_layout()

    scope = "de" if de_only else "all-system"
    plt.savefig(
        f"{CAPACITY_DIR}/electricity-capacity-{scope}-{scenario_name}-{'-'.join(map(str, years))}.png",
        bbox_inches="tight",
    )
    plt.close(fig)


if __name__ == "__main__":
    if "snakemake" not in globals():
        snakemake = mock_snakemake(
            "system_plots",
            clusters=27,
            opts="",
            sector_opts="none",
            run="KN2045_Mix",
        )

    configure_logging(snakemake)
    planning_horizons = list(snakemake.params.planning_horizons)
    # For the stochastic-grid variant of this rule (system_plots_grid_scenario
    # in stochastic_grid.smk), `run` is identical across all network_ids (see
    # export_ariadne_variables_grid_scenario for the same issue/fix) - label
    # by network_id there instead so titles/filenames stay distinguishable.
    scenario = (
        snakemake.wildcards.network_id
        if "network_id" in snakemake.wildcards.keys()
        else snakemake.wildcards.run
    )

    output_dir = Path(snakemake.params.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    PLOT_DIR = output_dir

    # Keyed by year (not parsed from the filename tail): input.networks is
    # expanded over config["scenario"]["planning_horizons"] in that same
    # order, so zipping is robust to whatever the network filenames look like.
    networks = {y: pypsa.Network(fn) for y, fn in zip(planning_horizons, snakemake.input.networks)}
    # Pre-solve counterpart of `networks` (grid CSV override already
    # applied, if any, but before capacity/dispatch optimization) - see
    # plot_grid_topology_maps/_grid_prenetwork_path (stochastic_grid.smk).
    prenetworks = {
        y: pypsa.Network(fn) for y, fn in zip(planning_horizons, snakemake.input.prenetworks)
    }

    ####### PLOTTING #########
    kwargs = {
        "groupby": pypsa.statistics.groupers["bus", "carrier"],
        "nice_names": False,
    }

    ### Storage Capacity Metrics ###
    for year in planning_horizons:
        n = networks[year]
        calculate_storage_capacity(
            n,
            scenario=scenario,
            year=year,
            plot_dir=PLOT_DIR
        )

    ### Capacities Map ###

    base = colormaps.get_cmap("RdPu")
    cmap = colors.LinearSegmentedColormap.from_list("", base(np.linspace(0.15, 1, 256)))
    onshore_regions = gpd.read_file(snakemake.input.regions_onshore[0]).set_index("name")

    ### Grid Topology Maps (prenetwork vs solved, AC/DC, DE-internal/+interconnectors) ###
    for year in planning_horizons:
        plot_grid_topology_maps(
            prenetworks.get(year), networks[year], onshore_regions, PLOT_DIR, scenario, year
        )

    for year in planning_horizons:
        n = networks[year]
        stores_capa = n.stores[n.stores.bus.str.contains("DE")][["bus", "carrier", "e_nom_opt"]].groupby(["bus", "carrier"]).sum() / 1e3

        su_capa = (n.storage_units[n.storage_units.bus.str.contains("DE") & (n.storage_units.carrier == "PHS")]
                .assign(e_nom_opt=lambda x: x["p_nom_opt"] * x["max_hours"] / 1e3)
                .groupby(["bus", "carrier"])["e_nom_opt"].sum())

        stores_capa = pd.concat([stores_capa, su_capa])

        techs = ["PHS",
                "battery",
                "home battery",
                "EV battery",
                "rural water tanks",
                "urban central water pits",
                "urban central water tanks",
                "urban decentral water tanks",
                "H2 Store",
                "iron-air battery",  # dropped automatically if not present in this network
        ]

        plot_storage_map(
            n,
            techs,
            onshore_regions,
            stores_capa,
            output_path=f"{PLOT_DIR}/storage_map_{scenario}_{year}.png",
            scenario_name=scenario,
            cmap=cmap,
        )

        # Capacity-by-region maps for the other capa_groups categories (Wind +
        # Solar / Backup / Storage Discharge / Demand-Side Flex) - beyond the
        # storage-only maps above.
        for group_name, tech_dict in capa_groups.items():
            group_slug = group_name.lower().replace(" + ", "_").replace(" ", "_")
            plot_group_capacity_maps(
                n,
                group_name,
                tech_dict,
                onshore_regions,
                output_path=f"{PLOT_DIR}/capacity_map_{group_slug}_{scenario}_{year}.png",
                scenario_name=scenario,
            )

    ### ENERGY BALANCES ###

    balances = defaultdict(dict)
    for year in planning_horizons:
        network = networks[year]
        ct = "DE"
        buses = network.buses.index[(network.buses.index.str[:2] == ct)].drop("DE")
        balance = (
            network.statistics.energy_balance(
                aggregate_time=False,
                nice_names=False,
                groupby=["bus", "carrier", "bus_carrier"],
            )
            .loc[:, buses, :, :]
            .droplevel("bus")
        )
        balances[scenario][year] = balance

    carriers_sets = [["AC", "low voltage"],
                     ["urban central heat", "rural heat", "urban decentral heat"],
                     ["H2"],
                     ["co2 stored"],
                     ["oil"],
                     ["renewable oil"],
                     ["gas"],
                     ["renewables gas"],
                     ["solid biomass"],
                     ]

    for year in planning_horizons:

        for carriers in carriers_sets:
            mask = balances[scenario][year].index.get_level_values("bus_carrier").isin(carriers)
            nb = balances[scenario][year][mask].groupby("carrier").sum().div(1e3).T

            plot_balance(
                nb,
                title=f"Energy Balance of - '{', '.join(carriers)}' ({scenario}, {year})",
                tech_colors=tech_colors,
                save_path=f"{PLOT_DIR}/energy_balance_{'-'.join(carriers)}_{scenario}_{year}.png",
            )

    ### CAPACITY COMPARISON ###

    # DE only
    plot_electricity_capacity(
        networks=networks,
        years=planning_horizons,
        c1_groups=c1_groups,
        c1_groups_name=c1_groups_name,
        tech_colors=tech_colors,
        CAPACITY_DIR=PLOT_DIR,
        de_only=True,
        scenario_name=scenario,
        ylim=500,
        kwargs=kwargs,
    )

    # Whole system
    plot_electricity_capacity(
        networks=networks,
        years=planning_horizons,
        c1_groups=c1_groups,
        c1_groups_name=c1_groups_name,
        tech_colors=tech_colors,
        CAPACITY_DIR=PLOT_DIR,
        de_only=False,
        scenario_name=scenario,
        ylim=1500,
    )

    logger.info(f"All system plots for {scenario} saved to: {PLOT_DIR}")
