# SPDX-FileCopyrightText: Contributors to PyPSA-DE <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
General-purpose helpers for ad-hoc PyPSA-DE analysis notebooks/scripts - not
part of the main pipeline DAG, and not tied to any one notebook. The idea is
a shared landing spot for reusable notebook helpers as they get factored out,
instead of every notebook re-implementing its own copy.

Home to grid-topology plotting (`plot_grid_network`, factored out of
notebooks/grid_topology.ipynb) and flexibility/scenario-comparison plotting
helpers (tech_colors, sector_colors, df_to_png, aggregate_by_keywords, etc.,
ported as-is from a sibling repo's flexibility_utils.py/flexibility_analysis.py
for use by plot_grid_scenario_comparison.py); more sections can be added
alongside them as other notebooks/scripts contribute their own reusable
pieces.

Imported interactively, e.g.:

    import sys
    sys.path.append(str(REPO / "scripts" / "pypsa-de"))
    from utils import plot_grid_network

    plot_grid_network(n, components="Line", scope="de_and_interconnectors")
"""

import colorsys
import hashlib
import re

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import numpy as np
import pandas as pd
from numpy import isclose

# ---------------------------------------------------------------------------
# Grid-topology plotting
# ---------------------------------------------------------------------------

NOM_ATTR = {"Line": "s_nom", "Link": "p_nom"}

# Static map extents in [lon_min, lon_max, lat_min, lat_max] - deliberately
# fixed rather than fit to the data, so a given scope always frames the same
# geography regardless of which branches happen to be in the plotted table.
# "de_and_interconnectors" is kept nearly as tight as "de_internal" - just
# enough margin to show interconnectors leaving the frame, not a Europe-wide
# view - since annotations that would fall outside it get clipped to the
# edge anyway (see `_annotate_branches`).
EXTENT = {
    "whole_network": [-11, 30, 34, 65],
    "de_and_interconnectors": [3, 17, 46, 57],
    "de_internal": [4, 16, 47, 56],
}

COMPONENT_COLORS = {"Line": "steelblue", "Link": "darkorange"}
COMPONENT_LABELS = {"Line": "Line (AC)", "Link": "Link (DC)"}


def _branches(n, component):
    return n.lines if component == "Line" else n.links.query("carrier == 'DC'")


def branch_table(n, component, horizon=None):
    '''One row per branch: effective capacity (opt if solved, else nom), length, build_year, endpoints.

    `horizon` only matters for an *unsolved* Link: capacity is floored to 0
    for any branch whose `build_year` is still in the future relative to
    `horizon` (mirrors `enforce_transmission_project_build_years` in
    `modify_prenetwork.py`, which gates Links only, never Lines). Pass
    `horizon=None` (the default) to skip that gating entirely and use the
    network's own capacity as-is - the right choice unless you specifically
    need "what's built by year X".
    '''
    nom = NOM_ATTR[component]
    df = _branches(n, component)
    opt = df[f"{nom}_opt"]
    if (opt > 0).any():
        capacity = opt
    elif component == "Link" and horizon is not None:
        not_yet_due = df.build_year > horizon
        capacity = df[nom].where(~not_yet_due, 0.0)
    else:
        capacity = df[nom]
    return pd.DataFrame(
        {
            "capacity_mw": capacity,
            "length_km": df["length"],
            "build_year": df["build_year"],
            "bus0": df["bus0"],
            "bus1": df["bus1"],
            "country0": n.buses.loc[df.bus0, "country"].to_numpy(),
            "country1": n.buses.loc[df.bus1, "country"].to_numpy(),
        },
        index=df.index,
    )


def scope_mask(table, scope, home_country="DE"):
    '''Which branches a scope keeps, based on `home_country` bus membership:
    "whole_network" - everything; "de_internal" - both ends in `home_country`;
    "de_and_interconnectors" - at least one end in `home_country`.'''
    c0, c1 = table.country0 == home_country, table.country1 == home_country
    if scope == "de_internal":
        return c0 & c1
    if scope == "de_and_interconnectors":
        return c0 | c1
    if scope == "whole_network":
        return pd.Series(True, index=table.index)
    raise ValueError(f"unknown scope {scope!r}, expected one of {list(EXTENT)}")


def _basemap(ax, extent, regions=None):
    ax.add_feature(cfeature.OCEAN, facecolor="#d5e6f0", zorder=0)
    ax.add_feature(cfeature.LAND, facecolor="#f7f4ee", zorder=0)
    if regions is not None:
        regions.to_crs(epsg=4326).plot(ax=ax, transform=ccrs.PlateCarree(), facecolor="#eae6da",
                                        edgecolor="white", linewidth=0.5, zorder=0.5)
    ax.coastlines(resolution="50m", linewidth=0.3, color="gray", zorder=1)
    ax.set_extent(extent, crs=ccrs.PlateCarree())


def _draw_branches(ax, n, table, color, width_scale_mw):
    for _, row in table.iterrows():
        x = [n.buses.at[row.bus0, "x"], n.buses.at[row.bus1, "x"]]
        y = [n.buses.at[row.bus0, "y"], n.buses.at[row.bus1, "y"]]
        width = 0.3 + row.capacity_mw / width_scale_mw
        ax.plot(x, y, transform=ccrs.PlateCarree(), color=color,
                linewidth=width, zorder=2, solid_capstyle="round")


def _annotate_branches(ax, n, table, extent, name_len=5, max_names=3, edge_margin=0.03):
    '''Label parallel branches once, at their midpoint, with summed capacity
    (GW) above and up to `max_names` branch names below. The extent is
    static (see `EXTENT`), so a branch's true midpoint can fall outside the
    visible frame (e.g. a DE-neighbour interconnector under the tight
    `"de_and_interconnectors"`/`"de_internal"` extents) - rather than
    silently dropping that label, it's clipped to just inside the extent's
    edge, so every drawn branch still gets a visible annotation.'''
    lon_min, lon_max, lat_min, lat_max = extent
    lon_pad = (lon_max - lon_min) * edge_margin
    lat_pad = (lat_max - lat_min) * edge_margin

    key = table.apply(lambda r: tuple(sorted((r.bus0, r.bus1))), axis=1)
    for (bus0, bus1), grp in table.groupby(key):
        mx = (n.buses.at[bus0, "x"] + n.buses.at[bus1, "x"]) / 2
        my = (n.buses.at[bus0, "y"] + n.buses.at[bus1, "y"]) / 2
        mx = np.clip(mx, lon_min + lon_pad, lon_max - lon_pad)
        my = np.clip(my, lat_min + lat_pad, lat_max - lat_pad)

        total_gw = grp.capacity_mw.sum() / 1e3
        shown = [str(idx)[:name_len] for idx in grp.index[:max_names]]
        if len(grp.index) > max_names:
            shown.append(f"+{len(grp.index) - max_names}")

        ax.text(mx, my, f"{total_gw:.1f}\n({','.join(shown)})", transform=ccrs.PlateCarree(),
                fontsize=8, ha="center", va="center", color="black", zorder=4, linespacing=1.1,
                path_effects=[pe.withStroke(linewidth=2, foreground="white")])


def _capacity_legend(ax, width_scale_mw, sizes_gw=(5, 15, 30)):
    handles = [plt.Line2D([0], [0], color="gray", lw=0.3 + s * 1e3 / width_scale_mw, label=f"{s} GW")
               for s in sizes_gw]
    ax.legend(handles=handles, loc="lower right", fontsize=7, framealpha=0.9, title="capacity")


def _summarize(table, label):
    total_gw = table.capacity_mw.sum() / 1e3
    total_twkm = (table.capacity_mw * table.length_km).sum() / 1e6
    stats = {"n_branches": len(table), "total_gw": total_gw, "total_twkm": total_twkm}
    print(f"{label}: {stats['n_branches']} branches, {stats['total_gw']:.1f} GW, {stats['total_twkm']:.3f} TWkm")
    return stats


def _stats_box(ax, tables):
    '''Per-component (+ combined, if more than one) branch count/GW/TWkm -
    printed to stdout via `_summarize` and rendered directly on the plot so
    it travels with any exported PNG.'''
    lines = []
    for component, table in tables.items():
        stats = _summarize(table, component)
        lines.append(f"{component}: {stats['n_branches']} branches, {stats['total_gw']:.1f} GW, {stats['total_twkm']:.2f} TWkm")
    if len(tables) > 1:
        stats = _summarize(pd.concat(tables.values()), "combined")
        lines.append(f"combined: {stats['n_branches']} branches, {stats['total_gw']:.1f} GW, {stats['total_twkm']:.2f} TWkm")

    ax.text(0.98, 0.98, "\n".join(lines), transform=ax.transAxes, ha="right", va="top", fontsize=7.5,
            zorder=5, bbox=dict(facecolor="white", alpha=0.85, edgecolor="gray", boxstyle="round,pad=0.4"))


def plot_grid_network(
    n,
    components=("Line", "Link"),
    scope="de_and_interconnectors",
    horizon=None,
    annotate=True,
    figsize=(10, 10),
    export_path=None,
    regions=None,
    home_country="DE",
    width_scale_mw=6000,
    title=None,
):
    '''Plot a PyPSA network's AC lines and/or DC links on a static basemap.

    n:            a `pypsa.Network`.
    components:   "Line", "Link", or a tuple/list of both.
    scope:        "whole_network" | "de_internal" | "de_and_interconnectors"
                  (see `scope_mask`); controls both which branches are kept
                  and the static map extent (`EXTENT`).
    horizon:      planning year for Link build-year gating; see `branch_table`.
    annotate:     label each branch with "<capacity GW> (<name(s)>)".
    figsize:      passed straight to `plt.subplots`.
    export_path:  if given, save the figure there (format inferred from the
                  path's suffix) in addition to displaying it.
    regions:      optional GeoDataFrame of region shapes to draw under the
                  branches (e.g. `gpd.read_file(...regions_onshore....geojson)`).
    home_country: ISO-2 code defining "de_internal"/"de_and_interconnectors".

    Per-component and (if both components are plotted) combined summary
    statistics - branch count, total GW, total TWkm - for the chosen scope
    are printed and rendered directly on the plot. Returns `(fig, ax)`.
    '''
    if isinstance(components, str):
        components = (components,)

    proj = ccrs.EqualEarth()
    fig, ax = plt.subplots(figsize=figsize, subplot_kw={"projection": proj})
    extent = EXTENT[scope]
    _basemap(ax, extent, regions)

    legend_handles = []
    all_buses = []
    tables = {}
    for component in components:
        table = branch_table(n, component, horizon)
        table = table.loc[scope_mask(table, scope, home_country)]
        tables[component] = table

        color = COMPONENT_COLORS[component]
        _draw_branches(ax, n, table, color, width_scale_mw)
        if annotate:
            _annotate_branches(ax, n, table, extent)
        legend_handles.append(plt.Line2D([0], [0], color=color, lw=2.5, label=COMPONENT_LABELS[component]))
        all_buses.append(table[["bus0", "bus1"]].to_numpy().ravel())

    buses = pd.Index(np.concatenate(all_buses)).unique() if all_buses else pd.Index([])
    bx = n.buses.loc[buses]
    ax.scatter(bx.x, bx.y, transform=ccrs.PlateCarree(), s=6, color="firebrick", zorder=3)

    if len(components) > 1:
        component_legend = ax.legend(handles=legend_handles, loc="upper left", fontsize=7,
                                      framealpha=0.9, title="branch type")
        ax.add_artist(component_legend)
    _capacity_legend(ax, width_scale_mw)
    _stats_box(ax, tables)

    ax.set_title(title or f"{getattr(n, 'name', None) or 'network'}\n{'+'.join(components)}, scope={scope}")
    plt.tight_layout()

    if export_path is not None:
        fig.savefig(export_path, dpi=150, bbox_inches="tight")
        print(f"saved {export_path}")
    plt.show()
    return fig, ax


# ---------------------------------------------------------------------------
# Flexibility / scenario-comparison plotting (ported as-is from a sibling
# repo's flexibility_utils.py, plus aggregate_by_keywords from that repo's
# flexibility_analysis.py)
# ---------------------------------------------------------------------------

tech_colors = {
    "AC": "#70af1d",
    "DC": "#8a1caf",
    "solar": "#f9d002",
    "onwind": "#235ebc",
    "hydro": "#298c81",
    "offwind-dc": "#74c6f2",
    "solar-hsat": "#fdb915",
    "offwind-ac": "#6895dd",
    "ror": "#3dbfb0",
    "PHS": "#51dbcc",
    "lignite": "#826837",
    "lignite (+CHP)": "#826837",
    "coal": "black",
    "coal (+CHP)": "black",
    "oil": "#c9c9c9",
    "Oil": "black",
    "uranium": "#ff8c00",
    "none": "",
    "co2": "#f29dae",
    "co2 stored": "#f2385a",
    "co2 sequestered": "#f2682f",
    "gas": "#e0986c",
    "H2": "#bf13a0",
    "battery": "fuchsia",
    "EV battery": "#baf238",
    "urban central heat": "#d15959",
    "urban central water tanks": "#e9977d",
    "urban central solar thermal": "#d7b24c",
    "biogas": "#e3d37d",
    "solid biomass": "#baa741",
    "biomass (+CHP)": "#baa741",
    "methanol": "#FF7B00",
    "home battery": "#80c944",
    "unsustainable biogas": "",
    "unsustainable solid biomass": "#998622",
    "unsustainable bioliquids": "#32CD32",
    "H2 Store": "#bf13a0",
    "gas for industry CC": "#692e0a",
    "urban central air heat pump": "#6cfb6b",
    "battery discharger": "palegreen",
    "urban central resistive heater": "#8cdf85",
    "biomass to liquid CC": "#32DDD2",
    "urban central solid biomass CHP": "#9d9042",
    "battery charger": "#88a75b",
    "home battery discharger": "#3c5221",
    "H2 pipeline retrofitted": "#ba99b5",
    "DAC": "#ff5270",
    "solid biomass for industry CC": "#47411c",
    "land transport oil": "#afafaf",
    "BEV charger": "#baf238",
    "BEV charging": "#baf238",
    "industry methanol": "#468c8b",
    "shipping oil": "#808080",
    "OCGT": "#e0986c",
    "urban central gas boiler": "#b0904f",
    "gas compressing": "#e05b09",
    "kerosene for aviation": "#a1ffe6",
    "HVC to air": "k",
    "SMR": "#870c71",
    "home battery charger": "#5e8032",
    "coal for industry": "#343434",
    "electricity distribution grid": "#97ad8c",
    "urban central water tanks discharger": "#b9816e",
    "urban central solid biomass CHP CC": "#6c5d28",
    "methanolisation": "#00FFBF",
    "biogas to gas": "#e36311",
    "waste CHP CC": "#e3d3ff",
    "Fischer-Tropsch": "#25c49a",
    "process emissions": "#222222",
    "biomass to liquid": "#32CD32",
    "biogas to gas CC": "#e51245",
    "waste CHP": "#e3d37d",
    "H2 pipeline": "#f081dc",
    "Sabatier": "#9850ad",
    "naphtha for industry": "#57ebc4",
    "agriculture machinery oil": "#949494",
    "H2 Fuel Cell": "#c251ae",
    "urban central gas CHP": "#8d5e56",
    "oil refining": "#a6a6a6",
    "solid biomass for industry": "#7a6d26",
    "SMR CC": "#4f1745",
    "H2 Electrolysis": "#ff29d9",
    "gas for industry": "#853403",
    "urban central gas CHP CC": "#6e4e4c",
    "urban central water tanks charger": "#b57a67",
    "process emissions CC": "#000000",
    "urban central heat vent": "#a74747",
    "oil primary": "#d2d2d2",
    "gas primary": "#e05b09",
    "solar rooftop": "#ffea80",
    "land transport EV": "#baf238",
    "industry electricity": "violet",
    "industrial demand": "violet",
    "low-temperature heat for industry": "#8f2727",
    "H2 for industry": "#f073da",
    "agriculture heat": "#d9a5a5",
    "land transport fuel cell": "#6b3161",
    "electricity": "#110d63",
    "demand": "#110d63",
    "agriculture electricity": "#494778",
    "low voltage": "#97ad8c",
    "non-sequestered HVC": "",
    "rural air heat pump": "#5af95d",
    "rural biomass boiler": "#a1a066",
    "rural gas boiler": "#d4722e",
    "rural ground heat pump": "#48f74f",
    "rural heat": "#ff7c7c",
    "rural resistive heater": "#a5ed9d",
    "rural solar thermal": "#f1c069",
    "rural water tanks": "#f7b7a3",
    "rural water tanks charger": "#b3abb0",
    "rural water tanks discharger": "#ba9685",
    "urban decentral air heat pump": "#5af95d",
    "urban decentral biomass boiler": "#b0b87b",
    "urban decentral gas boiler": "#ba8947",
    "urban decentral heat": "#a33c3c",
    "urban decentral resistive heater": "#98e991",
    "urban decentral solar thermal": "#e5bc5a",
    "urban decentral water tanks": "#f2b2a3",
    "urban decentral water tanks charger": "#b3becc",
    "urban decentral water tanks discharger": "#baac9e",
    "urban central lignite CHP": "#f7a889",
    "urban central oil CHP": "#e26952",
    "urban decentral oil boiler": "",
    "urban central coal CHP": "#b40426",
    "rural oil boiler": "",
    "CCGT": "purple",
    "gas (+CHP)": "purple",
    "nuclear": "lime",
    "renewable oil": "darkviolet",
    "renewable gas": "gold",
    "Other": "gray",
    "H2 OCGT": "#3b4cc0",
    "H2 pipeline (Kernnetz)": "#6788ee",
    "H2 retrofit OCGT": "#9abbff",
    "urban central H2 CHP": "#c9d7f0",
    "urban central H2 retrofit CHP": "#edd1c2",
    "interconnectors": "darkorange",
    "coal CHP": "#b40426",
    "H2 CHP": "pink",
    "H2 (+CHP)": "pink",
    "biomass CHP": "#9d9042",
    "Fuel Cell": "#c251ae",
    "import": "orange",
    "export": "purple",
    "interconnectors supply": "orange",
    "interconnectors demand": "purple",
    "PHS charging": "darkgreen",
    "PHS discharging": "darkgreen",
    "heat pump": "red",
    "heat pump (central)": "red",
    "heat pump (decentral)": "firebrick",
    "resistive heater": "khaki",
    "resistive heater (central)": "khaki",
    "resistive heater (decentral)": "darkkhaki",
    "Import/Export": "dimgrey",
    "Sonstige": "gainsboro",
    "industry DSM": "navy",
    "Industry DSM": "navy",
    "industry DSM ramp down": "mediumslateblue",
    "industry DSM compensate": "cornflowerblue",
    "Onshore wind": "#235ebc",
    "Offshore wind": "#6788ee",
    "Solar": "#ffea80",
    "Gas": "#e0986c",
    "Gas CHP": "#e0986c",
    "Pumped storage": "#51dbcc",
    "Battery": "palegreen",
    "Power-to-heat": "crimson",
    "Power-to-heat (central)": "crimson",
    "Power-to-heat (decentral)": "crimson",
    "Electrolysis": "#ff29d9",
    "Iron-air battery": "#c9954d",
    "Others": "grey",
    "V2G": "tomato",
    "Vehicle-to-grid": "tomato",
    "iron-air battery": "#f5e6b3",
    "iron-air battery storage": "#daa520",
    "iron-air battery charger": "#c9954d",
    "iron-air battery discharger": "#8b6f47",
    "electricity distribution grid losses": "#97ad8c",
    "agriculture machinery electricity": "#6b3161",
    "eFuels": "seagreen",
    "Renewable Gas": "gold",
    "Renewable Oil": "darkviolet",
    "Methanol": "#FF7B00",
    "Coal incl. CHP": "#826837",
}

sector_colors = {
    'Electricity': '#110d63',
    'Heat': '#d15959',
    'H2': '#bf13a0',
    'Fuels': '#1abc9c',
    'Gas': '#e0986c',
    'Biomass': '#baa741',
    'Other': 'lightgrey'
    }

# Capacity grouping for scenario-comparison plots (plot_grid_scenario_comparison.py
# capacity_comparison.png) and per-scenario capacity-by-region maps (system_plots.py).
capa_groups = {
    "Wind + Solar": {
        "Onshore wind": ["onwind"],
        "Offshore wind": ["offwind-ac", "offwind-dc"],
        "Solar": ["solar", "solar rooftop", "solar-hsat"],
    },
    "Backup": {
        "Gas": ["OCGT", "CCGT"],
        "Gas CHP": ["urban central gas CHP", "urban central gas CHP CC"],
        "H2": ["H2 turbine", "H2 OCGT", "H2 CCGT", "H2 retrofit OCGT", "H2 retrofit CCGT"],
        "H2 CHP": ["urban central H2 CHP", "urban central H2 retrofit CHP"],
        "Coal incl. CHP": ["coal", "lignite", "urban central coal CHP", "urban central lignite CHP"],
        "Others": [
            "solid biomass",
            "urban central solid biomass CHP",
            "urban central solid biomass CHP CC",
            "waste CHP",
            "waste CHP CC",
            "oil",
            "urban central oil CHP",
            "H2 Fuel Cell",
        ],
    },
    "Storage Discharge": {
        "Pumped storage": ["PHS"],
        "Battery": ["battery discharger", "home battery discharger"],
        "Iron-air battery": ["iron-air battery discharger"],
        "Vehicle-to-grid": ["V2G"],
    },
    "Demand-Side Flex": {
        "Power-to-heat (central)": [
            "urban central air heat pump",
            "urban central resistive heater",
        ],
        "Power-to-heat (decentral)": [
            "rural air heat pump",
            "rural ground heat pump",
            "rural resistive heater",
            "urban decentral air heat pump",
            "urban decentral resistive heater",
        ],
        "Electrolysis": ["H2 Electrolysis"],
        "BEV charging": ["BEV charger"],
        "Industry DSM": ["industry DSM ramp down"],
    },
}

# Energy (not power) capacity - Stores' e_nom_opt plus PHS storage_units'
# implied energy capacity (p_nom_opt * max_hours), grouped the same
# {tech_label: [carriers]} way as capa_groups so it can share the same
# region-map plotting machinery (see system_plots.py's
# plot_group_capacity_maps/plot_group_capacity_diff_maps), despite coming
# from a different statistics call (n.stores/n.storage_units, not
# n.statistics.optimal_capacity - these carriers mostly aren't on
# AC/low-voltage buses). Carrier list matches the `techs` passed to
# system_plots.py's plot_storage_map.
storage_energy_capa_groups = {
    "Storage Energy": {
        "Pumped storage": ["PHS"],
        "Battery": ["battery", "home battery", "EV battery", "iron-air battery"],
        "Heat storage": [
            "rural water tanks",
            "urban central water pits",
            "urban central water tanks",
            "urban decentral water tanks",
        ],
        "H2 Store": ["H2 Store"],
    },
}

# Technology aggregation for plot_electricity_capacity in system_plots.py
# (supply/demand bar chart, distinct from capa_groups' broader categories).
resistive_heater = [
    "urban central resistive heater",
    "rural resistive heater",
    "urban decentral resistive heater",
]
gas_boiler = [
    "urban central gas boiler",
    "rural gas boiler",
    "urban decentral gas boiler",
]
heat_pump = [
    "urban central air heat pump",
    "rural air heat pump",
    "rural ground heat pump",
    "urban decentral air heat pump",
]
solar = ["solar", "solar-hsat"]
offwind = ["offwind-ac", "offwind-dc"]
h2_ocgt = ["H2 OCGT", "H2 retrofit OCGT"]

c1_groups = [resistive_heater, gas_boiler, heat_pump, solar, offwind, h2_ocgt]

c1_groups_name = [
    "resistive heater",
    "gas boiler",
    "heat pump",
    "solar",
    "offwind",
    "H2 OCGT",
]

def df_to_png(df, filename="table.png"):

    if df.empty:
        return

    fig, ax = plt.subplots(figsize=(10, 0.5 + 0.25*len(df)))
    ax.axis('off')
    tbl = ax.table(cellText=df.values,
                   colLabels=df.columns,
                   rowLabels=df.index,
                   loc='center')
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(10)
    fig.savefig(filename, bbox_inches='tight')
    plt.close()


# scenario_colors is keyed by `run` name for a normal cross-run comparison,
# or by `network_id` for the stochastic-grid LT/ST comparison (see
# plot_grid_scenario_comparison.py) - both can appear as plot legend/table
# entries, so both live in the same dict. Grid-topology variants share one
# base hue per "identity" (eev/reduced_10/reduced_40/stochastic); ST
# (portfolio-{X}_on-{Y}_st) entries use X's hue since that's the capacity
# decision being evaluated, shaded by which canvas Y they're dispatched on.
scenario_colors = {
    # Ariadne scenarios (this repo's `run` names)
    'ExPol': '#0d9488',
    'KN2045_Mix': '#065f46',
    'KN2045_Elek': '#059669',
    'KN2045_H2': '#10b981',
    'KN2045_NFniedrig': '#34d399',
    'KN2045_NFhoch': '#6ee7b7',

    # Stochastic-grid LT (topology-*) network variants
    'topology-eev': '#2563eb',
    'topology-reduced_10': '#f59e0b',
    'topology-reduced_40': '#dc2626',
    'topology-stochastic__reduced_10': '#7c3aed',
    'topology-stochastic__reduced_40': '#4c1d95',

    # Stochastic-grid ST (portfolio-{portfolio}_on-{grid_scenario}_st) network variants
    'portfolio-eev_on-reduced_10_st': '#60a5fa',
    'portfolio-eev_on-reduced_40_st': '#1e40af',
    'portfolio-reduced_10_on-reduced_10_st': '#fbbf24',
    'portfolio-reduced_10_on-reduced_40_st': '#b45309',
    'portfolio-reduced_40_on-reduced_10_st': '#f87171',
    'portfolio-reduced_40_on-reduced_40_st': '#991b1b',
    'portfolio-stochastic_on-reduced_10_st': '#a78bfa',
    'portfolio-stochastic_on-reduced_40_st': '#5b21b6',
}

# Short labels for the same network_ids/run names as scenario_colors, for
# x-axis ticks etc. where the full name is too long to fit. "SP" =
# stochastic program (the joint two-stage solve, matching the WS/SP/EEV/...
# terminology in compute_stochastic_metrics.py); "@" reads as "dispatched
# on" for the ST portfolio-{X}_on-{Y}_st networks.
scenario_abbrev = {
    'topology-eev': 'EEV',
    'topology-reduced_10': 'R10',
    'topology-reduced_40': 'R40',
    'topology-stochastic__reduced_10': 'SP-R10',
    'topology-stochastic__reduced_40': 'SP-R40',

    'portfolio-eev_on-reduced_10_st': 'EEV@R10',
    'portfolio-eev_on-reduced_40_st': 'EEV@R40',
    'portfolio-reduced_10_on-reduced_10_st': 'R10@R10',
    'portfolio-reduced_10_on-reduced_40_st': 'R10@R40',
    'portfolio-reduced_40_on-reduced_10_st': 'R40@R10',
    'portfolio-reduced_40_on-reduced_40_st': 'R40@R40',
    'portfolio-stochastic_on-reduced_10_st': 'SP@R10',
    'portfolio-stochastic_on-reduced_40_st': 'SP@R40',
}


# Per-grid_scenario-VALUE abbreviation (i.e. one entry of GRID_SCENARIO_VALUES
# = list(stochastic_grid_scenarios.scenarios.keys()) + ["eev", "stochastic"]
# in stochastic_grid.smk - not a full network_id). scenario_label/
# scenario_color compose full network_id labels/colors out of these parts
# (see _RE_TOPOLOGY* / _RE_PORTFOLIO below) so a *new* grid_scenario entry
# added to config (data/pypsa-de/grid_scenarios has several - e.g.
# "exogenous2030", "optimal2025" - not yet wired into any config's
# stochastic_grid_scenarios.scenarios) gets working abbreviations/colors for
# every network_id it appears in (topology-*, topology-stochastic__*,
# portfolio-*_on-*_st) automatically, without editing this dict - only add
# an entry here if the auto-generated fallback (_auto_value_abbrev) isn't
# good enough.
_grid_value_abbrev = {
    "reduced_10": "R10",
    "reduced_40": "R40",
    "eev": "EEV",
    "stochastic": "SP",
}

_RE_TOPOLOGY_STOCHASTIC = re.compile(r"^topology-stochastic__(?P<gs>.+)$")
_RE_TOPOLOGY = re.compile(r"^topology-(?P<gs>.+)$")
_RE_PORTFOLIO = re.compile(r"^portfolio-(?P<p>.+?)_on-(?P<gs>.+)_st$")


def _auto_value_abbrev(value):
    """Deterministic fallback abbreviation for a grid_scenario/portfolio
    value with no entry in _grid_value_abbrev: leading letters (max 4,
    uppercased) + trailing digits (max 2) - e.g. "exogenous2030" -> "EXOG30",
    "optimal2025" -> "OPTI25" (kept 4 letters rather than "OPT" specifically
    so it can't be confused with the unrelated "optimal"/"optimal_YYYY"
    pathway-network abbreviation, see scenario_label below).
    """
    letters = re.match(r"[a-zA-Z]*", value).group(0) or re.sub(r"[^a-zA-Z]", "", value)
    digits = re.search(r"\d+", value)
    digits = digits.group(0)[-2:] if digits else ""
    return f"{letters[:4].upper()}{digits}"


def _grid_value_label(value):
    return _grid_value_abbrev.get(value, _auto_value_abbrev(value))


def _auto_value_color(seed):
    """Deterministic fallback color for a string with no entry in
    scenario_colors - hashed to a hue so distinct strings (reliably,
    including across repeated runs/machines) get distinct, stable colors
    without hand-picking one."""
    digest = int(hashlib.md5(seed.encode()).hexdigest(), 16)
    hue = (digest % 360) / 360
    r, g, b = colorsys.hls_to_rgb(hue, 0.45, 0.55)
    return f"#{int(r * 255):02x}{int(g * 255):02x}{int(b * 255):02x}"


def scenario_label(network_id):
    """Short display label for a network_id/run name, for table columns and
    axis ticks where the full name is too long to fit.

    Checks scenario_abbrev first (verbatim, hand-picked entries - Ariadne
    `run` names and any grid-scenario network_id worth a specific override),
    then "optimal"/"optimal_YYYY" (the suffix-less pathway-reference
    networks - network_id "optimal" for the stochastic-grid target year and
    "optimal_YYYY" for earlier myopic-pathway years, see
    GRID_EXPORT_NETWORK_IDS_LT in stochastic_grid.smk), then decomposes any
    other topology-*/portfolio-*_st network_id into its grid_scenario/
    portfolio parts and composes their labels (see _grid_value_label) -
    covering every grid_scenario in stochastic_grid_scenarios.scenarios
    automatically, including ones with no manual abbreviation. Only a
    network_id matching none of these falls back to the raw id.
    """
    if network_id in scenario_abbrev:
        return scenario_abbrev[network_id]
    if network_id == "optimal":
        return "OPT"
    if network_id.startswith("optimal_"):
        return f"OPT-{network_id.removeprefix('optimal_')[-2:]}"
    m = _RE_TOPOLOGY_STOCHASTIC.match(network_id)
    if m:
        return f"SP-{_grid_value_label(m.group('gs'))}"
    m = _RE_PORTFOLIO.match(network_id)
    if m:
        return f"{_grid_value_label(m.group('p'))}@{_grid_value_label(m.group('gs'))}"
    m = _RE_TOPOLOGY.match(network_id)
    if m:
        return _grid_value_label(m.group("gs"))
    return network_id


def scenario_color(network_id):
    """Plot color for a network_id/run name, with the same "optimal"/
    "optimal_YYYY" fallback as scenario_label (see there). Anything else not
    in scenario_colors gets a deterministic hash-based color (see
    _auto_value_color) instead of a flat grey, so new grid_scenario entries
    (see scenario_label) stay visually distinguishable in comparison plots
    without a manual color pick."""
    if network_id in scenario_colors:
        return scenario_colors[network_id]
    if network_id == "optimal" or network_id.startswith("optimal_"):
        return "black"
    return _auto_value_color(network_id)


def aggregate_by_keywords(df, groups, exclude_patterns=None):
    """
    Aggregate rows in df according to keyword groups using substring matching.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with row index as technology names.
    groups : dict
        Keys = new aggregated name,
        Values = list of substrings to match in the index (case-insensitive).
        Example: {"battery": ["battery charger", "battery discharger"]}
    exclude_patterns : list of str, optional
        Substrings that protect a row from being matched by any group.
        Rows whose index contains any of these patterns are skipped entirely
        and remain in the output under their original name.
        Example: ["iron-air"] prevents "iron-air battery charger" from being
        absorbed into a "battery charger" group.
        If None or empty, no rows are excluded. Safe to use if patterns are
        not present in the index.

    Returns
    -------
    pd.DataFrame
        DataFrame with grouped rows summed under new_name, ungrouped rows
        retained as-is, and excluded rows preserved under their original names.
    """
    df_out = df.copy()
    exclude_patterns = exclude_patterns or []
    for new_name, keywords in groups.items():
        pattern = "|".join(keywords)
        mask = df_out.index.to_series().str.contains(pattern, case=False)
        for excl in exclude_patterns:
            mask &= ~df_out.index.to_series().str.contains(excl, case=False)
        if mask.any():
            summed = df_out.loc[mask].sum()
            df_out = df_out.drop(df_out.index[mask])
            df_out.loc[new_name] = summed
    return df_out
