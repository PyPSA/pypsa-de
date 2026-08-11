# SPDX-FileCopyrightText: Contributors to PyPSA-DE <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT
"""
General-purpose helpers for ad-hoc PyPSA-DE analysis notebooks/scripts - not
part of the main pipeline DAG, and not tied to any one notebook. The idea is
a shared landing spot for reusable notebook helpers as they get factored out,
instead of every notebook re-implementing its own copy.

Currently home to grid-topology plotting (`plot_grid_network`, factored out
of notebooks/grid_topology.ipynb); more sections can be added alongside it as
other notebooks contribute their own reusable pieces.

Imported interactively, e.g.:

    import sys
    sys.path.append(str(REPO / "scripts" / "pypsa-de"))
    from utils import plot_grid_network

    plot_grid_network(n, components="Line", scope="de_and_interconnectors")
"""

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import numpy as np
import pandas as pd

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
