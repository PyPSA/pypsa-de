# SPDX-FileCopyrightText: Contributors to PyPSA-DE
# SPDX-License-Identifier: MIT
"""
Audit processed exogenous transmission projects without running Snakemake.

See doc/transmission_project_report.md for usage and accounting conventions.
"""

import argparse
import json
import logging
import re
from colorsys import hsv_to_rgb
from pathlib import Path

import geopandas as gpd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pypsa
from matplotlib.font_manager import FontProperties
from matplotlib.patches import FancyArrowPatch
from matplotlib.textpath import TextToPath
from shapely import wkt

CATEGORIES = ("AC_internal", "AC_interconnectors", "DC_internal", "DC_interconnectors")
PROJECT_FILES = ("new_lines", "new_links", "upgraded_lines", "upgraded_links")
DEFAULT_INTERVALS = ((2025, 2030), (2030, 2035))


def read_csv(path):
    """Read workflow CSVs, including the empty placeholder files it writes."""
    frame = pd.read_csv(path, index_col=0, dtype={"bus0": str, "bus1": str})
    frame.index = frame.index.astype(str)
    if not frame.index.is_unique:
        raise ValueError(f"Duplicate asset/project IDs in {path}")
    return frame


def source_catalog(directories):
    frames = []
    for directory in directories:
        for action in PROJECT_FILES:
            path = Path(directory) / f"{action}.csv"
            if path.exists():
                frame = read_csv(path)
                frames.append(
                    frame.assign(
                        project_id=frame.index,
                        source_file=str(path),
                        project_action=action,
                    )
                )
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def identify_source(asset, row, action, catalog):
    """Recover legacy provenance only when the source match is unique."""
    if pd.notna(row.get("source_file")) and pd.notna(row.get("project_id")):
        return str(row.source_file), str(row.project_id), "recorded"
    if catalog.empty:
        return "", asset, "unresolved"
    candidates = catalog.loc[catalog.project_action.eq(action)]
    if action.startswith("new"):
        candidates = candidates.loc[candidates.project_id.eq(asset)]
    else:
        # Old upgrade outputs replaced the CSV project ID with a network ID.
        # Tags + build year are retained; never infer identity by proximity.
        tags = str(row.get("tags", ""))
        if not tags or tags == "nan":
            return "", asset, "unresolved"
        candidates = candidates.loc[
            candidates.tags.fillna("").eq(tags)
            & candidates.build_year.eq(row.build_year)
        ]
    if len(candidates) == 1:
        match = candidates.iloc[0]
        return match.source_file, match.project_id, "unique legacy match"
    return "", asset, "unresolved" if candidates.empty else "ambiguous"


def active_capacity(row, nominal, year):
    build = float(row.get("build_year", 0))
    lifetime = float(row.get("lifetime", np.inf))
    if pd.isna(build):
        return np.nan
    if pd.isna(lifetime):
        raise ValueError("Missing lifetime in network or project table")
    return float(row[nominal]) if build <= year < build + lifetime else 0.0


def endpoint_order(geometry, buses):
    """Compare source endpoints with physical bus coordinates before region snapping."""
    if not isinstance(geometry, str) or not geometry:
        return "unknown"
    line = wkt.loads(geometry)
    if line.is_empty:
        return "unknown"
    coords = np.vstack([np.asarray(line.coords)[[0, -1]], buses[["x", "y"]].to_numpy()])
    if not np.isfinite(coords).all():
        return "unknown"
    points = gpd.GeoSeries(
        gpd.points_from_xy(coords[:, 0], coords[:, 1]), crs=4326
    ).to_crs(3035)
    xy = np.column_stack([points.x, points.y])
    direct = np.linalg.norm(xy[:2] - xy[2:], axis=1).sum()
    reverse = np.linalg.norm(xy[:2] - xy[2:][::-1], axis=1).sum()
    if abs(direct - reverse) < 1:
        return "unknown"
    return "reversed" if reverse < direct else "same"


def split_base(base, added):
    """
    Move the pieces of split base branches from the project table to the base.

    Pieces replace the base branch named in ``split_of``. They are base assets,
    not additions, and carry its capacity from the start.
    """
    if "split_of" not in added:
        return base, added
    pieces = added[added.split_of.notna()]
    base = pd.concat(
        [base.drop(pieces.split_of.unique()), pieces.reindex(columns=base.columns)]
    )
    return base, added.drop(pieces.index)


def retirement_row(
    asset, row, project, kind, nominal, buses, catalog, start, end, country
):
    """Ledger row of a branch that a project retires without successor (``retired_by``)."""
    year = float(row.build_year) + float(row.lifetime)
    event = start < year <= end
    loss = active_capacity(row, nominal, np.nextafter(year, -np.inf)) if event else 0.0
    hit = catalog.loc[catalog.project_id.eq(project)] if not catalog.empty else catalog
    c0, c1 = buses.loc[[str(row.bus0), str(row.bus1)], "country"].fillna("").astype(str)
    category = ("AC" if kind == "lines" else "DC") + (
        "_internal" if c0 == c1 == country else "_interconnectors"
    )
    return dict(
        asset=asset,
        project_id=project,
        source_file=hit.source_file.iloc[0] if len(hit) == 1 else "",
        project_action=f"retired_{kind}",
        provenance="recorded",
        replaced_asset=asset,
        build_year=year,
        project_status="",
        tags="",
        category=category,
        in_scope=country in (c0, c1),
        bus0=str(row.bus0),
        bus1=str(row.bus1),
        endpoint_order="unknown",
        country0=c0,
        country1=c1,
        capacity_start=active_capacity(row, nominal, start),
        capacity_end=active_capacity(row, nominal, end),
        gross_added=0.0,
        retired=loss,
        net_added=-loss,
        nominal_capacity=float(row[nominal]),
        unit="MVA" if kind == "lines" else "MW",
        length_km=row.get("length", np.nan),
    )


def build_ledger(n, directory, catalog, start, end, country="DE"):
    """Account for additions and retirement adjustments at physical branch level."""
    if end <= start:
        raise ValueError("The end year must be greater than the start year")
    buses = n.buses.copy()
    new_buses = read_csv(Path(directory) / "new_buses.csv")
    if not new_buses.empty:
        buses = pd.concat([buses, new_buses])
    if not buses.index.is_unique:
        raise ValueError("Project buses overlap existing bus IDs")
    rows = []
    background = []
    for kind, base, nominal in (
        ("lines", n.lines, "s_nom"),
        ("links", n.links.query("carrier == 'DC'"), "p_nom"),
    ):
        added = read_csv(Path(directory) / f"new_{kind}.csv")
        adjusted = read_csv(Path(directory) / f"adjust_{kind}.csv")
        base, added = split_base(base, added)
        missing = adjusted.index.difference(base.index)
        if len(missing):
            raise ValueError(
                f"Retirement adjustments reference missing {kind}: {list(missing)}"
            )
        if len(base.index.intersection(added.index)):
            raise ValueError(f"New {kind} overlap base IDs; use the unextended base.nc")
        adjusted_base = base.copy()
        adjusted_base.update(adjusted)
        for asset, row in adjusted_base.iterrows():
            if (
                max(
                    active_capacity(row, nominal, start),
                    active_capacity(row, nominal, end),
                )
                > 0
            ):
                background.append(
                    dict(row, asset=asset, carrier="AC" if kind == "lines" else "DC")
                )
        for asset, row in added.iterrows():
            replacement = row.get("replaced_asset", "")
            if pd.isna(replacement):
                replacement = ""
            if not replacement and asset.endswith("_upgraded"):
                replacement = asset.removesuffix("_upgraded")
            replacement = str(replacement)
            action = ("upgraded_" if replacement else "new_") + kind
            source, project_id, provenance = identify_source(
                asset, row, action, catalog
            )
            before = active_capacity(row, nominal, start)
            after = active_capacity(row, nominal, end)
            # Attribute changes to commissioning events, not the occupancy of
            # a short-lived stage at interval endpoints. A successor retires
            # its predecessor; the predecessor does not also book that loss.
            event = start < float(row.build_year) <= end
            undated = pd.isna(row.build_year)
            addition = float(row[nominal]) if event else np.nan if undated else 0.0
            retirement = 0.0
            if replacement:
                if replacement in added.index:
                    old = added.loc[replacement]
                else:
                    if replacement not in adjusted.index:
                        raise ValueError(
                            f"Upgrade {asset}: missing retirement adjustment for {replacement}"
                        )
                    old = adjusted_base.loc[replacement]
                if event:
                    retirement = active_capacity(
                        old, nominal, np.nextafter(float(row.build_year), -np.inf)
                    )
            # Genuine lifetime expiry without a successor is still a loss.
            expires = float(row.build_year) + float(row.get("lifetime", np.inf))
            successors = added.get("replaced_asset", pd.Series(dtype=object))
            if start < expires <= end and not successors.eq(asset).any():
                retirement += float(row[nominal])
            endpoints = buses.loc[[str(row.bus0), str(row.bus1)]]
            c0, c1 = endpoints.country.fillna("").astype(str)
            if not c0 or not c1:
                raise ValueError(f"Missing endpoint country for {asset}")
            category = ("AC" if kind == "lines" else "DC") + (
                "_internal" if c0 == c1 == country else "_interconnectors"
            )
            rows.append(
                dict(
                    asset=asset,
                    project_id=project_id,
                    source_file=source,
                    project_action=action,
                    provenance=provenance,
                    replaced_asset=replacement,
                    build_year=row.build_year,
                    project_status=row.get("project_status", ""),
                    tags=row.get("tags", ""),
                    category=category,
                    in_scope=(c0 == country or c1 == country),
                    bus0=str(row.bus0),
                    bus1=str(row.bus1),
                    endpoint_order=endpoint_order(row.get("geometry"), endpoints),
                    country0=c0,
                    country1=c1,
                    capacity_start=before,
                    capacity_end=after,
                    gross_added=addition,
                    retired=retirement,
                    net_added=addition - retirement,
                    nominal_capacity=float(row[nominal]),
                    unit="MVA" if kind == "lines" else "MW",
                    length_km=row.get("length", np.nan),
                )
            )
        for asset in adjusted.index[
            adjusted.get(
                "retired_by", pd.Series(index=adjusted.index, dtype=object)
            ).notna()
        ]:
            rows.append(
                retirement_row(
                    asset,
                    adjusted_base.loc[asset],
                    adjusted.at[asset, "retired_by"],
                    kind,
                    nominal,
                    buses,
                    catalog,
                    start,
                    end,
                    country,
                )
            )
        # Standalone retirement adjustments must not silently disappear.
        represented = {
            r["replaced_asset"]
            for r in rows
            if r["category"].startswith("AC" if kind == "lines" else "DC")
        }
        unmatched = adjusted.index.difference(pd.Index(list(represented)))
        if len(unmatched):
            raise ValueError(
                f"Retirement adjustments without replacement projects: {list(unmatched)}"
            )
    columns = """asset project_id source_file project_action provenance replaced_asset
        build_year project_status tags category in_scope bus0 bus1 country0 country1 endpoint_order
        capacity_start capacity_end gross_added retired net_added nominal_capacity unit length_km""".split()
    ledger = pd.DataFrame(rows, columns=columns)
    ledger["commissioning_known"] = ledger.build_year.notna()
    if not ledger.commissioning_known.all():
        logging.getLogger(__name__).warning(
            "Undated projects excluded from interval totals: %s",
            ledger.loc[~ledger.commissioning_known, "project_id"].tolist(),
        )
    ledger["changed"] = (ledger.gross_added.abs() + ledger.retired.abs()) > 1e-6
    return ledger, pd.DataFrame(background), buses


def load_regions(path, level=2):
    regions = gpd.read_file(path).to_crs(3035)
    if "LEVL_CODE" in regions:
        regions = regions.loc[regions.LEVL_CODE.astype(int).eq(level)]
    if "NUTS_ID" not in regions:
        raise ValueError("Region shapes must have a NUTS_ID column")
    if regions.empty or not regions.NUTS_ID.astype(str).str.len().eq(level + 2).all():
        raise ValueError(f"Supply NUTS{level} polygons for this resolution")
    regions = regions.dissolve(by="NUTS_ID")
    return regions


def locate_buses(buses, regions, resolution="nuts2"):
    points = gpd.GeoDataFrame(
        buses[["country"]].copy(),
        geometry=gpd.points_from_xy(buses.x, buses.y),
        crs=4326,
    ).to_crs(3035)
    if resolution == "grid":
        return pd.DataFrame(
            {"region": buses.index, "x": points.geometry.x, "y": points.geometry.y},
            index=buses.index,
        )
    joined = gpd.sjoin(
        points,
        regions[["geometry"]].rename_axis("matched_region"),
        how="left",
        predicate="intersects",
    )
    joined = joined.loc[~joined.index.duplicated()]
    # Offshore/outside-NUTS endpoints retain physical coordinates. Never snap
    # to an arbitrary inland region, which would imply a false model mapping.
    anchors = regions.geometry.representative_point()
    result = pd.DataFrame(index=points.index)
    result["region"] = joined.matched_region.fillna("outside regions")
    result["x"] = [
        anchors.loc[r].x if r in anchors.index else p.x
        for r, p in zip(result.region, points.geometry)
    ]
    result["y"] = [
        anchors.loc[r].y if r in anchors.index else p.y
        for r, p in zip(result.region, points.geometry)
    ]
    return result


def add_region_mapping(ledger, locations, busmap=None, resolution="nuts2"):
    for port in (0, 1):
        ledger[f"region{port}"] = ledger[f"bus{port}"].map(locations.region)
        if resolution != "grid":
            ledger[f"{resolution}_{port}"] = ledger[f"region{port}"]
        if busmap is not None:
            ledger[f"model_bus{port}"] = ledger[f"bus{port}"].map(busmap)
    ledger["within_region"] = ledger.region0.eq(ledger.region1) & ledger.region0.ne(
        "outside regions"
    )
    if resolution != "grid":
        ledger[f"within_{resolution}"] = ledger.within_region
    if busmap is not None:
        ledger["model_mapping"] = np.where(
            ledger.model_bus0.isna() | ledger.model_bus1.isna(),
            "unmapped",
            np.where(
                ledger.model_bus0.eq(ledger.model_bus1),
                "collapsed within model bus",
                "between model buses",
            ),
        )
    return ledger


def project_details(tags):
    """Extract readable references from JSON and PostgreSQL-style source tags."""
    fields = dict(re.findall(r'"([^"\n]+)"\s*(?:=>|:)\s*"([^"\n]*)"', str(tags)))
    keys = [
        k
        for k in fields
        if "proj_id" in k
        or "invest_id" in k
        or k in ("NEP-ID", "Maßnahmen-ID", "NOVA-Typ")
    ]
    return " | ".join(f"{k}: {fields[k]}" for k in keys if fields[k])


def project_metadata(row, catalog):
    """Use supplied source metadata; endpoint descriptions are not official names."""
    tags = str(row.tags)
    fields = dict(re.findall(r'"([^"\n]+)"\s*(?:=>|:)\s*"([^"\n]*)"', tags))
    for key, value in re.findall(r"(?:^|[,{])\s*([\w -]+):\s*([^,}\]]+)", tags):
        fields.setdefault(key.strip(), value.strip())
    dataset = Path(row.source_file).parent.name if row.source_file else "unresolved"
    tyndp = re.search(r"tyndp[ _-]?(\d{4})", f"{row.project_id} {dataset} {tags}", re.I)
    plan = fields.get("origin") or (f"TYNDP {tyndp[1]}" if tyndp else dataset.upper())
    if dataset == "manual":
        plan += " (manual CSV)"
    reference = fields.get("NEP-ID") or next(
        (v for k, v in fields.items() if "proj_id" in k), str(row.project_id)
    )
    name = fields.get("name") or fields.get("project_name") or fields.get("title")
    endpoints = [
        fields.get("location0") or fields.get("substation_from", ""),
        fields.get("location1") or fields.get("substation_to", ""),
    ]
    if not catalog.empty:
        source = catalog.loc[
            catalog.project_id.eq(str(row.project_id))
            & catalog.source_file.eq(row.source_file)
        ]
        if len(source) == 1:
            endpoints = [
                source.iloc[0].get(f"bus{port}", "")
                if pd.notna(source.iloc[0].get(f"bus{port}"))
                and source.iloc[0].get(f"bus{port}")
                else endpoints[port]
                for port in (0, 1)
            ]
    endpoints = [str(value).strip() if pd.notna(value) else "" for value in endpoints]
    source_endpoints = endpoints.copy()
    order = getattr(row, "endpoint_order", "same")
    if order == "reversed":
        endpoints.reverse()
    if not name:
        named = [value for value in endpoints if value]
        name = " → ".join(named) + " (endpoints)" if named else "Name not supplied"
        if order == "unknown" and named:
            name = (
                " ↔ ".join(source_endpoints) + " (source endpoints; order unresolved)"
            )
    if order == "unknown":
        endpoints = ["", ""]
    url = re.search(
        r'https?://[^\s"\'<>,}\]]+',
        fields.get("source_url") or fields.get("url") or tags,
    )
    page = re.search(r'"source_pdf_page"\s*:\s*"?(\d+)', tags)
    return dict(
        source_plan=plan,
        reference_id=reference,
        project_name=name,
        source_url=url[0] if url else "",
        source_page=int(page[1]) if page else None,
        source_endpoint0_name=source_endpoints[0],
        source_endpoint1_name=source_endpoints[1],
        endpoint0_name=endpoints[0] or getattr(row, "bus0", "unknown"),
        endpoint1_name=endpoints[1] or getattr(row, "bus1", "unknown"),
    )


def apply_source_links(metadata, ledger, path):
    """
    Override source URL and PDF page from a lookup table.

    Points datasets whose tags lack a working link or page (e.g. the default
    ``nep`` folder) to the exact source page without editing them.

    Parameters
    ----------
    metadata : pandas.DataFrame
        project_metadata per ledger row, indexed like ledger.
    ledger : pandas.DataFrame
        Project ledger with project_id and source_file.
    path : str or pathlib.Path
        CSV with columns dataset (source folder name), project_id,
        source_url and source_pdf_page.

    Returns
    -------
    pandas.DataFrame
        metadata with source_url and source_page overridden where listed.
    """
    links = pd.read_csv(path, dtype={"project_id": str})
    links = links.set_index(["dataset", "project_id"])
    keys = pd.MultiIndex.from_arrays(
        [
            ledger.source_file.fillna("").map(lambda f: Path(f).parent.name),
            ledger.project_id.astype(str),
        ]
    )
    found = keys.isin(links.index)
    metadata = metadata.copy()
    metadata.loc[found, "source_url"] = links.loc[keys[found], "source_url"].values
    metadata.loc[found, "source_page"] = links.loc[
        keys[found], "source_pdf_page"
    ].values
    return metadata


def visible_projects(ledger):
    """Exclude branches collapsed at display resolution or by the supplied busmap."""
    visible = ledger.in_scope & ledger.changed & ~ledger.within_region
    if "model_mapping" in ledger:
        visible &= ledger.model_mapping.ne("collapsed within model bus")
    return visible


def volume(selected, changed):
    """
    Count measures that do not depend on the region size.

    The displayed ``net_added`` sums every branch that crosses a region border, so an
    upgrade chain crossing k borders counts k times and branches within a region not at all.

    Parameters
    ----------
    selected : pandas.DataFrame
        Shown ledger rows of one category.
    changed : pandas.DataFrame
        All changed in-scope ledger rows of the same category, shown or not.

    Returns
    -------
    dict
        ``net_added_once``: shown capacity with each project counted once, as the
        length-weighted mean over its shown branches (MVA or MW).
        ``shown_branches_per_project``: shown branches per shown project.
        ``net_twkm`` and ``shown_net_twkm``: sum of net_added x length_km / 1e6 over all
        changed or only shown branches (TVA km for lines, TW km for links).
    """

    def twkm(rows):
        return float((rows.net_added * rows.length_km).sum()) / 1e6

    projects = selected.groupby("project_id")
    weighted = (
        (selected.net_added * selected.length_km).groupby(selected.project_id).sum()
    )
    once = (weighted / projects.length_km.sum()).fillna(projects.net_added.mean())
    return dict(
        net_added_once=float(once.sum()),
        shown_branches_per_project=len(selected) / len(once) if len(once) else 0.0,
        net_twkm=twkm(changed),
        shown_net_twkm=twkm(selected),
    )


def wrap_legend(text, width):
    """Wrap to the available page width in points, including long source URLs."""
    font = FontProperties(family="DejaVu Sans", size=7)
    measure = TextToPath()
    lines = []
    for paragraph in text.splitlines():
        while paragraph:
            low, high = 1, len(paragraph)
            while low <= high:
                middle = (low + high) // 2
                if (
                    measure.get_text_width_height_descent(
                        paragraph[:middle], font, False
                    )[0]
                    <= width
                ):
                    low = middle + 1
                else:
                    high = middle - 1
            end = max(1, high)
            if end < len(paragraph):
                space = paragraph.rfind(" ", 0, end + 1)
                end = space if space > 0 else end
            lines.append(paragraph[:end].rstrip())
            paragraph = paragraph[end:].lstrip()
    return "\n".join(lines)


def draw_report(
    ledger,
    background,
    locations,
    regions,
    category,
    start,
    end,
    country,
    output,
    resolution="nuts2",
    note="Before clustering; no endogenous expansion.",
):
    """Draw the map once, then paginate its project legend on A4 sheets."""
    from matplotlib.backends.backend_pdf import PdfPages

    selected = ledger.loc[
        visible_projects(ledger) & ledger.category.eq(category)
    ].copy()
    selected = selected.sort_values(["build_year", "project_id", "asset"])
    selected["key"] = [f"{i + 1:02d}" for i in range(len(selected))]
    selected["line_color"] = [
        hsv_to_rgb((i * 0.61803398875) % 1, 0.72, 0.65) for i in range(len(selected))
    ]
    color = "#126b83" if category.startswith("AC") else "#b1456b"
    unit = "MVA" if category.startswith("AC") else "MW"
    home = regions.loc[regions.index.str.startswith(country)]
    if home.empty:
        raise ValueError(f"No region shapes for country {country}")
    minx, miny, maxx, maxy = home.total_bounds
    if not selected.empty:
        points = locations.loc[pd.unique(selected[["bus0", "bus1"]].values.ravel())]
        minx, miny = min(minx, points.x.min()), min(miny, points.y.min())
        maxx, maxy = max(maxx, points.x.max()), max(maxy, points.y.max())
    padx, pady = (maxx - minx) * 0.12, (maxy - miny) * 0.08

    chunks, chunk, used = [], [], 0
    for row in selected.itertuples():
        lines = [
            f"[{row.key}] {row.project_id} | {int(row.build_year)} | {row.net_added:+,.0f} {unit} | "
            f"{row.source_plan} | Folder: {Path(row.source_file).parent.name if row.source_file else 'unknown'} | "
            f"Project ID: {row.reference_id} | {row.project_name}",
            f"{row.endpoint0_name} ({row.region0}) > {row.endpoint1_name} ({row.region1}) | {row.project_status} | "
            f"{row.project_action} | gross {row.gross_added:+,.0f} | retired {row.retired:,.0f}"
            + (f" | replaces {row.replaced_asset}" if row.replaced_asset else ""),
            project_details(row.tags),
        ]
        entry = wrap_legend(
            " | ".join(line for line in lines if line), width=0.76 * 210 / 25.4 * 72 - 2
        )
        link, label = row.source_url, "source \u2197" if row.source_url else ""
        if link and pd.notna(row.source_page):
            link += f"#page={int(row.source_page)}"
            label = f"p. {int(row.source_page)} \u2197"
        size = entry.count("\n") + 2
        if chunk and used + size > (24 if not chunks else 72):
            chunks.append(chunk)
            chunk, used = [], 0
        chunk.append((entry, row.line_color, link, label))
        used += size
    omitted = ledger.loc[
        ledger.in_scope
        & ledger.changed
        & ~visible_projects(ledger)
        & ledger.category.eq(category)
    ]
    counts = omitted.project_id.value_counts().sort_index()
    if not counts.empty:
        omitted_text = " | ".join(
            f"{project_id} ({count})" for project_id, count in counts.items()
        )
        entry = wrap_legend(
            "Omitted within-region/model-bus projects: " + omitted_text,
            width=0.86 * 210 / 25.4 * 72 - 2,
        )
        if chunk and used + entry.count("\n") + 2 > (24 if not chunks else 72):
            chunks.append(chunk)
            chunk = []
        chunk.append((entry, "#56656e", "", ""))
    chunks.append(chunk)

    style = {
        "font.family": "DejaVu Sans",
        "font.size": 8,
        "pdf.fonttype": 42,
        "figure.autolayout": False,
    }
    with plt.rc_context(style), PdfPages(output) as pdf:
        fig = plt.figure(figsize=(210 / 25.4, 297 / 25.4), facecolor="#fbfaf7")
        fig.text(
            0.07,
            0.955,
            f"{country} / TRANSMISSION PROJECT AUDIT",
            color=color,
            weight="bold",
        )
        fig.text(0.07, 0.923, category.replace("_", " "), fontsize=23, weight="bold")
        fig.text(
            0.07,
            0.897,
            f"{start} until {end} | nominal ratings | {resolution.upper()} endpoints",
            fontsize=10,
        )
        totals = (
            f"Net {selected.net_added.sum():+,.0f} {unit}     "
            f"Gross {selected.gross_added.sum():+,.0f}     Retired {selected.retired.sum():,.0f}"
        )
        fig.text(0.07, 0.868, totals, fontsize=12, color=color, weight="bold")
        ax = fig.add_axes((0.06, 0.37, 0.88, 0.47), facecolor="#eaf0f2")
        regions.cx[minx - padx : maxx + padx, miny - pady : maxy + pady].simplify(
            500
        ).plot(ax=ax, facecolor="#f1f1ed", edgecolor="#c9cfd0", linewidth=0.35)
        home.plot(ax=ax, facecolor="#e0e8e8", edgecolor="white", linewidth=0.6)
        branches = (
            background.loc[background.carrier.eq(category[:2])]
            if not background.empty
            else background
        )
        for branch in branches.itertuples():
            a, b = locations.loc[branch.bus0], locations.loc[branch.bus1]
            ax.plot(
                [a.x, b.x], [a.y, b.y], color="#afbdc3", lw=0.45, alpha=0.5, zorder=1
            )
        endpoints = regions.index.isin(selected.region0) | regions.index.isin(
            selected.region1
        )
        labels = regions.loc[
            endpoints
            | (
                regions.index.str.startswith(country)
                if resolution != "nuts3"
                else False
            )
        ]
        for code, region in labels.iterrows():
            point = region.geometry.representative_point()
            name = region.get("NAME_LATN", region.get("NUTS_NAME", ""))
            label = f"{code}\n{name}" if pd.notna(name) and name else code
            ax.text(
                point.x,
                point.y,
                label,
                fontsize=4.5,
                color="#73818a",
                ha="center",
                va="center",
                zorder=2,
            )
        pairs = {}
        for row in selected.itertuples():
            a, b = locations.loc[row.bus0], locations.loc[row.bus1]
            pair = tuple(sorted((row.region0, row.region1)))
            count = pairs.get(pair, 0)
            pairs[pair] = count + 1
            rad = (-1 if count % 2 else 1) * (0.08 + 0.12 * (count // 2))
            if row.within_region or (a.x == b.x and a.y == b.y):
                radius = (maxx - minx) * (0.022 + count * 0.012)
                theta = np.linspace(0.1, 2 * np.pi - 0.1, 80)
                ax.plot(
                    a.x + radius * np.sin(theta),
                    a.y + radius * (1 - np.cos(theta)),
                    color=row.line_color,
                    lw=1.8,
                )
                label = (a.x, a.y + 2 * radius)
            else:
                ax.add_patch(
                    FancyArrowPatch(
                        (a.x, a.y),
                        (b.x, b.y),
                        arrowstyle="-",
                        connectionstyle=f"arc3,rad={rad}",
                        linewidth=min(5, 1.2 + abs(row.net_added) / 2500),
                        color=row.line_color,
                        zorder=3,
                    )
                )
                label = (
                    (a.x + b.x) / 2 + rad * (b.y - a.y) / 2,
                    (a.y + b.y) / 2 - rad * (b.x - a.x) / 2,
                )
            ax.scatter(
                [a.x, b.x],
                [a.y, b.y],
                s=12,
                color=row.line_color,
                edgecolor="white",
                lw=0.4,
                zorder=4,
            )
            ax.annotate(
                row.key,
                label,
                fontsize=6,
                ha="center",
                zorder=5,
                color=row.line_color,
                weight="bold",
                bbox=dict(
                    boxstyle="round,pad=.2", fc="white", ec=row.line_color, lw=0.4
                ),
            )
        ax.set(
            xlim=(minx - padx, maxx + padx),
            ylim=(miny - pady, maxy + pady),
            aspect="equal",
        )
        ax.set_axis_off()
        if selected.empty:
            ax.text(
                0.5,
                0.5,
                "No project capacity changes in this interval",
                transform=ax.transAxes,
                ha="center",
            )
        legend_title = fig.text(
            0.07,
            0.345,
            f"PROJECT KEY / numbered connections; net change in {unit}",
            color=color,
            weight="bold",
        )
        legend = []
        footer = fig.text(0.07, 0.035, "", fontsize=7, color="#56656e")
        for page, entries in enumerate(chunks, 1):
            ax.set_visible(page == 1)
            legend_title.set_y(0.345 if page == 1 else 0.84)
            for text in legend:
                text.remove()
            legend = []
            y = 0.325 if page == 1 else 0.82
            for entry, entry_color, link, label in entries or [
                ("No additions or retirements.", "#223441", "", "")
            ]:
                legend.append(
                    fig.text(
                        0.07,
                        y,
                        entry,
                        va="top",
                        fontsize=7,
                        linespacing=1.25,
                        color=entry_color,
                    )
                )
                if link:
                    legend.append(
                        fig.text(
                            0.93,
                            y,
                            label,
                            va="top",
                            ha="right",
                            fontsize=7,
                            color=entry_color,
                            url=link,
                        )
                    )
                y -= (entry.count("\n") + 2) * 7 * 1.25 / (297 / 25.4 * 72)
            footer.set_text(
                "Within-region / within-model-bus projects omitted; retained in project_ledger.csv.\n"
                f"{note} Page {page}/{len(chunks)}"
            )
            pdf.savefig(fig)
        plt.close(fig)
    return selected[["asset", "key"]]


def generate_report(
    base_network,
    projects,
    regions_path,
    start,
    end,
    output,
    country="DE",
    source_dirs=(),
    busmap=None,
    resolution="nuts2",
    source_links=None,
):
    if resolution not in ("nuts2", "nuts3", "grid"):
        raise ValueError("Resolution must be nuts2, nuts3, or grid")
    logging.getLogger("fontTools.subset").setLevel(logging.WARNING)
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    catalog = source_catalog(source_dirs)
    n = pypsa.Network(base_network)
    ledger, background, buses = build_ledger(n, projects, catalog, start, end, country)
    regions = load_regions(regions_path, level=3 if resolution == "nuts3" else 2)
    locations = locate_buses(buses, regions, resolution)
    mapping = None
    if busmap:
        mapping = pd.read_csv(busmap, index_col=0, dtype=str).squeeze("columns")
        if not isinstance(mapping, pd.Series):
            raise ValueError(
                "--busmap must contain one value column mapping physical buses to final model buses"
            )
        mapping.index = mapping.index.astype(str)
    ledger = add_region_mapping(ledger, locations, mapping, resolution)
    metadata = pd.DataFrame(
        [project_metadata(row, catalog) for row in ledger.itertuples()],
        index=ledger.index,
        columns=[
            "source_plan",
            "reference_id",
            "project_name",
            "source_url",
            "source_page",
            "source_endpoint0_name",
            "source_endpoint1_name",
            "endpoint0_name",
            "endpoint1_name",
        ],
    )
    if source_links:
        metadata = apply_source_links(metadata, ledger, source_links)
    ledger = pd.concat([ledger, metadata], axis=1)
    ledger["shown"] = visible_projects(ledger)
    summary = {}
    for category in CATEGORIES:
        keys = draw_report(
            ledger,
            background,
            locations,
            regions,
            category,
            start,
            end,
            country,
            output / f"{category}_{start}_until_{end}.pdf",
            resolution,
        )
        ledger.loc[
            ledger.asset.isin(keys.asset) & ledger.category.eq(category), "map_key"
        ] = ledger.asset.map(keys.set_index("asset").key)
        selected = ledger.loc[ledger.shown & ledger.category.eq(category)]
        excluded = ledger.loc[
            ledger.in_scope
            & ledger.changed
            & ~ledger.shown
            & ledger.category.eq(category)
        ]
        summary[category] = {
            "unit": "MVA" if category.startswith("AC") else "MW",
            **{
                col: float(selected[col].sum())
                for col in ("gross_added", "retired", "net_added")
            },
            "changed_projects": int(selected.changed.sum()),
            "undated_projects": int(
                (
                    ledger.in_scope
                    & ledger.category.eq(category)
                    & ~ledger.commissioning_known
                ).sum()
            ),
            "excluded_projects": len(excluded),
            "excluded_net_added": float(excluded.net_added.sum()),
            **volume(selected, pd.concat([selected, excluded])),
        }
        if mapping is not None:
            summary[category]["net_by_model_mapping"] = {
                str(key): float(value)
                for key, value in selected.groupby("model_mapping")
                .net_added.sum()
                .items()
            }
    ledger.to_csv(output / "project_ledger.csv", index=False)
    ledger["corridor0"] = ledger[["region0", "region1"]].min(axis=1)
    ledger["corridor1"] = ledger[["region0", "region1"]].max(axis=1)
    corridors = (
        ledger.loc[ledger.shown]
        .groupby(["category", "corridor0", "corridor1", "unit"], dropna=False)
        .agg(
            net_added=("net_added", "sum"),
            gross_added=("gross_added", "sum"),
            retired=("retired", "sum"),
            projects=("project_id", lambda x: " | ".join(x)),
        )
        .reset_index()
    )
    corridors.to_csv(output / f"{resolution}_corridors.csv", index=False)
    if not catalog.empty:
        accepted = set(zip(ledger.source_file, ledger.project_id))
        catalog["processing_result"] = [
            "included"
            if (r.source_file, r.project_id) in accepted
            else "not found in processed output; inspect build log/config"
            for _, r in catalog.iterrows()
        ]
    if catalog.empty:
        catalog = pd.DataFrame(
            columns=["source_file", "project_id", "processing_result"]
        )
    catalog.to_csv(output / "source_audit.csv", index=False)
    manifest = dict(
        start_year=start,
        end_year=end,
        country=country,
        base_network=str(base_network),
        processed_projects=str(projects),
        regions=str(regions_path),
        resolution=resolution,
        source_dirs=list(map(str, source_dirs)),
        busmap=str(busmap) if busmap else None,
        source_links=str(source_links) if source_links else None,
        accounting="Commissioning-event changes within start < year <= end; each successor retires its immediate predecessor once. Endpoint capacities retain build_year <= year < build_year + lifetime. Displayed totals exclude within-region and known within-model-bus branches; full physical accounting remains in project_ledger.csv.",
        categories=summary,
    )
    (output / "summary.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return ledger, summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--base-network", required=True, help="Unextended physical networks/base.nc"
    )
    parser.add_argument(
        "--projects",
        required=True,
        help="Processed transmission_projects directory (five workflow CSVs)",
    )
    parser.add_argument(
        "--regions",
        "--nuts2",
        dest="regions",
        help="NUTS polygons; grid mode uses NUTS2 as backdrop",
    )
    parser.add_argument(
        "--resolution", choices=("nuts2", "nuts3", "grid"), default="nuts2"
    )
    parser.add_argument(
        "--start-year", type=int, help="Custom interval; requires --end-year"
    )
    parser.add_argument(
        "--end-year", type=int, help="Custom interval; requires --start-year"
    )
    parser.add_argument("--country", default="DE")
    parser.add_argument(
        "--source-dir",
        action="append",
        default=[],
        help="Raw source directory; repeat for NEP, manual, TYNDP",
    )
    parser.add_argument(
        "--busmap",
        help="Optional composed physical-to-model busmap CSV, one value column",
    )
    parser.add_argument(
        "--source-links",
        help="Optional CSV (dataset, project_id, source_url, source_pdf_page) "
        "overriding source links without editing the project data",
    )
    parser.add_argument(
        "--output-dir",
        help="Output root for default intervals; custom interval writes here directly",
    )
    args = parser.parse_args()
    if (args.start_year is None) != (args.end_year is None):
        parser.error(
            "Supply both --start-year and --end-year, or omit both for the default intervals"
        )
    custom = args.start_year is not None
    intervals = ((args.start_year, args.end_year),) if custom else DEFAULT_INTERVALS
    root = Path(
        args.output_dir or f"results/transmission_project_audit/{args.resolution}"
    )
    for start, end in intervals:
        output = root if custom and args.output_dir else root / f"{start}_until_{end}"
        generate_report(
            args.base_network,
            args.projects,
            args.regions
            or f"data/nuts/NUTS_RG_01M_2021_4326_LEVL_{3 if args.resolution == 'nuts3' else 2}.geojson",
            start,
            end,
            output,
            args.country,
            args.source_dir,
            args.busmap,
            args.resolution,
            args.source_links,
        )


if __name__ == "__main__":
    if "snakemake" in globals():
        snakemake = globals()["snakemake"]
        generate_report(
            snakemake.input.base_network,
            Path(snakemake.input.new_lines).parent,
            snakemake.input.regions,
            int(snakemake.wildcards.report_start),
            int(snakemake.wildcards.report_end),
            Path(snakemake.output.summary).parent,
            source_dirs=snakemake.input.source_dirs,
            resolution=snakemake.wildcards.report_resolution,
        )
    else:
        main()
