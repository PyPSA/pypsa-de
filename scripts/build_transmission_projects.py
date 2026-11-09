# SPDX-FileCopyrightText: Contributors to PyPSA-Eur <https://github.com/pypsa/pypsa-eur>
#
# SPDX-License-Identifier: MIT


"""
Gets the transmission projects defined in the config file, concatenates and
deduplicates them. Projects are later included in [add_electricity][].

Inputs
------

- `networks/base_network.nc`:  Base network topology for the electricity grid. This is processed in [base_network][].
- `data/transmission_projects/"project_name"/`: Takes the transmission projects from the subfolder of data/transmission_projects. The subfolder name is the project name.
- `offshore_shapes.geojson`: Shapefile containing the offshore regions. Used to determine if a new bus should be added for a new line or link.
- `europe_shape.geojson`: Shapefile containing the shape of Europe. Used to determine if a project is within the considered countries.
- `data/transmission_projects/line_types.csv`: Reviewed line types missing from the PyPSA catalogue (e.g. 3.6 kA 380 kV, 220-kV HTLS).

Outputs
-------

- `transmission_projects/new_lines.csv`: New project lines to be added to the network. This includes new lines and upgraded lines, and the pieces of base branches split at reviewed junctions (`split_of` names the branch they replace).
- `transmission_projects/new_links.csv`: New project links to be added to the network. This includes new links and upgraded links.
- `transmission_projects/adjust_lines.csv`: For lines which are upgraded, the decommissioning year of the existing line is adjusted to the build year of the upgraded line. Lines a project retires without successor carry the project in `retired_by`.
- `transmission_projects/adjust_links.csv`: For links which are upgraded, the decommissioning year of the existing link is adjusted to the build year of the upgraded link.
- `transmission_projects/new_buses.csv`: For some links, we have to add new buses (e.g. North Sea Wind Power Hub).
- `transmission_projects/line_types.csv`: The reviewed line types, registered by [add_transmission_projects_and_dlr][] before the lines.
"""

import json
import logging
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import pypsa
import shapely
from pypsa.descriptors import nominal_attrs
from scipy import sparse, spatial
from scipy.sparse.csgraph import dijkstra
from shapely import ops, wkt
from shapely.geometry import LineString, Point

from scripts._helpers import PYPSA_V1, configure_logging, set_scenario_config
from scripts._transmission_overlaps import reconcile

logger = logging.getLogger(__name__)


def select_project_path(project, paths):
    matching = [p for p in paths if Path(p).name == project]
    if len(matching) != 1:
        raise ValueError(f"Expected exactly one folder named {project}, got {matching}")
    return matching[0]


if PYPSA_V1:
    pypsa.options.params.add.return_names = True


def add_new_buses(n, new_ports, suffix=" bus"):
    # Add new buses for the ports which do not have an existing bus close by. If there are multiple ports at the same location, only one bus is added.
    duplicated = new_ports.duplicated(subset=["x", "y"], keep="first")
    to_add = new_ports[~duplicated]
    names = to_add.index.astype(str) + suffix
    if names.isin(n.buses.index).any():
        raise ValueError(
            "Offshore bus name already exists; check project IDs and endpoint coordinates"
        )
    added_buses = n.add(
        "Bus",
        name=to_add.index,
        suffix=suffix,
        x=to_add.x,
        y=to_add.y,
        v_nom=380,
        under_construction=True,
        symbol="substation",
        substation_off=True,
        substation_lv=False,
        carrier="AC",
    )
    new_buses = n.buses.loc[added_buses].copy().dropna(axis=1, how="all")
    new_ports.loc[to_add.index, "neighbor"] = added_buses
    new_ports["neighbor"] = new_ports.groupby(["x", "y"])["neighbor"].transform("first")
    return new_ports, new_buses


def find_country_for_bus(bus, shapes):
    """
    Find the country of a bus based on its coordinates and the provided
    shapefile.

    Shapefile must contain a column "country" with the country names.
    """
    point = Point(bus.x, bus.y)
    country = shapes.loc[shapes.contains(point), "country"]
    return country.values[0]


def connect_new_lines(
    lines,
    n,
    new_buses_df,
    offshore_shapes=None,
    distance_upper_bound=np.inf,
    bus_carrier="AC",
):
    """
    Find the closest existing bus to the port of each line.

    If closest bus is further away than distance_upper_bound and is
    inside an offshore region, a new bus is created. and the line is
    connected to it.
    """
    bus_carrier = np.atleast_1d(bus_carrier)
    for port in [0, 1]:
        # Explicit reviewed onshore transition points (e.g. cable/overhead
        # junctions) must not snap across a regional boundary to a nearby bus.
        junction_column = f"junction{port}"
        junctions = {}
        if junction_column in lines:
            for ident, row in lines.iterrows():
                junction = row[junction_column]
                if pd.isna(junction) or not str(junction).strip():
                    continue
                name = f"project_junction:{junction}"
                x, y = get_bus_coords_from_port(row["geometry"], port)
                if name in n.buses.index:
                    if not np.allclose(
                        n.buses.loc[name, ["x", "y"]].astype(float),
                        [x, y],
                        atol=1e-9,
                        rtol=0,
                    ):
                        raise ValueError(
                            f"Inconsistent coordinates for junction {junction}"
                        )
                else:
                    country = row.get(f"junction{port}_country")
                    if pd.isna(country) or not str(country).strip():
                        raise ValueError(
                            f"Missing reviewed country for junction {junction}"
                        )
                    n.add(
                        "Bus",
                        name,
                        x=x,
                        y=y,
                        v_nom=float(row.v_nom),
                        carrier="AC",
                        country=country,
                        under_construction=True,
                        symbol="substation",
                        substation_lv=False,
                        substation_off=False,
                    )
                    new_buses_df = pd.concat([new_buses_df, n.buses.loc[[name]].copy()])
                junctions[ident] = name
        # Include buses created for the previous endpoint when matching this one.
        buses = n.buses.query("carrier in @bus_carrier").copy()
        bus_tree = spatial.KDTree(buses[["x", "y"]])
        lines_port = lines["geometry"].apply(
            lambda x: pd.Series(
                get_bus_coords_from_port(x, port=port), index=["x", "y"]
            )
        )
        distances, indices = bus_tree.query(lines_port)
        # Series of lines with closest bus in the existing network and whether they match the distance criterion
        lines_port["neighbor"] = buses.iloc[indices].index
        lines_port["match_distance"] = distances < distance_upper_bound
        for ident, name in junctions.items():
            lines_port.loc[ident, ["neighbor", "match_distance"]] = [name, True]
        # For buses which are not close to any existing bus, only add a new bus if the line is going offshore (e.g. North Sea Wind Power Hub)
        if not lines_port.match_distance.all() and offshore_shapes.union_all():
            potential_new_buses = lines_port[~lines_port.match_distance]
            is_offshore = potential_new_buses.apply(
                lambda x: offshore_shapes.union_all().contains(Point(x.x, x.y)), axis=1
            )
            new_buses = potential_new_buses[is_offshore]
            if not new_buses.empty:
                new_port, new_buses = add_new_buses(n, new_buses, suffix=f" bus{port}")
                new_buses["country"] = new_buses.apply(
                    lambda bus: find_country_for_bus(bus, offshore_shapes), axis=1
                )
                lines_port.loc[new_port.index, "match_distance"] = True
                lines_port.loc[new_port.index, "neighbor"] = new_port.neighbor
                new_buses_df = pd.concat([new_buses_df, new_buses])

        if not lines_port.match_distance.all():
            logger.warning(
                "Could not find bus close enough to connect the the following lines:\n"
                + str(lines_port[~lines_port.match_distance].index.to_list())
                + "\n Lines will be ignored."
            )
            lines.drop(lines_port[~lines_port.match_distance].index, inplace=True)
            lines_port = lines_port[lines_port.match_distance]

        lines.loc[lines_port.index, f"bus{port}"] = lines_port["neighbor"]

    lines = lines.assign(under_construction=True)

    return lines, new_buses_df


def get_branch_coords_from_geometry(linestring, reversed=False):
    """
    Reduces a linestring to its start and end points. Used to simplify the
    linestring which can have more than two points.

    Parameters
    ----------
    linestring : Shapely linestring
        The linestring to reduce.
    reversed : bool, optional
        If True, returns the end and start points instead of the start and end points.
        Defaults to False.

    Returns
    -------
    numpy.ndarray
        Flattened array of start and end coordinates.
    """
    coords = np.asarray(linestring.coords)
    ind = [0, -1] if not reversed else [-1, 0]
    start_end_coords = coords[ind]
    return start_end_coords.flatten()


def filter_junction_duplicates(duplicates, projects, existing):
    """
    Reviewed project junctions require matching topology for deduplication.

    A short new cable beside an existing crossing is not a duplicate merely
    because all endpoints lie within the geometric search radius.
    """
    keep = []
    for project, asset in duplicates.items():
        row = projects.loc[project]
        explicit_junction = any(
            pd.notna(row.get(f"junction{port}"))
            and str(row.get(f"junction{port}")).strip()
            for port in (0, 1)
        )
        if not explicit_junction or {row.bus0, row.bus1} == {
            existing.at[asset, "bus0"],
            existing.at[asset, "bus1"],
        }:
            keep.append(project)
    return duplicates.loc[keep]


def filter_parallel_duplicates(duplicates, projects, existing):
    """
    Reviewed parallel new builds are duplicates only if the corridor is complete.

    The corridor comprises all existing lines on the matched line's bus pair and
    voltage. With the reviewed circuit count before the project, the project is
    built; with the count after it, it is already in the network. Any other count
    keeps the duplicate and is logged for review.
    """
    if "parallel_circuits" not in projects:
        return duplicates
    build = []
    for project, asset in duplicates.items():
        check = projects.at[project, "parallel_circuits"]
        if pd.isna(check):
            continue
        check = json.loads(check)
        pair = {existing.at[asset, "bus0"], existing.at[asset, "bus1"]}
        corridor = existing[
            (existing.v_nom == existing.at[asset, "v_nom"])
            & existing.apply(lambda line: {line.bus0, line.bus1} == pair, axis=1)
        ]
        circuits = corridor.num_parallel.sum()
        if circuits == check["before"]:
            build.append(project)
        elif circuits != check["after"]:
            logger.warning(
                f"Parallel circuit check mismatch for {project}: corridor {sorted(corridor.index)} "
                f"has {circuits} circuits, expected {check['before']} before or {check['after']} after; "
                "kept as duplicate for review"
            )
    return duplicates.drop(build)


def get_branch_coords_from_buses(line):
    """
    Gets line string for branch component in an pypsa network.

    Parameters
    ----------
    line : pandas.Series
        A row from a branch component DataFrame with bus0 and bus1 attributes.

    Returns
    -------
    numpy.ndarray: Flattened array of start and end coordinates.
    """
    start_coords = n.buses.loc[line.bus0, ["x", "y"]].values
    end_coords = n.buses.loc[line.bus1, ["x", "y"]].values
    return np.array([start_coords, end_coords]).flatten()


def get_branch_coords(line):
    """
    Get start and end coordinates of a network branch.

    Project branches added earlier keep their source geometry, so a later
    plan's project is compared with the same source coordinates instead of
    the bus positions the earlier project was snapped to.

    Parameters
    ----------
    line : pandas.Series
        A row from a branch component DataFrame.

    Returns
    -------
    numpy.ndarray: Flattened array of start and end coordinates.
    """
    if isinstance(line.get("geometry"), LineString):
        return get_branch_coords_from_geometry(line["geometry"])
    return get_branch_coords_from_buses(line)


def get_bus_coords_from_port(linestring, port=0):
    """
    Extracts the coordinates of a specified port from a given linestring.

    Parameters
    ----------
    linestring: The shapely linestring.
    port (int): The index of the port to extract coordinates from. Default is 0.

    Returns
    -------
    tuple: The coordinates of the specified port as a tuple (x, y).
    """
    coords = np.asarray(linestring.coords)
    ind = [0, -1]
    coords = coords[ind]
    coords = coords[port]
    return coords


def find_closest_lines(lines, new_lines, distance_upper_bound=0.1, type="new"):
    """
    Find the closest lines in the existing set of lines to a set of new lines.

    Parameters
    ----------
    lines (pandas.DataFrame): DataFrame of the existing lines.
    new_lines (pandas.DataFrame): DataFrame with column geometry containing the new lines.
    distance_upper_bound (float, optional): Maximum distance to consider a line as a match. Defaults to 0.1 which corresponds to approximately 15 km.

    Returns
    -------
    pandas.Series: Series containing with index the new lines and values providing closest existing line.
    """

    # get coordinates of start and end points of all lines, for new lines we need to check both directions
    treelines = lines.apply(get_branch_coords, axis=1)
    querylines = pd.concat(
        [
            new_lines["geometry"].apply(get_branch_coords_from_geometry),
            new_lines["geometry"].apply(get_branch_coords_from_geometry, reversed=True),
        ]
    )
    treelines = np.vstack(treelines)
    querylines = np.vstack(querylines)
    tree = spatial.KDTree(treelines)
    dist, ind = tree.query(querylines, distance_upper_bound=distance_upper_bound)
    found_b = ind < len(lines)
    # since the new lines are checked in both directions, we need to find the correct index of the new line
    found_i = np.arange(len(querylines))[found_b] % len(new_lines)
    # create a DataFrame with the distances, new line and its closest existing line
    line_map = pd.DataFrame(
        dict(D=dist[found_b], existing_line=lines.index[ind[found_b] % len(lines)]),
        index=new_lines.index[found_i].rename("new_lines"),
    )
    if type == "new":
        if len(found_i) != 0:
            # compare if attribute of new line and existing line is similar
            attr = "p_nom" if "p_nom" in lines else "v_nom"
            # potential duplicates
            duplicated = line_map["existing_line"]
            # only if lines are similar in terms of p_nom or v_nom they are kept as duplicates
            to_keep = is_similar(
                new_lines.loc[duplicated.index, attr],
                duplicated.map(lines[attr]),
                percentage=10,
            )
            line_map = line_map[to_keep]
            if not line_map.empty:
                logger.warning(
                    "Found new lines similar to existing lines:\n"
                    + str(line_map["existing_line"].to_dict())
                    + "\n Lines are assumed to be duplicated and will be ignored."
                )
    # only keep the closer line of the new line pair (since lines are checked in both directions)
    # equidistant parallel branches (one stored reversed) go to the first name
    line_map = line_map.sort_values(by=["D", "existing_line"])[
        lambda ds: ~ds.index.duplicated(keep="first")
    ].sort_index()["existing_line"]
    if type == "upgraded":
        line_map = add_branch_chains(lines, new_lines, line_map, distance_upper_bound)
    return line_map


def find_branch_chain(lines, linestring, distance_upper_bound):
    """
    Find the shortest chain of existing branches between the ends of a line.

    The base network splits corridors at intermediate buses, so an upgrade
    may span several branches. Branches of one voltage level form a graph;
    the chain may be at most twice as long as the straight distance.

    Parameters
    ----------
    lines : pandas.DataFrame
        Existing branches with bus0, bus1 and length (km).
    linestring : shapely.geometry.LineString
        Geometry of the upgraded line.
    distance_upper_bound : float
        Maximum distance between a line end and its nearest bus in degrees.

    Returns
    -------
    list
        Names of the chained branches; empty if no chain is found.
    """
    ends = get_branch_coords_from_geometry(linestring).reshape(2, 2)
    straight = gpd.GeoSeries([linestring], crs=4326).to_crs(3035).length.iloc[0]
    chain, limit = [], 2 * straight / 1e3
    level = lines["v_nom"] if "v_nom" in lines else pd.Series(0, lines.index)
    for _, group in lines.groupby(level):
        buses = pd.Index(pd.unique(group[["bus0", "bus1"]].values.ravel()))
        dist, ind = spatial.KDTree(n.buses.loc[buses, ["x", "y"]]).query(ends)
        if (dist > distance_upper_bound).any() or ind[0] == ind[1]:
            continue
        i, j = buses.get_indexer(group.bus0), buses.get_indexer(group.bus1)
        edges = group.assign(u=np.minimum(i, j), v=np.maximum(i, j))
        edges = edges.sort_values("length").drop_duplicates(["u", "v"])
        weight = edges.length.clip(lower=1e-3)
        graph = sparse.coo_matrix((weight, (edges.u, edges.v)), shape=(len(buses),) * 2)
        length, pred = dijkstra(
            graph, directed=False, indices=ind[0], return_predecessors=True
        )
        if length[ind[1]] > limit:
            continue
        limit, chain, k = length[ind[1]], [], ind[1]
        names = edges.reset_index().set_index(["u", "v"]).iloc[:, 0]
        while k != ind[0]:
            chain.append(names[(min(k, pred[k]), max(k, pred[k]))])
            k = pred[k]
    return chain


def add_branch_chains(lines, new_lines, line_map, distance_upper_bound):
    """
    Map upgrades that match no single existing branch to a chain of branches.

    Parameters
    ----------
    lines : pandas.DataFrame
        Existing branches.
    new_lines : pandas.DataFrame
        Upgrades with column geometry.
    line_map : pandas.Series
        Upgrades matched to a single existing branch.
    distance_upper_bound : float
        Maximum distance between a line end and its nearest bus in degrees.

    Returns
    -------
    pandas.Series
        line_map with one entry per chained branch; unmatched upgrades are
        logged and left out.
    """
    not_found = new_lines.index.difference(line_map.index)
    chains = pd.Series(
        [
            find_branch_chain(lines, new_lines.at[i, "geometry"], distance_upper_bound)
            for i in not_found
        ],
        index=not_found,
        dtype=object,
    )
    found = chains.map(len) > 0
    if found.any():
        logger.info(
            "Upgraded lines matched to a chain of existing lines:\n"
            + str(chains[found].to_dict())
        )
    if not found.all():
        logger.warning(
            "Could not find upgraded lines close enough to existing lines:\n"
            + str(chains.index[~found].to_list())
            + "\n Lines will be ignored."
        )
    return pd.concat([line_map, chains[found].explode()])


def check_reviewed_branch_scopes(line_map, projects):
    """
    Check reviewed branch scopes against the match and drop excluded branches.

    Every reviewed branch must still be matched. Branches with the action
    ``exclude`` lie outside the project, e.g. beyond its end substation.
    Branches with the action ``retire`` are old lines the project removes
    without successor, e.g. when the base already holds the new line beside it.

    Returns
    -------
    tuple of pandas.Series
        line_map without the excluded and retired branches, and the retired
        branches indexed by project.
    """
    retired = pd.Series(dtype=object)
    if "branch_scope_rules" not in projects:
        return line_map, retired
    keep = pd.Series(True, index=range(len(line_map)))
    for ident, raw in projects["branch_scope_rules"].dropna().items():
        rules = json.loads(raw)
        own = line_map.index == ident
        missing = set(rules) - set(line_map[own])
        if missing:
            raise ValueError(
                f"Reviewed branch scope no longer matches {ident}: {missing}"
            )
        dropped = {
            asset
            for asset, rule in rules.items()
            if rule["action"] in ("exclude", "retire")
        }
        keep &= ~(own & line_map.isin(dropped).to_numpy())
        gone = [asset for asset, rule in rules.items() if rule["action"] == "retire"]
        retired = pd.concat(
            [retired, pd.Series(gone, index=[ident] * len(gone), dtype=object)]
        )
    return line_map[keep.to_numpy()], retired


def retirements(retired, projects):
    """Decommission branches that projects retire without successor in their build year."""
    out = pd.DataFrame(index=pd.Index(retired.to_numpy(), dtype=object))
    out["build_year"] = 1990
    out["lifetime"] = projects.loc[retired.index, "build_year"].to_numpy() - 1990
    out["retired_by"] = retired.index.to_numpy()
    return out


def expand_corridor_upgrades(line_map, candidates, projects):
    """
    Extend each upgrade to all branches on its matched bus pairs and voltages.

    OpenStreetMap often splits one multi-circuit line into parallel branches.
    An upgrade of the line then affects every branch of the corridor. Each
    branch keeps its own circuit count when the project gives none, or gives the
    matched branch's or the corridor's count. Reviewed partial scopes and
    unexplained circuit counts are not expanded; the latter are logged.

    Returns
    -------
    tuple of pandas.Series
        Extended line_map and, aligned with it, whether a row inherits the
        circuit count of its own branch.
    """
    pairs = candidates.apply(
        lambda line: (frozenset((line.bus0, line.bus1)), line.v_nom), axis=1
    )
    projects_out, assets_out, inherit_out = [], [], []
    for project, matched in line_map.groupby(level=0, sort=False):
        matched = matched.tolist()
        circuits = projects.at[project, "num_parallel"]
        reviewed = "branch_scope_rules" in projects and pd.notna(
            projects.at[project, "branch_scope_rules"]
        )
        siblings = (
            []
            if reviewed
            else [
                asset
                for asset in candidates.index[pairs.isin(set(pairs[matched]))]
                if asset not in matched
            ]
        )
        totals = {
            pairs[asset]: candidates.num_parallel[pairs == pairs[asset]].sum()
            for asset in matched
        }
        fits = pd.isna(circuits) or all(
            circuits in (candidates.at[asset, "num_parallel"], totals[pairs[asset]])
            for asset in matched
        )
        if siblings and not fits:
            logger.warning(
                f"Corridor upgrade not expanded for {project}: {circuits} circuits fit neither "
                f"{matched} nor its corridor; parallel branches {siblings} stay unchanged for review"
            )
            siblings = []
        assets = matched + siblings
        projects_out += [project] * len(assets)
        assets_out += assets
        inherit_out += [bool(siblings)] * len(assets)
    index = pd.Index(projects_out, name=line_map.index.name)
    return pd.Series(assets_out, index=index), pd.Series(inherit_out, index=index)


def add_line_types(n, path):
    """
    Add reviewed line types that the PyPSA catalogue lacks.

    Parameters
    ----------
    n : pypsa.Network
        Network whose ``line_types`` are extended in place.
    path : str or Path
        CSV with the PyPSA ``line_types`` columns, indexed by type name.

    Returns
    -------
    pandas.DataFrame
        The added types, for downstream networks that must register them too.
    """
    types = pd.read_csv(path, index_col=0)
    clash = types.index.intersection(n.line_types.index)
    if not clash.empty:
        raise ValueError(
            f"Reviewed line types already exist in the catalogue: {clash.tolist()}"
        )
    n.add("LineType", types.index, **types)
    return types


def resolve_explicit_connections(rows, buses, candidates):
    """
    Build reviewed explicit connections instead of geographic matching.

    Each component connects given model buses with the reviewed countries. It
    either adds a branch or replaces the named existing branches, which must
    span the same countries. Reviewed fields override the row; ``None`` clears
    a field (e.g. length inherited from the replaced branch).

    Returns
    -------
    tuple
        Added branches, replacing branches (integer index) and, aligned with
        them, the replaced existing branches.
    """
    added, upgraded, replaced = [], [], []
    for ident, row in rows.iterrows():
        parts = json.loads(row["explicit_connection"])["components"]
        for part in parts:
            ends = [part["bus0"], part["bus1"]]
            if not set(ends) <= set(buses.index):
                raise ValueError(f"Explicit connection {ident}: unknown bus in {ends}")
            countries = buses.loc[ends, "country"].tolist()
            if countries != part["countries"] or "DE" not in countries:
                raise ValueError(
                    f"Explicit connection {ident}: bus countries {countries}, reviewed {part['countries']}"
                )
            item = row.copy()
            item["bus0"], item["bus1"] = ends
            item["connection_part"] = part["part"]
            for field in ("v_nom", "num_parallel", "length", "type"):
                if field in part:
                    item[field] = np.nan if part[field] is None else part[field]
            for asset in part["replaces"]:
                if asset not in candidates.index:
                    raise ValueError(
                        f"Explicit connection {ident}: unknown branch {asset}"
                    )
                spans = sorted(
                    buses.loc[candidates.loc[asset, ["bus0", "bus1"]], "country"]
                )
                if spans != sorted(countries):
                    raise ValueError(
                        f"Explicit connection {ident}: {asset} does not span {countries}"
                    )
                upgraded.append(item)
                replaced.append(asset)
            if not part["replaces"]:
                item.name = ident if len(parts) == 1 else f"{ident}_{part['part']}"
                added.append(item)
    return (
        pd.DataFrame(added),
        pd.DataFrame(upgraded).reset_index(drop=True),
        pd.Series(replaced, dtype=object),
    )


SPLIT_COLUMNS = [
    "bus0",
    "bus1",
    "type",
    "s_nom",
    "s_max_pu",
    "length",
    "num_parallel",
    "v_nom",
    "i_nom",
    "dc",
    "underground",
    "under_construction",
    "carrier",
    "build_year",
    "lifetime",
    "geometry",
]


def split_branches(n, frames, max_offset=0.05):
    """
    Split base branches at the reviewed junctions named in explicit connections.

    A junction marks a station that a base branch passes without a bus. The
    branch is replaced by consecutive pieces from bus0, with length shares
    from the projection of each junction onto its geometry; all other values
    are kept, so the split alone changes no capacity. Splits apply whatever
    the status of the naming project, so every selection sees one topology.

    Parameters
    ----------
    n : pypsa.Network
        Base network, changed in place.
    frames : list of pandas.DataFrame
        All project rows of the included plans.
    max_offset : float
        Largest distance (degrees) between a junction and its branch.

    Returns
    -------
    tuple of pandas.DataFrame
        Pieces, with the replaced branch in ``split_of``, and junction buses.
    """
    junctions = {}
    for frame in frames:
        for raw in frame.get("explicit_connection", pd.Series(dtype=object)).dropna():
            for spec in json.loads(raw).get("splits", []):
                if junctions.setdefault(spec["junction"], spec) != spec:
                    raise ValueError(f"Inconsistent split junction {spec['junction']}")
    cuts = {}
    for spec in junctions.values():
        name = f"project_junction:{spec['junction']}"
        voltages = {n.lines.at[branch, "v_nom"] for branch in spec["branches"]}
        if len(voltages) != 1:
            raise ValueError(
                f"Split junction {name} joins branches of different voltage"
            )
        n.add(
            "Bus",
            name,
            x=spec["x"],
            y=spec["y"],
            v_nom=float(voltages.pop()),
            carrier="AC",
            country=spec["country"],
            symbol="substation",
            substation_lv=False,
            substation_off=False,
        )
        for branch in spec["branches"]:
            cuts.setdefault(branch, []).append((name, Point(spec["x"], spec["y"])))
    pieces = []
    for branch, points in sorted(cuts.items()):
        line = n.lines.loc[branch]
        geometry = wkt.loads(line.geometry)
        if (
            Point(geometry.coords[0]).distance(
                Point(n.buses.loc[line.bus0, ["x", "y"]])
            )
            > 1e-6
        ):
            geometry = geometry.reverse()
        if any(geometry.distance(point) > max_offset for _, point in points):
            raise ValueError(f"Split junction too far from {branch}")
        shares = sorted(
            (geometry.project(point, normalized=True), name) for name, point in points
        )
        bounds = [0.0, *(share for share, _ in shares), 1.0]
        ends = [line.bus0, *(name for _, name in shares), line.bus1]
        if min(np.diff(bounds)) <= 0:
            raise ValueError(f"Split junctions do not cut {branch} into pieces")
        for i in range(len(ends) - 1):
            piece = line.reindex(SPLIT_COLUMNS)
            piece["bus0"], piece["bus1"] = ends[i], ends[i + 1]
            piece["length"] = line.length * (bounds[i + 1] - bounds[i])
            piece["geometry"] = ops.substring(
                geometry, bounds[i], bounds[i + 1], normalized=True
            ).wkt
            piece["split_of"] = branch
            pieces.append(piece.rename(f"{branch}#{i}"))
    pieces = pd.DataFrame(pieces)
    if not pieces.empty:
        n.remove("Line", list(cuts))
        n.add("Line", pieces.index, **pieces.drop(columns="split_of"))
        logger.info(f"Split branches at reviewed junctions: {sorted(cuts)}")
    buses = n.buses.loc[[f"project_junction:{j}" for j in junctions]].copy()
    return pieces, buses


def adjust_decommissioning(upgraded_lines, line_map):
    """
    Adjust the decommissioning year of the existing lines to the built year of
    the upgraded lines.
    """
    to_update = pd.DataFrame(index=line_map)
    to_update["build_year"] = (
        1990  # dummy build_year to make existing lines decommissioned when upgraded lines are built
    )
    to_update["lifetime"] = (
        upgraded_lines.rename(line_map)["build_year"] - 1990
    )  # set lifetime to the difference between build year of upgraded line and existing line
    return to_update


def get_upgraded_lines(branch_component, n, upgraded_lines, line_map):
    """
    Get upgraded lines by merging information of existing line and upgraded
    line.
    """
    upgraded_lines = upgraded_lines.copy()
    explicit = {}
    if branch_component == "Line":
        # Remember project intent before filling from the original base asset.
        # Later stages must distinguish a specified conductor from a fallback.
        for field in ("type", "v_nom", "num_parallel"):
            if field in upgraded_lines:
                upgraded_lines[field] = upgraded_lines[field].replace(
                    r"^\s*$", np.nan, regex=True
                )
                explicit[field] = upgraded_lines[field].notna().to_numpy()
            else:
                explicit[field] = np.zeros(len(upgraded_lines), dtype=bool)
    # existing branch per project row; several projects may share one branch
    existing = n.df(branch_component).loc[line_map].set_axis(line_map.index)
    # the line type of the existing branch does not fit a new voltage level
    if "type" in existing and "v_nom" in upgraded_lines:
        new_voltage = upgraded_lines["v_nom"].notna() & upgraded_lines["v_nom"].ne(
            existing["v_nom"]
        )
        same_380_400_level = upgraded_lines["v_nom"].isin([380, 400]) & existing[
            "v_nom"
        ].isin([380, 400])
        new_voltage &= ~same_380_400_level
        existing.loc[new_voltage, "type"] = np.nan
    # project values take precedence, gaps are filled from the existing branch
    columns = list(dict.fromkeys(["bus0", "bus1", *upgraded_lines.columns]))
    if branch_component == "Line":
        columns += [
            field for field in explicit if field not in columns and field in existing
        ]
    lines_to_add = upgraded_lines.combine_first(existing)[columns]
    for field, provided in explicit.items():
        lines_to_add[f"_project_{field}_explicit"] = provided
    if branch_component == "Line":
        lines_to_add["matched_base_asset"] = line_map.to_numpy()
    # fall back to the existing asset's capacity only when the project gives
    # none (lines carry no s_nom - it is recomputed from the line type later;
    # links carry their upgraded p_nom and must keep it)
    attr = nominal_attrs[branch_component]
    if attr not in upgraded_lines.columns:
        lines_to_add[attr] = existing[attr]
    # change index of new lines to avoid duplicates
    lines_to_add.index = line_map.astype(str).values + "_upgraded"
    return lines_to_add


def stage_upgrades(new_df, adjust_df):
    """
    Turn several upgrades of one existing branch into successive stages.

    Upgrades from all plans are matched against the same existing branches.
    Per branch, upgrades are ordered by build year; each stage replaces the
    previous one in its build year, so capacity is never counted twice. For
    equal build years only the first upgrade (plan include order) is kept.
    Undated upgrades are dropped if the branch has a dated upgrade.

    Parameters
    ----------
    new_df : pandas.DataFrame
        New and upgraded branches; upgrades are indexed ``<asset>_upgraded``
        and carry the original branch in ``replaced_asset``.
    adjust_df : pandas.DataFrame
        Retirement adjustments of the original branches.

    Returns
    -------
    tuple of pandas.DataFrame
        Staged branches and one retirement adjustment per original branch.
    """
    if new_df.empty or "replaced_asset" not in new_df:
        return new_df, adjust_df
    is_upgrade = new_df.replaced_asset.notna()
    ups = new_df[is_upgrade].sort_values("build_year", kind="stable")
    asset = ups.replaced_asset
    drop = ups.duplicated(["replaced_asset", "build_year"]) | (
        ups.build_year.isna() & asset.duplicated()
    )
    if drop.any():
        logger.warning(
            "Upgrades of an already upgraded branch in the same (or no) build "
            "year are ignored:\n"
            + str(ups.loc[drop, ["project_id", "replaced_asset"]].values.tolist())
        )
    ups = ups[~drop]
    asset = ups.replaced_asset.to_numpy()
    step = ups.groupby(asset).cumcount().to_numpy()
    names = pd.Series(asset + "_upgraded") + [str(k + 1) if k else "" for k in step]
    previous = names.groupby(asset).shift(1).to_numpy()
    next_year = ups.groupby(asset).build_year.shift(-1).to_numpy()
    if "_project_type_explicit" in ups:
        ups["line_type_origin"] = ""
        ups["line_type_inherited_from"] = ""
        col = {field: ups.columns.get_loc(field) for field in ups.columns}
        inherited = []
        for original, positions in ups.groupby(
            "replaced_asset", sort=False
        ).indices.items():
            preceding = None
            for pos in positions:
                # Unspecified electrical fields retain the preceding stage's
                # state. Explicit changes, including circuit scope, are kept.
                if preceding is not None:
                    for field in ("v_nom", "num_parallel"):
                        if not ups.iat[pos, col[f"_project_{field}_explicit"]]:
                            ups.iat[pos, col[field]] = ups.iat[preceding, col[field]]
                if ups.iat[pos, col["_project_type_explicit"]]:
                    origin, parent = "project", ""
                elif preceding is None:
                    origin, parent = "base_network", original
                else:
                    voltage = ups.iat[pos, col["v_nom"]]
                    old_voltage = ups.iat[preceding, col["v_nom"]]
                    same_level = (
                        pd.notna(voltage)
                        and pd.notna(old_voltage)
                        and (
                            voltage == old_voltage
                            or voltage in (380, 400)
                            and old_voltage in (380, 400)
                        )
                    )
                    if same_level:
                        ups.iat[pos, col["type"]] = ups.iat[preceding, col["type"]]
                        origin, parent = "previous_stage", names.iloc[preceding]
                        inherited.append((ups.iat[pos, col["project_id"]], parent))
                    else:
                        # A conductor for another voltage is not a valid
                        # inherited design. Keep the existing fallback policy.
                        ups.iat[pos, col["type"]] = np.nan
                        origin, parent = "unspecified_voltage_change", ""
                ups.iat[pos, col["line_type_origin"]] = origin
                ups.iat[pos, col["line_type_inherited_from"]] = parent
                preceding = pos
        if inherited:
            logger.info(
                "Line types inherited from preceding upgrade stages:\n%s", inherited
            )
    ups = ups.set_axis(names.to_numpy())
    ups["replaced_asset"] = np.where(step == 0, asset, previous)
    ups["lifetime"] = np.nan_to_num(next_year - ups.build_year, nan=np.inf)
    staged = ups.replaced_asset.str.contains("_upgraded")
    if staged.any():
        logger.info(
            "Upgrades staged on previous upgrades:\n"
            + str(ups.loc[staged, ["project_id", "replaced_asset"]].values.tolist())
        )
    new_df = pd.concat([new_df[~is_upgrade], ups])
    new_df["lifetime"] = new_df["lifetime"].fillna(np.inf)
    adjust_df = adjust_df.groupby(level=0).min()
    return new_df, adjust_df


def apply_reviewed_branch_scopes(lines, base_lines):
    """
    Apply opt-in, asset-specific source reviews to terminal AC stages.

    Partial upgrades replace an aggregate with an upgraded circuit and a
    separately retained circuit. Both use the unchanged catalogue rating path.
    Explicit guards reject stale reviews if matching or preceding design changes.
    """
    if "branch_scope_rules" not in lines:
        return lines
    lines = lines.copy()
    retained = []
    for name, row in lines.iterrows():
        raw = row.get("branch_scope_rules")
        if pd.isna(raw) or not str(raw).strip():
            continue
        rules = json.loads(raw)
        rule = rules.get(row.get("matched_base_asset"))
        if rule is None:
            continue
        if lines.get("replaced_asset", pd.Series(dtype=object)).eq(name).any():
            raise ValueError(
                f"Reviewed circuit scope must be rechecked for a later stage: {name}"
            )
        previous_name = row["replaced_asset"]
        previous = (
            lines.loc[previous_name]
            if previous_name in lines.index
            else base_lines.loc[previous_name]
        )
        for field, expected_key in (
            ("num_parallel", "expected_previous_num_parallel"),
            ("v_nom", "expected_previous_voltage"),
            ("type", "expected_previous_type"),
        ):
            actual, expected = previous[field], rule[expected_key]
            # The normal catalogue fallback runs after staging. Compare its
            # effective value when an earlier voltage conversion cleared type.
            if field == "type" and (pd.isna(actual) or not str(actual).strip()):
                actual = "Al/St 240/40 4-bundle 380.0"
            if (
                str(actual) != str(expected)
                if field == "type"
                else not np.isclose(float(actual), float(expected))
            ):
                raise ValueError(
                    f"Stale circuit-scope review for {name}: {field}={actual}, expected {expected}"
                )
        if rule["action"] not in ("replace", "partial"):
            raise ValueError(f"Unknown circuit-scope action: {rule['action']}")
        for field, value in rule["values"].items():
            if field not in ("num_parallel", "v_nom", "type"):
                raise ValueError(f"Unsupported scoped field {field}")
            lines.loc[name, field] = value
        lines.loc[name, "circuit_role"] = "upgraded"
        lines.loc[name, "scope_basis"] = rule["basis"]
        if rule["action"] == "partial":
            affected = float(lines.loc[name, "num_parallel"])
            total = float(previous["num_parallel"])
            if (
                not (0 < affected < total)
                or not affected.is_integer()
                or not total.is_integer()
            ):
                raise ValueError(f"Invalid affected/retained circuit count for {name}")
            if float(row.v_nom) != float(previous.v_nom):
                raise ValueError(
                    f"Partial circuit review cannot change voltage class: {name}"
                )
            residual = lines.loc[name].copy()
            residual.name = f"{name}_retained"
            if residual.name in lines.index:
                raise ValueError(f"Duplicate retained branch {residual.name}")
            for field in ("type", "v_nom", "underground"):
                if field in previous:
                    residual[field] = previous[field]
            residual["num_parallel"] = total - affected
            residual["replaced_asset"] = np.nan
            residual["retained_from_asset"] = previous_name
            residual["scope_parent_stage"] = name
            residual["circuit_role"] = "retained"
            residual["line_type_origin"] = "retained_previous_circuit"
            residual["line_type_inherited_from"] = previous_name
            residual["branch_scope_rules"] = np.nan
            if pd.notna(residual.get("tags")):
                tags = json.loads(residual["tags"])
                tags["retained_circuit_scope"] = dict(
                    parent_stage=name,
                    previous_asset=previous_name,
                    circuits=total - affected,
                    basis=rule["basis"],
                )
                residual["tags"] = json.dumps(tags, ensure_ascii=False, sort_keys=True)
            retained.append(residual)
    if retained:
        lines = pd.concat([lines, pd.DataFrame(retained)])
    return lines


def deduplicate_new_asset_ids(branches):
    """
    Match the network import's first-included-plan precedence for new IDs.

    Network.add already skips a repeated new asset ID. CSV concatenation must
    make the same choice; repeated pre-staging upgrade IDs are intentional.
    """
    is_new = (
        branches.get("replaced_asset", pd.Series(np.nan, index=branches.index))
        .isna()
        .to_numpy()
    )
    drop = np.zeros(len(branches), dtype=bool)
    drop[np.flatnonzero(is_new)] = branches.index[is_new].duplicated(keep="first")
    if drop.any():
        logger.warning(
            "Repeated new asset IDs: retaining first included plan, matching Network.add: %s",
            branches.index[drop].tolist(),
        )
    return branches.iloc[~drop]


def get_project_files(path, skip=[]):
    path = Path(path)
    lines = {}
    files = [
        p
        for p in path.iterdir()
        if p.is_file()
        and p.suffix == ".csv"
        and not any(substring in p.name for substring in skip)
    ]
    if not files:
        logger.warning(f"No projects found for {path.parent.name}")
        return lines
    for file in files:
        df = pd.read_csv(file, index_col=0)
        df["geometry"] = df.apply(
            lambda x: LineString([[x.x0, x.y0], [x.x1, x.y1]]), axis=1
        )
        df.drop(columns=["x0", "y0", "x1", "y1"], inplace=True)
        lines[file.stem] = df
    return lines


def remove_projects_outside_countries(lines, europe_shape):
    """
    Remove projects which are not in the considered countries.
    """
    europe_shape_prepped = shapely.prepared.prep(europe_shape)
    is_within_covered_countries = lines["geometry"].apply(
        lambda x: europe_shape_prepped.contains(x)
    )

    if not is_within_covered_countries.all():
        logger.warning(
            "Project lines outside of the covered area (skipping): "
            + ", ".join(str(i) for i in lines.loc[~is_within_covered_countries].index)
        )

    lines = lines.loc[is_within_covered_countries]
    return lines


def is_similar(ds1, ds2, percentage=10):
    """
    Check if values in series ds2 are within a specified percentage of series
    ds1.

    Returns:
    - A boolean series where True indicates ds2 values are within the percentage range of ds2.
    """
    lower_bound = ds1 * (1 - percentage / 100)
    upper_bound = ds1 * (1 + percentage / 100)
    return np.logical_and(ds2 >= lower_bound, ds2 <= upper_bound)


def set_underwater_fraction(new_links, offshore_shapes):
    new_links_gds = gpd.GeoSeries(new_links["geometry"])
    new_links["underwater_fraction"] = (
        new_links_gds.intersection(offshore_shapes.union_all()).length
        / new_links_gds.length
    ).round(2)


def add_projects(
    n,
    new_lines_df,
    new_links_df,
    adjust_lines_df,
    adjust_links_df,
    new_buses_df,
    europe_shape,
    offshore_shapes,
    path,
    plan,
    status=["confirmed", "under construction"],
    skip=[],
    upgrade_candidates=None,
):
    if upgrade_candidates is None:
        upgrade_candidates = {
            "Line": n.lines.copy(),
            "Link": n.links.query("carrier == 'DC'").copy(),
        }
    lines_dict = get_project_files(path, skip=skip)
    for key, lines in lines_dict.items():
        # Keep the source identity even when an upgrade is renamed to the
        # matched network branch. Used by the standalone project audit.
        lines["project_id"] = lines.index.astype(str)
        lines["source_file"] = str(Path(path) / f"{key}.csv")
        lines["project_action"] = key
        logger.info(f"Processing {key.replace('_', ' ')} projects from {plan}.")
        lines = remove_projects_outside_countries(lines, europe_shape)
        if isinstance(status, dict):
            status = status[plan]
        lines = lines.loc[lines.project_status.isin(status)]
        if "explicit_connection" in lines and lines.explicit_connection.notna().any():
            explicit = lines.explicit_connection.notna()
            added, replacing, line_map = resolve_explicit_connections(
                lines[explicit], n.buses, upgrade_candidates["Line"]
            )
            lines = lines[~explicit]
            logger.info(
                f"Explicit connections from {plan}: added {added.index.tolist()}, replaced {line_map.tolist()}"
            )
            if not added.empty:
                new_lines_df = pd.concat(
                    [new_lines_df, added.assign(under_construction=True)]
                )
                n.add("Line", added.index, **added.assign(under_construction=True))
            if not replacing.empty:
                replacing = replacing.assign(replaced_asset=line_map)
                adjust_lines_df = pd.concat(
                    [adjust_lines_df, adjust_decommissioning(replacing, line_map)]
                )
                new_lines_df = pd.concat(
                    [new_lines_df, get_upgraded_lines("Line", n, replacing, line_map)]
                )
        if lines.empty:
            continue
        if key == "new_lines":
            new_lines, new_buses_df = connect_new_lines(
                lines, n, new_buses_df, bus_carrier="AC"
            )
            duplicate_lines = find_closest_lines(
                n.lines, new_lines, distance_upper_bound=0.10, type="new"
            )
            duplicate_lines = filter_junction_duplicates(
                duplicate_lines, new_lines, n.lines
            )
            duplicate_lines = filter_parallel_duplicates(
                duplicate_lines, new_lines, n.lines
            )
            new_lines = new_lines.drop(duplicate_lines.index, errors="ignore")
            new_lines_df = pd.concat([new_lines_df, new_lines])
            # add new lines to network to be able to find added duplicates
            n.add("Line", new_lines.index, **new_lines)
        elif key == "new_links":
            new_links, new_buses_df = connect_new_lines(
                lines,
                n,
                new_buses_df,
                offshore_shapes=offshore_shapes,
                distance_upper_bound=0.4,
                bus_carrier=["AC", "DC"],
            )
            duplicate_links = find_closest_lines(
                n.links, new_links, distance_upper_bound=0.10, type="new"
            )
            new_links = new_links.drop(duplicate_links.index, errors="ignore")
            set_underwater_fraction(new_links, offshore_shapes)
            new_links_df = pd.concat([new_links_df, new_links])
            # add new links to network to be able to find added duplicates
            n.add("Link", new_links.index, **new_links)
        elif key == "upgraded_lines":
            line_map = find_closest_lines(
                upgrade_candidates["Line"],
                lines,
                distance_upper_bound=0.30,
                type="upgraded",
            )
            line_map, retired = check_reviewed_branch_scopes(line_map, lines)
            if not retired.empty:
                adjust_lines_df = pd.concat(
                    [adjust_lines_df, retirements(retired, lines)]
                )
            line_map, inherit = expand_corridor_upgrades(
                line_map, upgrade_candidates["Line"], lines
            )
            # chained upgrades repeat a project row once per existing branch
            chained = line_map.index.duplicated(keep=False)
            upgraded_lines = lines.loc[line_map.index].reset_index(drop=True)
            upgraded_lines.loc[chained, "length"] = np.nan
            upgraded_lines.loc[inherit.to_numpy(), "num_parallel"] = np.nan
            line_map = line_map.reset_index(drop=True)
            upgraded_lines = upgraded_lines.assign(replaced_asset=line_map)
            lines_to_adjust = adjust_decommissioning(upgraded_lines, line_map)
            adjust_lines_df = pd.concat([adjust_lines_df, lines_to_adjust])
            upgraded_lines = get_upgraded_lines("Line", n, upgraded_lines, line_map)
            new_lines_df = pd.concat([new_lines_df, upgraded_lines])
        elif key == "upgraded_links":
            line_map = find_closest_lines(
                upgrade_candidates["Link"],
                lines,
                distance_upper_bound=0.30,
                type="upgraded",
            )
            chained = line_map.index.duplicated(keep=False)
            upgraded_links = lines.loc[line_map.index].reset_index(drop=True)
            upgraded_links.loc[chained, "length"] = np.nan
            line_map = line_map.reset_index(drop=True)
            upgraded_links = upgraded_links.assign(replaced_asset=line_map)
            links_to_adjust = adjust_decommissioning(upgraded_links, line_map)
            adjust_links_df = pd.concat([adjust_links_df, links_to_adjust])
            upgraded_links = get_upgraded_lines("Link", n, upgraded_links, line_map)
            new_links_df = pd.concat([new_links_df, upgraded_links])
            set_underwater_fraction(new_links_df, offshore_shapes)
        else:
            logger.warning(f"Unknown project type {key}")
            continue
    return new_lines_df, new_links_df, adjust_lines_df, adjust_links_df, new_buses_df


def fill_length_from_geometry(line, line_factor=1.2):
    if not pd.isna(line.length):
        return line.length
    length = gpd.GeoSeries(line["geometry"], crs=4326).to_crs(3035).length.values[0]
    length = length / 1000 * line_factor
    return round(length, 1)


if __name__ == "__main__":
    if "snakemake" not in globals():
        from scripts._helpers import mock_snakemake

        snakemake = mock_snakemake("build_transmission_projects")
    configure_logging(snakemake)
    set_scenario_config(snakemake)

    line_factor = snakemake.params.line_factor
    s_max_pu = snakemake.params.s_max_pu

    n = pypsa.Network(snakemake.input.base_network)
    reviewed_types = add_line_types(n, snakemake.input.line_types)

    new_lines_df = pd.DataFrame()
    new_links_df = pd.DataFrame()
    adjust_lines_df = pd.DataFrame()
    adjust_links_df = pd.DataFrame()
    new_buses_df = pd.DataFrame()

    europe_shape = gpd.read_file(snakemake.input.europe_shape).loc[0, "geometry"]
    offshore_shapes = gpd.read_file(snakemake.input.offshore_shapes).rename(
        {"name": "country"}, axis=1
    )

    transmission_projects = snakemake.params.transmission_projects
    projects = [
        project
        for project, include in transmission_projects["include"].items()
        if include
    ]
    paths = snakemake.input.transmission_projects
    frames = [
        frame
        for project in projects
        for frame in get_project_files(
            select_project_path(project, paths), transmission_projects["skip"]
        ).values()
    ]
    split_lines_df, new_buses_df = split_branches(n, frames)
    upgrade_candidates = {
        "Line": n.lines.copy(),
        "Link": n.links.query("carrier == 'DC'").copy(),
    }
    for project in projects:
        path = select_project_path(project, paths)
        new_lines_df, new_links_df, adjust_lines_df, adjust_links_df, new_buses_df = (
            add_projects(
                n,
                new_lines_df,
                new_links_df,
                adjust_lines_df,
                adjust_links_df,
                new_buses_df,
                europe_shape,
                offshore_shapes,
                path=path,
                plan=project,
                status=transmission_projects["status"],
                skip=transmission_projects["skip"],
                upgrade_candidates=upgrade_candidates,
            )
        )
    new_lines_df = deduplicate_new_asset_ids(new_lines_df)
    new_links_df = deduplicate_new_asset_ids(new_links_df)
    # Resolve overlapping plan descriptions before interpreting dates as stages.
    review_path = transmission_projects.get("overlap_decisions")
    reviews = pd.read_csv(review_path, dtype=str).fillna("") if review_path else None
    overlap_audits = []
    for kind in ("Line", "Link"):
        frame, adjustments = (
            (new_lines_df, adjust_lines_df)
            if kind == "Line"
            else (new_links_df, adjust_links_df)
        )
        frame, adjustments, audit = reconcile(frame, adjustments, reviews)
        audit["component"] = kind
        overlap_audits.append(audit)
        for item in audit.to_dict("records"):
            if item["action"] == "shared_branch_precedence":
                logger.warning(
                    "Upgrade overlap resolved on shared model branch: %s",
                    json.dumps(item, sort_keys=True),
                )
        if kind == "Line":
            new_lines_df, adjust_lines_df = frame, adjustments
        else:
            new_links_df, adjust_links_df = frame, adjustments
    audit_path = Path(snakemake.output.new_lines).parent / "overlap_audit.csv"
    pd.concat(overlap_audits, ignore_index=True).to_csv(audit_path, index=False)
    new_lines_df, adjust_lines_df = stage_upgrades(new_lines_df, adjust_lines_df)
    new_links_df, adjust_links_df = stage_upgrades(new_links_df, adjust_links_df)
    new_lines_df = apply_reviewed_branch_scopes(new_lines_df, n.lines)
    if not new_lines_df.empty:
        line_type = "Al/St 240/40 4-bundle 380.0"
        missing_type = new_lines_df["type"].replace(r"^\s*$", np.nan, regex=True).isna()
        if "line_type_origin" not in new_lines_df:
            new_lines_df["line_type_origin"] = "project"
        new_lines_df["line_type_origin"] = new_lines_df["line_type_origin"].fillna(
            "project"
        )
        new_lines_df.loc[missing_type, "line_type_origin"] = "catalogue_default"
        # Add new line type for new lines
        new_lines_df["type"] = (
            new_lines_df["type"].replace(r"^\s*$", np.nan, regex=True).fillna(line_type)
        )
        new_lines_df["num_parallel"] = new_lines_df["num_parallel"].fillna(2)
        if "underground" in new_lines_df.columns:
            new_lines_df["underground"] = (
                new_lines_df["underground"].astype("bool").fillna(False)
            )
        # Add carrier types of lines
        new_lines_df["carrier"] = "AC"
        # Fill empty length values with length calculated from geometry
        new_lines_df["length"] = new_lines_df.apply(
            fill_length_from_geometry, args=(line_factor,), axis=1
        )
        # get s_nom from line type
        new_lines_df["s_nom"] = (
            np.sqrt(3)
            * n.line_types.loc[new_lines_df["type"], "i_nom"].values
            * new_lines_df["v_nom"]
            * new_lines_df["num_parallel"]
        ).round(2)
        # set s_max_pu
        new_lines_df["s_max_pu"] = s_max_pu
        new_lines_df.drop(
            columns=[
                f"_project_{field}_explicit"
                for field in ("type", "v_nom", "num_parallel")
            ],
            errors="ignore",
            inplace=True,
        )
    if not new_links_df.empty:
        # Add carrier types of lines and links
        new_links_df["carrier"] = "DC"
        # Fill empty length values with length calculated from geometry
        new_links_df["length"] = new_links_df.apply(
            fill_length_from_geometry, args=(line_factor,), axis=1
        )
        # Whether to keep existing link capacity or set to zero
        not_upgraded = ~new_links_df.index.str.contains("upgraded")
        if transmission_projects["new_link_capacity"] == "keep":
            new_links_df.loc[not_upgraded, "p_nom"] = new_links_df["p_nom"].fillna(0)
        elif transmission_projects["new_link_capacity"] == "zero":
            new_links_df.loc[not_upgraded, "p_nom"] = 0
    # export csv files for new buses, lines, links and adjusted lines and links
    new_lines_df = pd.concat([split_lines_df, new_lines_df])
    new_lines_df.to_csv(snakemake.output.new_lines)
    new_links_df.to_csv(snakemake.output.new_links)
    adjust_lines_df.to_csv(snakemake.output.adjust_lines)
    adjust_links_df.to_csv(snakemake.output.adjust_links)
    new_buses_df.to_csv(snakemake.output.new_buses)
    reviewed_types.to_csv(snakemake.output.line_types)
