# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2023 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Preprocess hydrogen kernnetz based on data from FNB Gas
(https://fnb-gas.de/wasserstoffnetz-wasserstoff-kernnetz/).
"""

import logging

logger = logging.getLogger(__name__)

import os
import sys
import uuid

import geopandas as gpd
import numpy as np
import pandas as pd
from pypsa.geo import haversine_pts
from shapely import wkt
from shapely.geometry import LineString, Point

paths = ["workflow/submodules/pypsa-eur/scripts", "../submodules/pypsa-eur/scripts"]
for path in paths:
    sys.path.insert(0, os.path.abspath(path))
from build_gas_network import diameter_to_capacity

MANUAL_ADDRESSES = {
    "Oude Statenzijl": (7.205108658430258, 53.20183834422634),
    "Helgoland": (7.882663327316698, 54.183393795580166),
    "SEN-1": (6.5, 55.0),
    "AWZ": (14.220711180456643, 54.429208831326804),
    "Bremen": (8.795818388451732, 53.077669699449594),
    "Bad Lauchstädt": (11.869106908389433, 51.38797498313352),
    "Großkugel": (12.151584743366769, 51.4166927585755),
    "Bobbau": (12.269345975889912, 51.69045938775995),
    "Visbeck": (8.310468203836264, 52.834518912466216),
    "Elbe-Süd": (9.608042769377906, 53.57422954537108),
    "Salzgitter": (10.386847343138689, 52.13861418123843),
    "Wefensleben": (11.15557835653467, 52.176005656180244),
    "Fessenheim": (7.5352027843079, 47.91300212650956),
    "Hittistetten": (10.09644829589717, 48.32667870548472),
    "Lindau": (9.690886766574819, 47.55387858107057),
    "Ludwigshafen": (8.444314472678961, 49.477207809634784),
    "Niederhohndorf": (12.466430165766688, 50.7532612203904),
    "Rückersdorf": (12.21941992347776, 50.822251899358236),
    "Bissingen": (10.6158383, 48.7177493),
    "Rehden": (8.476178919627396, 52.60675277527164),
    "Eynatten": (6.083339457526605, 50.69260916361823),
    "Vlieghuis": (6.8382504272201095, 52.66036497820981),
    "Kalle": (6.921180663621839, 52.573992586428425),
    "Carling": (6.713267207127634, 49.16738919353264),
    "Legden": (7.099754098013676, 52.03269789265483),
    "Ledgen": (7.099754098013676, 52.03269789265483),
    "Reiningen": (8.374879149975513, 52.50849502371421),
    "Buchholz": (12.929212986885771, 52.15737808332214),
    "Sandkrug": (8.257391972093515, 53.05387937393471),
}


def diameter_to_capacity_h2(pipe_diameter_mm):
    """
    Calculate pipe capacity in MW based on diameter in mm. Linear
    interpolation.

    20 inch (500 mm)  50 bar -> 1.2   GW H2 pipe capacity (LHV) 36 inch
    (900 mm)  50 bar -> 4.7   GW H2 pipe capacity (LHV) 48 inch (1200
    mm) 80 bar -> 16.9  GW H2 pipe capacity (LHV)

    Based on table 4 of
    https://ehb.eu/files/downloads/EHB-Analysing-the-future-demand-supply-and-transport-of-hydrogen-June-2021-v3.pdf
    """
    # slopes definitions
    m0 = (1200 - 0) / (500 - 0)
    m1 = (4700 - 1200) / (900 - 500)
    m2 = (16900 - 4700) / (1200 - 900)
    # intercepts
    a0 = 0
    a1 = 1200 - m1 * 500
    a2 = 4700 - m2 * 900

    if pipe_diameter_mm < 500:
        return a0 + m0 * pipe_diameter_mm
    elif pipe_diameter_mm < 900:
        return a1 + m1 * pipe_diameter_mm
    else:
        return a2 + m2 * pipe_diameter_mm


def load_and_merge_raw(fn1, fn2):
    # load, clean and merge

    df_fn1 = pd.read_excel(fn1, skiprows=2, skipfooter=2)
    df_fn2_retrofit = pd.read_excel(
        fn2, "Wasserstoff-Kernnetz Umstellung", skiprows=2, skipfooter=4
    )
    df_fn2_new = pd.read_excel(
        fn2, "Wasserstoff-Kernnetz Neubau", skiprows=3, skipfooter=2
    )

    for df in [df_fn1, df_fn2_retrofit, df_fn2_new]:
        df.columns = df.columns.str.replace("\n", "")

    # clean first dataset
    # drop lines not in Kernetz
    df_fn1 = df_fn1[df_fn1["Bestandteil des Wasserstoff-Kernnetzes"] == "ja"]

    to_keep = [
        "Name (Lfd.Nr.-Von-Nach)",
        "Umstellungsdatum/ Planerische Inbetriebnahme",
        "Anfangspunkt(Ort)",
        "Endpunkt(Ort)",
        "Nenndurchmesser (DN)",
        "Länge (km)",
        "Druckstufe (DP)[mind. 30 barg]",
        "Bundesland",
        "Umstellung/ Neubau",
        "IPCEI-Projekt(ja/ nein)",
    ]

    to_rename = {
        "Name (Lfd.Nr.-Von-Nach)": "name",
        "Umstellungsdatum/ Planerische Inbetriebnahme": "build_year",
        "Nenndurchmesser (DN)": "diameter_mm",
        "Länge (km)": "length",
        "Druckstufe (DP)[mind. 30 barg]": "max_pressure_bar",
        "Umstellung/ Neubau": "retrofitted",
        "IPCEI-Projekt(ja/ nein)": "ipcei",
    }

    df_fn1 = df_fn1[to_keep].rename(columns=to_rename)

    # extract info on retrofitted
    df_fn1["retrofitted"] = df_fn1.retrofitted != "Neubau"

    # clean second dataset
    # select only pipes
    df_fn2_new = df_fn2_new[df_fn2_new["Maßnahmenart"] == "Leitung"]

    to_keep = [
        "Name",
        "Planerische Inbetriebnahme",
        "Anfangspunkt(Ort)",
        "Endpunkt(Ort)",
        "Nenndurchmesser (DN)",
        "Länge (km)",
        "Druckstufe (DP)[mind. 30 barg]",
        "Bundesland",
        "retrofitted",
        "IPCEI-Projekt(Name/ nein)",
    ]

    to_rename = {
        "Name": "name",
        "Planerische Inbetriebnahme": "build_year",
        "Nenndurchmesser (DN)": "diameter_mm",
        "Länge (km)": "length",
        "Druckstufe (DP)[mind. 30 barg]": "max_pressure_bar",
        "IPCEI-Projekt(Name/ nein)": "ipcei",
    }

    df_fn2_new["retrofitted"] = False
    df_fn2_retrofit["retrofitted"] = True
    df_fn2 = pd.concat([df_fn2_retrofit, df_fn2_new])[to_keep].rename(columns=to_rename)
    df = pd.concat([df_fn1, df_fn2])
    df.reset_index(drop=True, inplace=True)

    return df


def prepare_dataset(df):

    # clean length
    df.length = df.length.astype(float)

    # clean diameter
    df.diameter_mm = (
        df.diameter_mm.astype(str)
        .str.extractall(r"(\d+)")
        .groupby(level=0)
        .last()
        .astype(int)
    )

    # clean max pressure
    df.max_pressure_bar = (
        df.max_pressure_bar.astype(str)
        .str.extractall(r"(\d+[.,]?\d*)")
        .groupby(level=0)
        .last()
        .squeeze()
        .str.replace(",", ".")
        .astype(float)
    )

    # clean build_year
    df.build_year = (
        df.build_year.astype(str).str.extract(r"(\b\d{4}\b)").astype(float).fillna(2032)
    )

    # create bidirectional and set true
    df["bidirectional"] = True

    df[["BL1", "BL2"]] = (
        df["Bundesland"]
        .apply(lambda bl: [bl.split("/")[0].strip(), bl.split("/")[-1].strip()])
        .apply(pd.Series)
    )

    # calc capa
    df["p_nom"] = df.diameter_mm.apply(diameter_to_capacity_h2)

    # eliminated gas capacity from retrofitted pipes
    df["removed_gas_cap"] = df.diameter_mm.apply(diameter_to_capacity)
    df[df.retrofitted == False]["removed_gas_cap"] == 0

    # eliminate leading and trailing spaces
    df["Anfangspunkt(Ort)"] = df["Anfangspunkt(Ort)"].str.strip()
    df["Endpunkt(Ort)"] = df["Endpunkt(Ort)"].str.strip()

    # drop pipes with same start and end
    df = df[df["Anfangspunkt(Ort)"] != df["Endpunkt(Ort)"]]

    # drop pipes with length smaller than 5 km
    df = df[df.length > 5]

    # reindex
    df.reset_index(drop=True, inplace=True)

    return df


def geocode_locations(df):

    try:
        from geopy.extra.rate_limiter import RateLimiter
        from geopy.geocoders import Nominatim
    except:
        raise ModuleNotFoundError(
            "Optional dependency 'geopy' not found."
            "Install via 'conda install -c conda-forge geopy'"
        )

    locator = Nominatim(user_agent=str(uuid.uuid4()))
    geocode = RateLimiter(locator.geocode, min_delay_seconds=2)
    # load state data for checking
    gdf_state = gpd.read_file(snakemake.input.gadm).set_index("GID_1")

    def get_location(row):
        def get_loc_A(loc="location", add_info=""):
            loc_A = Point(
                gpd.tools.geocode(row[loc] + ", " + add_info, timeout=7)["geometry"][0]
            )
            return loc_A

        def get_loc_B(loc="location", add_info=""):
            loc_B = geocode([row[loc], add_info], timeout=7)
            if loc_B is not None:
                loc_B = Point(loc_B.longitude, loc_B.latitude)
            else:
                loc_B = Point(0, 0)
            return loc_B

        def is_in_state(point, state="state"):
            if (row[state] in gdf_state.NAME_1.tolist()) & (point is not None):
                polygon_geometry = gdf_state[
                    gdf_state.NAME_1 == row[state]
                ].geometry.squeeze()
                return point.within(polygon_geometry)
            else:
                return False

        loc = get_loc_A("location", "Deutschland")

        # check if location is in Bundesland
        if not is_in_state(loc, "state"):
            # check if other loc is in Bundesland
            loc = get_loc_B("location", "Deutschland")
            # if both methods do not return loc in Bundesland, add Bundesland info
            if not is_in_state(loc, "state"):
                loc = get_loc_A("location", row["state"] + ", Deutschland")
                # if no location in Bundesland can be found
                if not is_in_state(loc, "state"):
                    loc = Point(0, 0)

        return loc

    # extract locations and state
    locations1, locations2 = (
        df[["Anfangspunkt(Ort)", "BL1"]],
        df[["Endpunkt(Ort)", "BL2"]],
    )
    locations1.columns, locations2.columns = ["location", "state"], [
        "location",
        "state",
    ]
    locations = pd.concat([locations1, locations2], axis=0)
    locations.drop_duplicates(inplace=True)

    # (3min)
    locations["point"] = locations.apply(lambda row: get_location(row), axis=1)

    # map manual locations (NOT FOUND OR WRONG)
    locations.point = locations.apply(
        lambda row: Point(
            MANUAL_ADDRESSES.get(row.location)
            if row.location in MANUAL_ADDRESSES.keys()
            else row.point
        ),
        axis=1,
    )

    return locations


def assign_locations(df, locations):

    df["point0"] = pd.merge(
        df,
        locations,
        left_on=["Anfangspunkt(Ort)", "BL1"],
        right_on=["location", "state"],
        how="left",
    )["point"]
    df["point1"] = pd.merge(
        df,
        locations,
        left_on=["Endpunkt(Ort)", "BL2"],
        right_on=["location", "state"],
        how="left",
    )["point"]

    # calc length of points
    length_factor = 1.0
    length_factor = 1.0
    df["length_haversine"] = df.apply(
        lambda p: length_factor
        * haversine_pts([p.point0.x, p.point0.y], [p.point1.x, p.point1.y]),
        axis=1,
    )

    # calc length ratio
    df["length_ratio"] = df.apply(
        lambda row: max(row.length, row.length_haversine)
        / (min(row.length, row.length_haversine) + 1),
        axis=1,
    )

    # only keep pipes with realistic length ratio
    df = df.query("retrofitted or length_ratio <= 2")

    # calc LineString
    df["geometry"] = df.apply(lambda x: LineString([x["point0"], x["point1"]]), axis=1)

    return df


if __name__ == "__main__":
    if "snakemake" not in globals():
        import os
        import sys

        path = "../submodules/pypsa-eur/scripts"
        sys.path.insert(0, os.path.abspath(path))
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_wasserstoff_kernnetz")

    logging.basicConfig(level=snakemake.config["logging"]["level"])

    wasserstoff_kernnetz = load_and_merge_raw(
        snakemake.input.wasserstoff_kernnetz_1,
        snakemake.input.wasserstoff_kernnetz_2,
    )

    wasserstoff_kernnetz = prepare_dataset(wasserstoff_kernnetz)

    if snakemake.params.reload_locations:
        locations = geocode_locations(wasserstoff_kernnetz)
    else:
        locations = pd.read_csv(snakemake.input.locations, index_col=0)
        locations["point"] = locations["point"].apply(wkt.loads)

    wasserstoff_kernnetz = assign_locations(wasserstoff_kernnetz, locations)

    wasserstoff_kernnetz.to_csv(snakemake.output.cleaned_wasserstoff_kernnetz)
