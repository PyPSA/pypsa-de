# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2024- The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Using BNetzA data to get a high resolution map of German CHP plants.

(https://open-mastr.readthedocs.io/en/latest/).
"""

import logging

logger = logging.getLogger(__name__)
import os
import sys

import geopandas as gpd
import numpy as np
import pandas as pd
import pypsa
from powerplantmatching.export import map_country_bus

from scripts._helpers import configure_logging, mock_snakemake


def clean_data(combustion, biomass, geodata):
    """
    Clean the data and return a dataframe with the relevant information.

    PLZ is translated to longitude and latitude using the pyGeoDb data.
    """
    biomass.dropna(subset="Postleitzahl", inplace=True)
    biomass.rename(
        columns={"NameStromerzeugungseinheit": "NameKraftwerk"}, inplace=True
    )
    biomass["Einsatzort"] = ""

    data = pd.concat([biomass, combustion], join="inner", ignore_index=True)

    data["IndustryStatus"] = data["Einsatzort"].str.contains("Industrie")
    data["IndustryStatus"] = data["IndustryStatus"].apply(
        lambda x: False if pd.isna(x) else x
    )

    # Get only CHP plants
    CHP_raw = data.query("ThermischeNutzleistung > 0").copy()
    CHP_raw.NameKraftwerk = CHP_raw.apply(
        lambda x: x.EinheitMastrNummer if pd.isna(x.NameKraftwerk) else x.NameKraftwerk,
        axis=1,
    )

    rename_columns = {
        "KwkMastrNummer": "ID",
        "NameKraftwerk": "Name",
        "Energietraeger": "Fueltype",
        "Technologie": "Technology",
        "ElektrischeKwkLeistung": "Capacity",
        "ThermischeNutzleistung": "Capacity_thermal",
        "Inbetriebnahmedatum": "DateIn",
        "DatumEndgueltigeStilllegung": "DateOut",
        "Postleitzahl": "Postleitzahl",
        "IndustryStatus": "Industry",
    }
    CHP_sel = CHP_raw[rename_columns.keys()].rename(columns=rename_columns)

    # change date format
    CHP_sel.DateIn = CHP_sel.DateIn.str[:4].astype(float)
    CHP_sel.DateOut = CHP_sel.DateOut.str[:4].astype(float)

    # delete duplicates identified by KwkMastrNummer
    strategies = {
        "Name": "first",
        "Fueltype": "first",
        "Technology": "first",
        "Capacity": "mean",  # dataset duplicates full KWK capacity for each block
        "Capacity_thermal": "mean",  # dataset duplicates full KWK capacity for each block
        "DateIn": "mean",
        "DateOut": "mean",
        "Postleitzahl": "first",
        "Industry": "first",
    }
    CHP_sel = CHP_sel.groupby("ID").agg(strategies).reset_index()

    # set missing information to match the powerplant data format
    CHP_sel[["Set", "Country", "Efficiency"]] = ["CHP", "DE", ""]
    CHP_sel[["lat", "lon"]] = [float("nan"), float("nan")]

    # get location from PLZ
    CHP_sel.fillna({"lat": CHP_sel.Postleitzahl.map(geodata.lat)}, inplace=True)
    CHP_sel.fillna({"lon": CHP_sel.Postleitzahl.map(geodata.lon)}, inplace=True)

    fueltype = {
        "Erdgas": "Natural Gas",
        "Mineralölprodukte": "Oil",
        "Steinkohle": "Coal",
        "Braunkohle": "Lignite",
        "andere Gase": "Natural Gas",
        "nicht biogenere Abfälle": "Waste",
        "nicht biogener Abfall": "Waste",
        "Wärme": "Other",
        "Biomasse": "Bioenergy",
        "Wasserstoff": "Hydrogen",
    }
    technology = {
        "Verbrennungsmotor": "",
        "Gasturbinen mit Abhitzekessel": "CCGT",
        "Brennstoffzelle": "Fuel Cell",
        "Strilingmotor": "",
        "Stirlingmotor": "",
        "Kondensationsmaschine mit Entnahme": "Steam Turbine",
        "Sonstige": "",
        "Gasturbinen ohne Abhitzekessel": "OCGT",
        "Dampfmotor": "Steam Turbine",
        "Gegendruckmaschine mit Entnahme": "Steam Turbine",
        "Gegendruckmaschine ohne Entnahme": "Steam Turbine",
        "Gasturbinen mit nachgeschalteter Dampfturbine": "CCGT",
        "ORC (Organic Rankine Cycle)-Anlage": "Steam Turbine",
        "Kondensationsmaschine ohne Entnahme": "Steam Turbine",
    }

    CHP_sel.replace({"Fueltype": fueltype, "Technology": technology}, inplace=True)

    def lookup_geodata(missing_plz):
        for i in range(10):
            plz = missing_plz[:-1] + str(i)
            if plz in geodata.index:
                return geodata.loc[plz]
        for i in range(100):
            prefix = "0" if i < 10 else ""
            plz = missing_plz[:-2] + prefix + str(i)
            if plz in geodata.index:
                return geodata.loc[plz]

        return pd.Series((pd.NA, pd.NA))

    missing_i = CHP_sel.lat.isna() | CHP_sel.lon.isna()
    CHP_sel.loc[missing_i, ["lat", "lon"]] = CHP_sel.loc[
        missing_i, "Postleitzahl"
    ].apply(lookup_geodata)

    cols = [
        "Name",
        "Fueltype",
        "Technology",
        "Set",
        "Country",
        "Capacity",
        "Efficiency",
        "DateIn",
        "DateOut",
        "lat",
        "lon",
        "Capacity_thermal",
        "Industry",
    ]

    # convert unit of capacities from kW to MW
    CHP_sel.loc[:, ["Capacity", "Capacity_thermal"]] /= 1e3

    # add missing Fueltype for plants > 100 MW
    fuelmap = {
        "GuD Mitte": "Natural Gas",
        "HKW Mitte": "Natural Gas",
        "GuD Süd": "Natural Gas",
        "HKW Lichterfelde": "Natural Gas",
        "GuD Niehl 2 RheinEnergie": "Natural Gas",
        "HKW Marzahn": "Natural Gas",
        "Gasturbinen Heizkraftwerk Nossener Brücke": "Natural Gas",
        "SEE916495905242": "Natural Gas",
        "HKW Leipzig Nord": "Natural Gas",
        "HKW Reuter": "Waste",
        "Solvay Rb Kraftwerk": "Lignite",
        "GuD Süd Wolfsburg": "Natural Gas",
        "GuD Erfurt Ost": "Natural Gas",
        "KW Nord": "Natural Gas",
        "SEE904887370686": "Oil",
        "GuD2": "Natural Gas",
        "Heizkrafwerk Hafen der Stadtwerke Münster GmbH": "Waste",
        "Kraftwerk Ha": "Natural Gas",
        "Kraftwerk HA": "Natural Gas",
        "PKV Dampfsammelschienen-KWK-Anlage": "Natural Gas",
    }
    CHP_sel["Fueltype"] = (
        CHP_sel["Name"].map(fuelmap).combine_first(CHP_sel["Fueltype"])
    )

    return CHP_sel[cols].copy()


def calculate_efficiency(CHP_de):
    """
    Calculate the efficiency of the CHP plants depending on Capacity and
    DateIn.

    Following Triebs et al. (
    https://doi.org/10.1016/j.ecmx.2020.100068)
    """

    def EXT(cap, year):
        # returns the efficiency for extraction condensing turbine
        return ((44 / 2400) * cap + 0.125 * year - 204.75) / 100

    def BP(cap, year):
        # returns the efficiency for back pressure turbine
        return ((5e-3) * cap + 0.325 * year - 611.75) / 100

    # TODO: differentiate between extraction condensing turbine and back pressure turbine
    CHP_de["Efficiency"] = CHP_de.apply(
        lambda row: BP(row["Capacity"], row["DateIn"]), axis=1
    )

    return CHP_de


if __name__ == "__main__":
    if "snakemake" not in globals():
        snakemake = mock_snakemake(
            "build_existing_chp_de",
            clusters=27,
            run="KN2045_Bal_v4",
        )

    configure_logging(snakemake)
    logger.info("Retrieving and cleaning CHP data from BNetzA")
    biomass = pd.read_csv(
        snakemake.input.mastr_biomass, dtype={"Postleitzahl": str}, low_memory=False
    )
    combustion = pd.read_csv(
        snakemake.input.mastr_combustion, dtype={"Postleitzahl": str}, low_memory=False
    )

    geodata = pd.read_csv(
        snakemake.input.plz_mapping,
        index_col="plz",
        dtype={"plz": str},
        names=["plz", "lat", "lon"],
        skiprows=1,
    )

    logger.info("Cleaning data")
    CHP_de = clean_data(combustion, biomass, geodata)

    logger.info(
        "Calculating efficiency of CHP plants depending on capacity and build year."
    )
    CHP_de = calculate_efficiency(CHP_de)

    logger.info("Mapping CHP plants to regions")
    regions = gpd.read_file(snakemake.input.regions).set_index("name")
    geometry = gpd.points_from_xy(CHP_de["lon"], CHP_de["lat"])
    gdf = gpd.GeoDataFrame(geometry=geometry, crs=4326)
    CHP_de["bus"] = gpd.sjoin_nearest(gdf, regions, how="left")["name"]

    CHP_de.to_csv(snakemake.output.german_chp, index=False)
