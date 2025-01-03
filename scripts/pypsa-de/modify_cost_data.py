# -*- coding: utf-8 -*-

import logging
import os
import re
import sys

import numpy as np
import pandas as pd

from scripts._helpers import configure_logging, mock_snakemake

logger = logging.getLogger(__name__)


def carbon_component_fossils(costs, co2_price):
    """
    Add carbon component to fossil fuel costs.
    """

    carriers = ["gas", "oil", "lignite", "coal"]
    # specific emissions in tons CO2/MWh according to n.links[n.links.carrier =="your_carrier].efficiency2.unique().item()
    specific_emisisons = {
        "oil": 0.2571,
        "gas": 0.198,  # OCGT
        "coal": 0.3361,
        "lignite": 0.4069,
    }

    for c in carriers:
        carbon_add_on = specific_emisisons[c] * co2_price
        costs.at[(c, "fuel"), "value"] += carbon_add_on
        add_str = f" (added carbon component of {round(carbon_add_on,4)} €/MWh according to co2 price of {co2_price} €/t co2 and carbon intensity of {specific_emisisons[c]} t co2/MWh)"
        if pd.isna(costs.at[(c, "fuel"), "further description"]):
            costs.at[(c, "fuel"), "further description"] = add_str
        else:
            costs.at[(c, "fuel"), "further description"] = (
                str(costs.at[(c, "fuel"), "further description"]) + add_str
            )

    return costs


if __name__ == "__main__":
    if "snakemake" not in globals():
        snakemake = mock_snakemake(
            "modify_cost_data",
            planning_horizons="2020",
            file_path="../data/costs/",
            file_name="costs_2020.csv",
            cost_horizon="mean",
            run="KN2045_Bal_v4",
        )
    configure_logging(snakemake)

    # read in cost data from technology-data library
    costs = os.path.join(
        snakemake.params.file_path,
        snakemake.params.cost_horizon,
        snakemake.params.file_name,
    )

    # cost_horizon is a setting for technology-data and specifies either
    # mean, pessimist or optimist cost scenarios
    # the cost modifications file contains specific cost assumptions for
    # germany, developed in the ARIADNE project
    # here pessimist and optimistic scenarios correspond to a delay or a
    # speed up in cost reductions

    costs = pd.read_csv(costs, index_col=[0, 1]).sort_index()

    matched_year = int(
        re.search(
            r"costs_(\d{4})-modifications\.csv", snakemake.input.modifications
        ).group(1)
    )

    if matched_year <= 2020 or snakemake.params.cost_horizon == "mean":
        logger.info(f"Mean cost scenario for {matched_year}.")
        new_year = matched_year
    elif snakemake.params.cost_horizon == "pessimist":
        logger.info(f"Pessimistic cost scenario for {matched_year}.")
        new_year = matched_year + 5
    elif snakemake.params.cost_horizon == "optimist":
        logger.info(f"Optimistic cost scenario for {matched_year}.")
        new_year = matched_year - 5
    else:
        logger.error(
            "Invalid specification of cost options. Please choose 'mean', 'pessimist' or 'optimist' as config[costs][horizon]."
        )
        raise ValueError("Invalid specification of cost options.")

    new_filename = re.sub(
        r"costs_\d{4}-modifications\.csv",
        f"costs_{new_year}-modifications.csv",
        snakemake.input.modifications,
    )
    modifications = pd.read_csv(new_filename, index_col=[0, 1]).sort_index()
    if snakemake.params.NEP == 2021:
        modifications = modifications.query("source != 'NEP2023'")
    elif snakemake.params.NEP == 2023:
        modifications = modifications.query("source != 'NEP2021'")
    else:
        logger.warning(
            f"NEP year {snakemake.params.NEP} is not in modifications file. Falling back to NEP2021."
        )
        modifications = modifications.query("source != 'NEP2023'")

    costs.loc[modifications.index] = modifications
    logger.info(
        f"Modifications to the following technologies are applied:\n{list(costs.loc[modifications.index].index.get_level_values(0))}."
    )

    # add carbon component to fossil fuel costs
    investment_year = int(snakemake.wildcards.planning_horizons[-4:])
    if (snakemake.params.co2_price_add_on_fossils is not None) and (
        investment_year in snakemake.params.co2_price_add_on_fossils.keys()
    ):
        co2_price = snakemake.params.co2_price_add_on_fossils[investment_year]
        logger.info(
            f"Adding carbon component according to a co2 price of {co2_price} €/t to fossil fuel costs."
        )
        costs = carbon_component_fossils(costs, co2_price)

    logger.info(
        f"Scaling onwind costs towards Fh-ISE for Germany: {costs.loc["onwind", "investment"].value} {costs.loc['onwind', 'investment'].unit}."
    )
    # https://github.com/PyPSA/pypsa-ariadne/issues/179
    # https://www.ise.fraunhofer.de/de/veroeffentlichungen/studien/studie-stromgestehungskosten-erneuerbare-energien.html
    costs.at[("onwind", "investment"), "value"] *= 1.12

    # Assumption based on doi:10.1016/j.rser.2019.109506
    costs.at[("biomass boiler", "pelletizing cost"), "value"] += 8.8
    logger.info(
        f"Adding transport costs of 8.8 EUR/MWh to solid biomass pelletizing costs. New value: {costs.loc['biomass boiler', 'pelletizing cost'].value} {costs.loc['biomass boiler', 'pelletizing cost'].unit}."
    )

    # Klimaschutz- und Energieagentur Baden-Württemberg (KEA) Technikkatalog

    costs.at[("central water tank storage", "investment"), "value"] *= (
        1.12 / 0.6133
    )  # KEA costs / 2020 costs
    logger.info(
        f"Scaling central water tank storage investment costs to KEA Technikkatalog: {costs.loc['central water tank storage', 'investment'].value} {costs.loc['central water tank storage', 'investment'].unit}."
    )

    # increase central gas CHP lifetime to 40 years
    costs.at[("central gas CHP", "lifetime"), "value"] = 40
    logger.info(
        f"Setting lifetime of central gas CHP to {costs.at[("central gas CHP" , "lifetime") , "value"]} {costs.at[("central gas CHP" , "lifetime") , "unit"]}."
    )

    costs.to_csv(snakemake.output[0])
