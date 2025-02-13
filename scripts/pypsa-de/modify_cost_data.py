# -*- coding: utf-8 -*-

import logging
import os
import re

import numpy as np
import pandas as pd


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
        import sys

        path = "../submodules/pypsa-eur/scripts"
        sys.path.insert(0, os.path.abspath(path))
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "modify_cost_data",
            planning_horizons="2020",
            file_path="../data/costs/",
            file_name="costs_2020.csv",
            cost_horizon="mean",
            run="KN2045_Bal_v4",
        )
    logger = logging.getLogger(__name__)

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
        logger.warning(f"Mean cost scenario for {matched_year}.")
        new_year = matched_year
    elif snakemake.params.cost_horizon == "pessimist":
        logger.warning(f"Pessimistic cost scenario for {matched_year}.")
        new_year = matched_year + 5
    elif snakemake.params.cost_horizon == "optimist":
        logger.warning(f"Optimistic cost scenario for {matched_year}.")
        new_year = matched_year - 5
    else:
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
    print(costs.loc[modifications.index])

    # add carbon component to fossil fuel costs
    investment_year = int(snakemake.wildcards.planning_horizons[-4:])
    if investment_year in snakemake.params.co2_price_add_on_fossils.keys():
        co2_price = snakemake.params.co2_price_add_on_fossils[investment_year]
        logger.warning(
            f"Adding carbon component according to a co2 price of {co2_price} €/t to fossil fuel costs."
        )
        costs = carbon_component_fossils(costs, co2_price)

    logger.warning("Scaling onwind costs towards Fh-ISE  for Germany.")
    # https://github.com/PyPSA/pypsa-ariadne/issues/179
    # https://www.ise.fraunhofer.de/de/veroeffentlichungen/studien/studie-stromgestehungskosten-erneuerbare-energien.html
    costs.at[("onwind", "investment"), "value"] *= 1.12
    print(costs.loc["onwind", "investment"])

    logger.warning("Adding transport costs of 8.8 EUR/MWh to pelletizing costs.")
    # Assumption based on doi:10.1016/j.rser.2019.109506
    costs.at[("biomass boiler", "pelletizing cost"), "value"] += 8.8
    print(costs.loc["biomass boiler", "pelletizing cost"])

    costs.to_csv(snakemake.output[0])
