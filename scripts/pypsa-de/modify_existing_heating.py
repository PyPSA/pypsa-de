# -*- coding: utf-8 -*-


import pandas as pd

leitmodell = "REMod v1.0"

existing_heating = pd.read_csv(snakemake.input.existing_heating, index_col=0)

ariadne = pd.read_csv(
    snakemake.input.ariadne,
    index_col=["model", "scenario", "region", "variable", "unit"],
).loc[
    leitmodell,
    snakemake.params.fallback_reference_scenario,
    "Deutschland",
    :,
    "million",
]

print("Heating demand before modification:", existing_heating.loc["Germany"], sep="\n")

mapping = {
    "gas boiler": "Gas Boiler",
    "oil boiler": "Oil Boiler",
    "air heat pump": "Heat Pump|Electrical|Air",
    "ground heat pump": "Heat Pump|Electrical|Ground",
    "biomass boiler": "Biomass Boiler",
}

year = "2020"
for tech in mapping:
    stock = ariadne.at[
        f"Stock|Space Heating|{mapping[tech]}",
        year,
    ]

    peak = (
        stock
        * existing_heating.loc["Germany"].sum()
        / ariadne.at[f"Stock|Space Heating", year]
    )
    existing_heating.at["Germany", tech] = peak


print(
    f"Heating demand after modification with {leitmodell}:",
    existing_heating.loc["Germany"],
    sep="\n",
)


existing_heating.to_csv(snakemake.output.existing_heating)
