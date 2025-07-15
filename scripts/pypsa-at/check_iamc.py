import pandas as pd
from export_iamc_variables import BC_ALIAS
from pyam import IamDataFrame


def check_balances(df: pd.DataFrame):
    """

    Parameters
    ----------
    df

    Returns
    -------
    :
    """
    # for bc in ("AC",):  # sorted(set(BC_ALIAS.values())):
    for bc in sorted(set(BC_ALIAS.values())):
        for (year, region), ds in df.groupby(["Year", "Region"]):
            ds = ds.drop("billion EUR2020", level="Unit").droplevel(
                ["Year", "Region", "Scenario", "Model"]
            )

            # ds = filter_by(df, Year="2040", Region="IT0").drop("billion EUR2020", level="Unit").droplevel(
            #     ["Year", "Region", "Scenario", "Model"])

            data = ds.filter(like=bc) / 1e6
            primary = data.filter(regex=rf"^\('Primary Energy\|{bc}")
            secondary_supply = data.filter(regex=rf"^\('Secondary Energy\|{bc}")
            secondary_demand = data.filter(
                regex=rf"^\('Secondary Energy\|[a-zA-Z0-9\s]*\|{bc}"
            )
            final_demand = data.filter(regex=rf"^\('Final Energy\|{bc}")

            # ambient heat must be excluded because it is not a demand
            secondary_demand = secondary_demand.drop(
                secondary_demand.filter(like="Ambient Heat").index
            )
            # except for Heat to satisfy Heat Loads
            if bc == "Heat":
                ambient_heat = data.filter(like="Ambient Heat")
                secondary_supply = secondary_supply.add(ambient_heat, fill_value=0)

            if bc == "Oil":
                # # need to drop refining losses, because `Oil Primary` is
                # # simplified to `Oil` and refining losses are for the
                # # conversion from oil primary -> oil
                # primary = primary.drop("Primary Energy|Oil|Refining Losses", level="Variable")

                # also drop unsustainable liquids oil demands. They are not
                # correctly found by the regex
                # todo: rename to Oil|Liquids and move to primary
                secondary_demand = secondary_demand.drop(
                    "Secondary Energy|Oil|Oil|Unsustainable Bioliquids",
                    level="Variable",
                )

            balance = (
                primary.sum()
                + secondary_supply.sum()
                - secondary_demand.sum()
                - final_demand.sum()
            )

            if abs(balance) >= 0.1:
                print(bc, year, region, ":\t", f"{balance:.4f}")

    # per carrier:
    #   primary energy
    #   + import
    #   + secondary production
    #   - secondary demand
    #   - secondary losses
    #   - final energy
    #   - export
    #   == 0 ?

    # all units in {units} ?
    # assert "" not in df.index.unique("unit"), (
    #     f"Need to assign all energy units with 'MWh_suff'.\n{df[df.index.get_level_values('unit') == '']}"
    # )
    # all bus_carrier in BCS ?


def connect_sankey(df):
    df = IamDataFrame(df)

    sankey_mapping = {
        "Primary Energy|Coal": ("Coal Mining", "Coal Trade & Power Generation"),
        "Primary Energy|Gas": (
            "Natural Gas Extraction",
            "Gas Network & Power Generation",
        ),
        "Secondary Energy|Electricity|Non-Biomass Renewables": (
            "Non-Biomass Renewables",
            "Electricity Grid",
        ),
        "Secondary Energy|Electricity|Nuclear": ("Nuclear", "Electricity Grid"),
        "Secondary Energy|Electricity|Coal": (
            "Coal Trade & Power Generation",
            "Electricity Grid",
        ),
        "Secondary Energy|Electricity|Gas": (
            "Gas Network & Power Generation",
            "Electricity Grid",
        ),
        "Final Energy|Electricity": ("Electricity Grid", "Electricity Demand"),
        "Final Energy|Solids|Coal": (
            "Coal Trade & Power Generation",
            "Non-Electricity Coal Demand",
        ),
        "Final Energy|Gases": ("Gas Network & Power Generation", "Gas Demand"),
    }

    fig = df.filter(year=2050).plot.sankey(mapping=sankey_mapping)
    fig.show()


if __name__ == "__main__":
    xls = pd.read_excel(
        "/IdeaProjects/pypsa-at/results/v2025.02/KN2045_Mix/evaluation/exported_iamc_variables.xlsx",
        index_col=[0, 1, 2, 3, 4],
    )
    xls.columns.name = "Year"
    check_balances(xls.stack())
