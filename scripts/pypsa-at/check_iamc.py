import pandas as pd
from export_iamc_variables import BC_ALIAS


def check_balances(df: pd.DataFrame):
    """

    Parameters
    ----------
    df

    Returns
    -------
    :
    """
    for bc in sorted(set(BC_ALIAS.values())):
        for (year, region), ds in df.groupby(["Year", "Region"]):
            ds = ds.drop("billion EUR2020", level="Unit").droplevel(
                ["Year", "Region", "Scenario", "Model"]
            )

            # ds = filter_by(df, Year="2040", Region="IT0").drop("billion EUR2020", level="Unit").droplevel(
            #     ["Year", "Region", "Scenario", "Model"])

            data = ds.filter(like=bc) / 1e6
            primary = data.filter(regex=rf"^\('Primary Energy\|{bc}\|")
            _primary_agg = data.filter(regex=rf"^\('Primary Energy\|{bc}'")
            assert abs(primary.sum() - _primary_agg.sum()) < 1e-5, (
                f"Sum category mismatch: {year} {region} {bc}"
            )

            secondary_supply = data.filter(regex=rf"^\('Secondary Energy\|{bc}")
            secondary_demand = data.filter(
                regex=rf"^\('Secondary Energy\|[a-zA-Z0-9\s]*\|{bc}"
            )
            final_demand = data.filter(regex=rf"^\('Final Energy\|{bc}\|")

            # assert bypass + secondary supply == Final

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

            if abs(balance) >= 0.5:
                print(bc, year, region, ":\t", f"{balance:.4f}")


if __name__ == "__main__":
    xls = pd.read_excel(
        "/IdeaProjects/pypsa-at/results/v2025.02/KN2045_Mix/evaluation/exported_iamc_variables.xlsx",
        index_col=[0, 1, 2, 3, 4],
    )
    xls.columns.name = "Year"
    check_balances(xls.stack())
