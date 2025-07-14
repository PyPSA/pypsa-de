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
    # for bc in ("AC",):  # sorted(set(BC_ALIAS.values())):
    for bc in sorted(set(BC_ALIAS.values())):
        for (year, region), ds in df.groupby(["Year", "Region"]):
            ds = ds.drop("billion EUR2020", level="Unit").droplevel(
                ["Year", "Region", "Scenario", "Model"]
            )

            # ds = filter_by(df, year="2020", region="PL").drop("billion EUR2020", level="unit").droplevel(
            #     ["year", "region", "scenario", "model"])

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

            balance = (
                primary.sum()
                + secondary_supply.sum()
                - secondary_demand.sum()
                - final_demand.sum()
            )

            # meth = ds.filter(like="Methanolisation")

            if abs(balance) > 0.1:
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
    assert "" not in df.index.unique("unit"), (
        f"Need to assign all energy units with 'MWh_suff'.\n{df[df.index.get_level_values('unit') == '']}"
    )
    # all bus_carrier in BCS ?


if __name__ == "__main__":
    xls = pd.read_excel(
        "/IdeaProjects/pypsa-at/results/v2025.02/KN2045_Mix/evaluation/exported_iamc_variables.xlsx",
        index_col=[0, 1, 2, 3, 4],
    )
    xls.columns.name = "Year"
    check_balances(xls.stack())
