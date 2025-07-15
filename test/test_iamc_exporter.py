import pandas as pd
import pytest


@pytest.fixture(scope="module")
def exported_variables(results) -> pd.Series:
    """
    The exported variables table.

    Parameters
    ----------
    results
        The path to the run results folder.

    Returns
    -------
    :
        The iamc data with index level names and in long format.
    """
    iamc = pd.read_excel(
        results / "evaluation" / "exported_iamc_variables.xlsx",
        index_col=[0, 1, 2, 3, 4],
    )
    iamc.columns.name = "Year"
    return iamc.stack()  # long format data series


@pytest.mark("integration")
def test_biogas_balances(df: pd.DataFrame, threshold: float = 0.1) -> None:
    """
    Test if Solid Biomass supply and demand are balanced.

    Parameters
    ----------
    df
        The IAMC variables table.
    threshold
        The allowed imbalance in TWh.

    Returns
    -------
    :
    """
    bc = "Biomass"

    for (year, region), ds in df.groupby(["Year", "Region"]):
        ds = ds.drop("billion EUR2020", level="Unit").droplevel(
            ["Year", "Region", "Scenario", "Model"]
        )
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

        assert abs(balance) <= threshold, f"Imbalances detected: {balance:.3f} TWh"
