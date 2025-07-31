from itertools import product

import pandas as pd
import pytest

from test.test_evals import Constants


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
    return iamc.stack()


@pytest.mark("integration")
@pytest.mark.parametrize(
    ("year", "region"),
    [
        product(Constants.YEARS, Constants.LOCATIONS),
    ],
)
def test_biogas_balances(
    exported_variables: pd.DataFrame, year, region, threshold: float = 0.1
) -> None:
    """
    Test if Solid Biomass supply and demand are balanced.

    Parameters
    ----------
    exported_variables
        The IAMC variables table.
    year
    region
    threshold
        The allowed imbalance in TWh.

    Returns
    -------
    :
    """
    ds = (
        exported_variables.query("Year == @year and Region == @region")
        .drop("billion EUR2020", level="Unit")
        .droplevel(["Year", "Region", "Scenario", "Model"])
        .div(1e6)  # to TWh for readability
    )

    bc = "Biomass"
    primary = ds.filter(regex=rf"^\('Primary Energy\|{bc}")
    secondary_supply = ds.filter(regex=rf"^\('Secondary Energy\|{bc}")
    secondary_demand = ds.filter(regex=rf"^\('Secondary Energy\|[a-zA-Z0-9\s]*\|{bc}")
    final_demand = ds.filter(regex=rf"^\('Final Energy\|{bc}")

    # ambient heat must be excluded or balances cannot amount to zero
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
