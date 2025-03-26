# numpydoc ignore=PR02
"""Module for Energy System Modeling evaluations.

ESM evaluation functions use solved PyPSA networks to produce
human-readable output files and beautiful HTML graphs.
Every evaluation answers common questions in the energy sector.

Every evaluation function is an orthogonal unit of work, i.e. one
evaluation function does not want or require one another.

Any evaluation function exposes the same interface:

Parameters
----------
result_path : pathlib.Path
    The path to the result folder, e.g. 'pypsa-eur-sec/results'.
networks : dict[str, pypsa.Network]
    The loaded networks patched with the ExtendedStatisticsAccessor.
    Dictionary keys must the pypsa.Network 4 digit planning horizon
    (year). The dictionary values are pypsa.Network objects. It is
    recommended to use the esmtools.fileio.read_networks function to
    obtain the dictionary.
subdir : str, default="esm_run/evaluation"
    The output directory relative to the result_path. All results
    will be written to results_path / subdirectory.
    Depending on the output artifact, subfolders will be created,
    e.g. all result files are written to
    results_path / subdirectory / html.

Notes
-----
Evaluation functions may be imported and called just like any other
python function, or collectively via the run_eval script. The
script is only accessible after package installation. See
[Installation](../../how-to-guides/installation.md) section
for instructions.
"""

from esmtools.evals.capacity import (
    eval_electricity_capacities,
    eval_fuel_production_capacities,
)
from esmtools.evals.carbon_emission import eval_co2
from esmtools.evals.curtailment import eval_curtailment
from esmtools.evals.electricity import (
    eval_electricity_balance,
    eval_electricity_balance_ts,
    eval_electricity_demand,
    eval_electricity_production,
    eval_residual_load,
)
from esmtools.evals.final_demand import (
    eval_fed,
    eval_fed_biomass,
    eval_fed_building_heat,
    eval_fed_hh_services,
    eval_fed_sectoral,
)
from esmtools.evals.heat import (
    eval_district_heat_balance,
    eval_district_heat_balance_ts,
)
from esmtools.evals.hydrogen import eval_h2_balance, eval_h2_balance_ts
from esmtools.evals.industry import eval_industry_sectoral, eval_industry_total
from esmtools.evals.methane import eval_ch4_balance, eval_ch4_balance_ts
from esmtools.evals.price import eval_capex, eval_market_value, eval_opex, eval_revenue
from esmtools.evals.primary_energy import eval_heat_primary_energy
from esmtools.evals.storage import (
    eval_electricity_storage,
    eval_gas_storage_capacities,
    eval_phs_hydro_operation,
)
from esmtools.evals.transmission import eval_grid_capactiy
from esmtools.evals.transport import eval_transport_sectoral, eval_transport_total
from esmtools.views.capacity.heat_capacity import view_heat_capacity

__all__ = [
    "eval_capex",
    "eval_ch4_balance",
    "eval_ch4_balance_ts",
    "eval_co2",
    "eval_curtailment",
    "eval_electricity_demand",
    "eval_electricity_production",
    "eval_electricity_balance",
    "eval_electricity_balance_ts",
    "eval_electricity_capacities",
    "eval_electricity_storage",
    "eval_fed",
    "eval_fed_biomass",
    "eval_fed_building_heat",
    "eval_fed_hh_services",
    "eval_fed_sectoral",
    "eval_fuel_production_capacities",
    "eval_residual_load",
    # "eval_frequency_analysis",
    "eval_gas_storage_capacities",
    "eval_h2_balance",
    "eval_h2_balance_ts",
    "eval_heat_primary_energy",
    "eval_district_heat_balance",
    "eval_district_heat_balance_ts",
    "eval_industry_sectoral",
    "eval_industry_total",
    "eval_market_value",
    "eval_opex",
    "eval_phs_hydro_operation",
    # "eval_sankey",
    # "eval_sankey_detail",
    "eval_revenue",
    "eval_grid_capactiy",
    "eval_transport_sectoral",
    "eval_transport_total",
    # new structure:
    "view_heat_capacity",
]
