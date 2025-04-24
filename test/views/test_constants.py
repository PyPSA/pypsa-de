import inspect
import itertools

import pytest
from esmtools import constants
from esmtools.utils import get_mapping

new_collection_names_internal = {
    "tech_capacity_e": "",
    "tech_capacity_f": "",
    "tech_capacity_net": "",
    "tech_demand_e": "",
    "tech_demand_inflexible_e": "",
    "tech_gas_demand": "",
    "tech_gas_demand_eu": "",
    "tech_gas_demand_ts": "",
    "tech_gas_prod": "",
    "tech_gas_prod_eu": "",
    "tech_gas_prod_ts": "",
    "tech_gas_store": "",
    "tech_h2_demand": "",
    "tech_h2_demand_eu": "",
    "tech_h2_demand_ts": "",
    "tech_h2_prod": "",
    "tech_h2_prod_eu": "",
    "tech_h2_prod_ts": "",
    "tech_heat_demand_ts": "",
    "tech_heat_prod_ts": "",
    "tech_production_e": "",
}


# @pytest.mark.parametrize(
#     "component",
#     [
#         "capacities",
#         "fuel_capacities",
#         "gas_store_capacities",
#         "fuel_net_capacities",
#         "elec_demand",
#         "elec_prod",
#         "gas_production_mapping_EU",
#         "gas_production_mapping",
#         "gas_demand_time_series",
#         "gas_demand_mapping_EU",
#         "gas_demand_mapping",
#         "gas_prod_time_series",
#         "h2_production_mapping_EU",
#         "h2_production_mapping",
#         "h2_demand_mapping_EU",
#         "h2_demand_mapping",
#         "co2_emissions_mapping",
#         "h2_demand_time_series",
#         "h2_prod_time_series",
#         "district_heat_prod_time_series",
#         "district_heat_dem_time_series",
#     ],
# )
# def test_mapping_internal(component):
#     """Ensure correct collections."""
#     from Toolbox.get_style_dictionaries import get_filter_dictionaries
#
#     old = get_filter_dictionaries(component, "EN", "int")
#     new = inspect.getmembers(CarrierFilter)
#
#
# @pytest.mark.parametrize(
#     "component",
#     [
#         "capacities",
#         "fuel_capacities",
#         "fuel_net_capacities",
#         "elec_demand",
#         "elec_prod",
#         "elec_times_series_dem",
#         "elec_times_series_prod",
#         "gas_production_mapping",
#         "h2_production_mapping",
#         "co2_emissions_mapping",
#     ],
# )
# def test_mapping_external():
#     """"""


@pytest.mark.unit()
@pytest.mark.parametrize("mapping_type", ["internal", "external"])
def test_check_all_groups_used_in_mapping(mapping_type):
    """"""
    members = inspect.getmembers(constants.Group)
    groups = [m[1] for m in members if not m[0].startswith("_")]
    assert len(set(groups)) == len(groups), "Duplicated group names found."
    mapping = get_mapping("default", map_type="external")
    mapping_groups = set(mapping.values())
    # assert


def all_mappings():
    old_new_components = {
        ("capacities", "default"),
        ("co2_emissions_mapping", "co2"),
        ("district_heat_dem_time_series", "district_heat"),
        ("district_heat_prod_time_series", "district_heat"),
        ("elec_demand", "e_demand"),
        ("elec_prod", "e_production"),
        ("elec_times_series_dem", "electricity"),
        ("elec_times_series_prod", "electricity"),
        ("fuel_capacities", "capacity"),
        ("fuel_net_capacities", "default"),
        ("gas_demand_mapping", "capacity"),
        ("gas_demand_mapping_EU", "capacity"),
        ("gas_demand_time_series", "capacity"),
        ("gas_prod_time_series", "gas_prod"),  # not used anymore
        ("gas_production_mapping", "gas_prod"),  # not used anymore
        ("gas_production_mapping_EU", "default"),
        ("gas_store_capacities", "default"),
        ("h2_demand_mapping", "h2_demand"),
        ("h2_demand_mapping_EU", "h2_demand_eu"),  # not used anymore
        ("h2_demand_time_series", "h2_demand_ts"),  # not used anymore
        ("h2_prod_time_series", "h2_prod"),  # not used anymore
        ("h2_production_mapping", "default"),
        ("h2_production_mapping_EU", "h2_prod"),  # not used anymore
    }

    yield from itertools.product(["external", "internal"], old_new_components)


@pytest.mark.unit
@pytest.mark.parametrize("mapping", all_mappings(), ids=lambda p: f"{p[1][0]}-{p[0]}")
# def test_mapping_migration(external, component, ported):
def test_mapping_migration(mapping):
    """Test if the simplified mapping contains the toolbox mapping."""
    from toolbox import get_filter_dictionaries

    map_type, (component, map_name) = mapping
    skip_components = ("elec_times_series_dem", "elec_times_series_prod")
    if component in skip_components and map_type == "internal":
        pytest.skip(f"No internal mappings exist for '{component}'.")
    # the following mappings are not in use anymore:
    #  - "gas_prod"
    #  - "gas_store"
    old = get_filter_dictionaries(
        component, int_ext=map_type[:3], language_parameter="EN"
    )
    new_full = get_mapping(map_name=map_name, map_type=map_type)
    new_reduced = {k: v for k, v in new_full.items() if k in old}

    # some carriers are not needed anymore in the new evaluation.
    redundant_carrier = (
        "H2 Import NAF",
        "H2 Import RU",
        "H2 from SMR",
        "H2 from SMR CC",
        "H2 domestic import",
        "H2 domestic export",
        "H2 foreign import",
        "H2 foreign export",
        "H2 retro domestic import",
        "H2 retro domestic export",
        "H2 retro foreign import",
        "H2 retro foreign export",
        "H2 underground",
        "gas",
        "gas foreign export",
        "gas domestic export",
    )
    [old.pop(k, None) for k in redundant_carrier]

    # CH4 demand is mapped to 'Electricity OCGT' instead of
    # 'Electricity (OCGT)' in the new mappings
    if old.get("OCGT") == "Electricity (OCGT)":
        old["OCGT"] = "Electricity OCGT"

    assert old == new_reduced


@pytest.mark.unit
def test_all_groups_have_a_color_assigned():
    """Verify that all groups have a color assigned."""


@pytest.mark.unit
def test_all_carrier_have_a_group_assigned():
    """Verify that all carriers have a mapped group name."""


@pytest.mark.unit
def test_all_bus_carrier_have_a_group_assigned():
    """Verify that all bus carriers have a mapped group name."""
