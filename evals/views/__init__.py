# -*- coding: utf-8 -*-
"""Contains the views module."""

from evals.views.capacity.electricity_capacity import view_electricity_capacities
from evals.views.capacity.heat_capacity import view_heat_capacity

__all__ = ["view_heat_capacity", "view_electricity_capacities"]
