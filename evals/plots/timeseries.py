# SPDX-FileCopyrightText: 2023-2025 Austrian Gas Grid Management AG
#
# SPDX-License-Identifier: MIT
# For license information, see the LICENSE.txt file in the project root.
"""ESM time series scatter plots."""

from functools import cached_property

import pandas as pd
from plotly import graph_objects as go

from evals.constants import DataModel
from evals.plots._base import ESMChart, empty_figure
from evals.utils import apply_cutoff


class ESMTimeSeriesChart(ESMChart):
    """
    A class that produces one time series chart.

    Parameters
    ----------
    *args
        Positional arguments of the base class.

    **kwargs
        Key word arguments of the base class.
    """

    def __init__(self, *args: tuple, **kwargs: dict) -> None:
        super().__init__(*args, **kwargs)
        self.fig = go.Figure()
        self.year = self._df.index.unique("year")[0]
        self.yaxes_showgrid = self.yaxes_visible = True
        self.location = self._df.index.unique(DataModel.LOCATION)[0]

    @cached_property
    def df(self) -> pd.DataFrame:
        """
        Plot data formatted for time series charts.

        Returns
        -------
        :
            The formatted data for creating bar charts.
        """
        df = apply_cutoff(self._df, limit=self.cfg.cutoff, drop=self.cfg.cutoff_drop)
        df = self.custom_sort(
            df, by=self.cfg.plot_category, values=self.cfg.category_orders
        )
        df = self.fix_snapshots(df, int(self.year))
        df = df.droplevel([DataModel.YEAR, DataModel.LOCATION])

        return df.T  # transpose to iterate column wise over categories

    def plot(self) -> None:
        """
        Plot the data to the chart.

        This function iterates over the data series, adds traces to the
        figure, styles the inflexible demand, sets the layout, styles
        the title, legend, x-axis label, time series axes, and appends
        footnotes.
        """
        title = self.cfg.title.format(
            location=self.location, year=self.year, unit=self.unit
        )
        if self.empty_input:
            self.fig = empty_figure(title)
            return

        stackgroup = None
        for i, (name, series) in enumerate(self.df.items()):
            if self.cfg.stacked:
                stackgroup = "supply" if series.sum() >= 0 else "withdrawal"
            legendrank = 1000 + i if stackgroup == "supply" else 1000 - i
            self.fig.add_trace(
                go.Scatter(
                    x=series.index,
                    y=series.values,
                    hovertemplate="%{y:.2f} " + self.unit,
                    name=name,
                    fill=self.cfg.fill.get(name, "tonexty"),
                    fillpattern_shape=self.cfg.pattern.get(name),
                    line_dash=self.cfg.line_dash.get(name, "solid"),
                    line_width=self.cfg.line_width.get(name, 1),
                    line_color=self.cfg.colors.get(name),
                    line_shape=self.cfg.line_shape,
                    fillcolor=self.cfg.colors.get(name),
                    stackgroup=stackgroup,
                    legendrank=legendrank,
                )
            )

        self._style_inflexible_demand()
        self._set_base_layout()
        self._style_title_and_legend_and_xaxis_label()
        self._style_time_series_axes_and_layout(title)
        self._append_footnotes()

    @staticmethod
    def fix_snapshots(df: pd.DataFrame, year: int) -> pd.DataFrame:
        """
        Correct the year in snapshot timestamp column labels.

        Parameters
        ----------
        df
            The DataFrame with timestamps to be adjusted.
        year
            The correct year to use in the data frame columns.

        Returns
        -------
        :
            The DataFrame with corrected timestamps.
        """
        if isinstance(df.columns, pd.DatetimeIndex):
            df.columns = [s.replace(year=year) for s in df.columns]
        return df

    def _style_inflexible_demand(self) -> None:
        """Set the inflexible demand style if it exists."""
        self.fig.update_traces(
            selector={"name": "Inflexible Demand"},
            fillcolor=None,
            fill=None,
            stackgroup=None,
            legendrank=2000,  # first entry in legend (from top)
        )

    def _style_time_series_axes_and_layout(self, title) -> None:
        """
        Update the layout and axes for time series charts.

        Parameters
        ----------
        title
            The figure title to show at the top of the graph.
        """
        self.fig.update_yaxes(
            tickprefix="<b>",
            ticksuffix="</b>",
            tickfont_size=15,
            color=self.cfg.yaxis_color,
            title_font_size=15,
            tickformat=".0f",  # if "TW" in self.unit else ".3f",
            gridwidth=1,
            gridcolor="gainsboro",
        )
        self.fig.update_xaxes(ticklabelmode="period")
        self.fig.update_layout(
            title=title,
            yaxis_title=self.unit,
            hovermode="x",
        )
