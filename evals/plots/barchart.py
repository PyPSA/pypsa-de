# -*- coding: utf-8 -*-
"""ESM bar charts."""

from functools import cached_property

import numpy as np
import pandas as pd
from plotly import express as px
from plotly import graph_objects as go

from evals.constants import DataModel
from evals.plots._base import ESMChart, empty_figure
from evals.utils import apply_cutoff, prettify_number


class ESMBarChart(ESMChart):
    """The ESM Bar Chart exports metrics as plotly HTML file.

    Parameters
    ----------
    *args
        Positional arguments of the base class.

    **kwargs
        Key word arguments of the base class.
    """

    def __init__(self, *args: object, **kwargs: object) -> None:
        super().__init__(*args, **kwargs)
        self.fig = go.Figure()

        self.location = self._df.index.unique(DataModel.LOCATION)[0]
        self.col_values = self._df.columns[0]

    @cached_property
    def barmode(self) -> str:
        """Determine the barmode for the bar chart.

        Returns
        -------
        :
            The barmode for the bar chart, either "relative" if there
            are both negative and positive values, or "stack" if not.
        """
        has_negatives = self._df.lt(0).to_numpy().any()
        has_positives = self._df.ge(0).to_numpy().any()
        return "relative" if has_negatives and has_positives else "stack"

    @cached_property
    def df(self) -> pd.DataFrame:
        """Plot data formatted for bar charts.

        Returns
        -------
        :
            The formatted data for creating bar charts.
        """
        df = apply_cutoff(self._df, limit=self.cfg.cutoff, drop=self.cfg.cutoff_drop)

        df = self.custom_sort(
            df.reset_index(),
            by=self.cfg.plot_category,
            values=self.cfg.category_orders,
            ascending=True,
        )
        df["display_value"] = df[self.col_values].apply(prettify_number)

        return df

    def plot(self) -> None:
        """Create the bar chart."""
        title = self.cfg.title.format(location=self.location, unit=self.unit)

        if self.empty_input or self.df[self.col_values].isna().all():
            self.fig = empty_figure(title)
            return

        pattern = {
            col: self.cfg.pattern.get(col, "")
            for col in self.df[self.cfg.plot_category].unique()
        }
        self.fig = px.bar(
            self.df,
            color_discrete_map=self.cfg.colors,
            pattern_shape=self.cfg.plot_category,
            pattern_shape_map=pattern,
            barmode=self.barmode,
            x=self.cfg.plot_xaxis,
            y=self.col_values,
            color=self.cfg.plot_category,
            text=self.col_values,
            title=title,
            labels={
                self.col_values: "<b>" + self.unit + "</b>",
                self.cfg.plot_category: self.cfg.legend_header,
            },
            custom_data=[self.cfg.plot_category, "display_value"],
        )

        self._set_base_layout()
        self._style_bars()
        self._style_title_and_legend_and_xaxis_label()
        self._append_footnotes()
        self.fig.for_each_trace(self._set_legend_rank, selector={"type": "bar"})

        # add total sum labels at the end of the bar(s)
        if self.barmode == "relative":
            self.fig.add_hline(y=0)  # visual separator between supply and withdrawal
            self._add_total_sum_trace("Lower Sum", orientation="down")
            self._add_total_sum_trace("Upper Sum", orientation="up")
        else:
            self._add_total_sum_trace("Sum")

    def _add_total_sum_trace(
        self,
        name_trace: str,
        orientation: str = None,
    ) -> None:
        """Create a scatter trace for total sum labels.

        The label will show the total sum as text at the end of the bar
        trace. The label will not be part of the legend.

        Parameters
        ----------
        name_trace
            The name of the trace useful to identify the trace in the
            JSON representation.
        orientation : optional, {'up', 'down', None}
            The orientation of the trace used to choose the
            text position and the sign of the values.
        """
        sign = 1
        if orientation == "up":
            values = self.df[self.df[self.col_values].gt(0)]
        elif orientation == "down":
            sign = -1
            values = self.df[self.df[self.col_values].le(0)]
        else:  # barmode = stacked
            values = self.df

        totals = values.groupby(self.cfg.plot_xaxis).sum(numeric_only=True)
        totals["display_value"] = totals[self.col_values].apply(prettify_number)
        y_offset = totals[self.col_values].abs().max() / 100 * sign

        scatter = go.Scatter(
            x=totals.index,
            y=totals[self.col_values] + y_offset,
            text=totals["display_value"],
            texttemplate="<b> %{text} " + self.unit + "</b>",
            mode="text",
            textposition=f"{'bottom' if orientation == 'down' else 'top'} center",
            showlegend=False,
            name=name_trace,
            textfont={"size": 18},
            hoverinfo="skip",
        )

        self.fig.add_trace(scatter)

    def _style_bars(self) -> None:
        """Update bar trace styles."""
        self.fig.update_traces(
            selector={"type": "bar"},
            width=0.6,
            textposition="inside",
            insidetextanchor="middle",
            texttemplate="<b>%{customdata[1]}</b>",
            insidetextfont={"size": 16},
            textangle=0,
            hovertemplate="%{customdata[0]}: %{customdata[1]} " + self.unit,
            hoverlabel={"namelength": 0},
        )

    def _set_legend_rank(self, trace: go.Bar) -> go.Bar:
        """Set the legendrank attribute for bar traces.

        Only traces listed in the category_order
        configuration item are considered.

        Parameters
        ----------
        trace
            The trace object to set the legendrank for.

        Returns
        -------
        :
            The updated bar trace, or the original bar trace.
        """
        if trace["name"] in self.cfg.category_orders:
            y = trace["y"]  # need to drop nan and inf or the sum may return NaN
            trace_sum = y[np.isfinite(y)].sum()
            sign = -1 if trace_sum < 0 else 1
            pos = self.cfg.category_orders.index(trace["name"])
            trace = trace.update(legendrank=1000 + pos * sign)  # 1000 = plotly default
        return trace
