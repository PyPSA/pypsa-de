# SPDX-FileCopyrightText: 2023-2025 Austrian Gas Grid Management AG
#
# SPDX-License-Identifier: MIT
# For license information, see the LICENSE.txt file in the project root.
"""ESM grouped barcharts."""

from functools import cached_property
from itertools import product

import numpy as np
import pandas as pd
from plotly import express as px
from plotly import graph_objects as go
from plotly.subplots import make_subplots

from evals.constants import DataModel
from evals.plots._base import ESMChart, empty_figure
from evals.utils import apply_cutoff, prettify_number


class ESMGroupedBarChart(ESMChart):
    """
    A class that produces multiple bar charts in subplots.

    Parameters
    ----------
    *args
        Positional arguments of the base class.

    **kwargs
        Key word arguments of the base class.
    """

    def __init__(self, *args: tuple, **kwargs: dict) -> None:
        super().__init__(*args, **kwargs)
        self.location = self._df.index.unique(DataModel.LOCATION)[0]
        self.col_values = self._df.columns[0]

        # self.df is accessed below. location and col_values must be
        # set before the first access to the df property.
        ncols = len(self.df[DataModel.BUS_CARRIER].unique())
        column_widths = [0.85 / ncols] * ncols
        self.fig = make_subplots(
            rows=1, cols=ncols, shared_yaxes=True, column_widths=column_widths
        )

    @cached_property
    def df(self) -> pd.DataFrame:
        """
        Plot data formatted for grouped bar charts.

        Returns
        -------
        :
            The formatted data for creating bar charts.
        """
        df = apply_cutoff(self._df, limit=self.cfg.cutoff, drop=False)
        df = df.reset_index()

        # need to add missing carrier for every sector to prevent
        # broken sort order. If a carrier is missing in, lets say the
        # first subplot, the sort order will be different (compared)
        # to a sector that has all carriers. To prevent this, we add
        # missing dummy carriers before sorting.
        fill_values = product(
            df[DataModel.YEAR].unique(),
            df[DataModel.LOCATION].unique(),
            df[DataModel.CARRIER].unique(),
            df[DataModel.BUS_CARRIER].unique(),
        )
        df_fill = pd.DataFrame(columns=DataModel.YEAR_IDX_NAMES, data=fill_values)
        df_fill[self.col_values] = np.nan
        df = pd.concat([df, df_fill], ignore_index=True)

        # sort every sector in alphabetical order and separately to
        # correctly align sectors and stacked carrier traces in bars
        df_list = []
        for _, df_sector in df.groupby(self.cfg.facet_column, sort=True):
            sorted_sector = self.custom_sort(
                df_sector,
                by=self.cfg.plot_category,
                values=self.cfg.category_orders,
                ascending=True,
            )
            df_list.append(sorted_sector)
        df = pd.concat(df_list)
        # df = df.groupby(self.cfg.facet_column, sort=True).apply(
        #     self.custom_sort,
        #     by=self.cfg.plot_category,
        #     values=self.cfg.category_orders,
        #     ascending=True,
        #     # include_groups=False,
        # )
        # remove NaN categories again after sorting with all categories
        df = df.dropna(how="all", subset=self.col_values)
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
            x=self.cfg.plot_xaxis,
            y=self.col_values,
            facet_col=self.cfg.facet_column,
            facet_col_spacing=0.04,
            pattern_shape=self.cfg.plot_category,
            pattern_shape_map=pattern,
            color=self.cfg.plot_category,
            color_discrete_map=self.cfg.colors,
            text=self.cfg.facet_column,  # needed to rename xaxis and dropped afterward
            title=title,
            custom_data=[self.cfg.plot_category, "display_value"],
        )

        self.fig.for_each_xaxis(self._rename_xaxis)
        self.fig.for_each_xaxis(self._add_total_sum_subplot_traces)
        self.fig.update_annotations(text="")  # remove text labels

        self._set_base_layout()
        self._style_grouped_bars()
        self._style_title_and_legend_and_xaxis_label()
        self._append_footnotes()
        self.fig.for_each_xaxis(self._style_inner_xaxis_labels)

    def _rename_xaxis(self, xaxis: go.XAxis) -> None:
        """
        Update the xaxis labels.

        The function iterates over subplot columns and looks for the
        sector name in the figure data where the xaxis index matches
        and updates the xaxis label and removes the upper text.

        Parameters
        ----------
        xaxis
            The subplot xaxis (a dictionary).

        Notes
        -----
        A better way to set the xaxis labels is desirable. However, I
        could not find a way to replace the 'year' string using the
        'label' argument in plotly.express.bar(), because all columns
        are named the same (='year') and the label argument only
        maps old to new names in a dictionary.
        """
        layout = self.fig["layout"]
        idx = xaxis["anchor"].lstrip("y")
        for data in self.fig["data"]:
            if data["xaxis"] == f"x{idx}":
                sector = data["text"][0]
                layout[f"xaxis{idx}"]["title"]["text"] = f"<b>{sector}"
                break

    def _style_inner_xaxis_labels(self, xaxis: go.XAxis) -> None:
        """
        Set the font size for the inner xaxis labels.

        Parameters
        ----------
        xaxis
            The subplot xaxis (a dictionary-like object).
        """
        xaxis.update(
            tickfont_size=self.cfg.xaxis_font_size, categoryorder="category ascending"
        )

    def _style_grouped_bars(self) -> None:
        """Style bar traces for grouped bar charts."""
        self.fig.update_traces(
            selector={"type": "bar"},
            width=0.8,
            textposition="inside",
            insidetextanchor="middle",
            texttemplate="<b>%{customdata[1]}</b>",
            textangle=0,
            insidetextfont={"size": 16},
            hovertemplate="%{customdata[0]}: %{customdata[1]} " + self.unit,
            hoverlabel={"namelength": 0},
        )

    def _add_total_sum_subplot_traces(self, xaxis: go.XAxis) -> None:
        """
        Add traces for total sum labels in every subplot.

        The Xaxis is needed to parse the subplot position dynamically.
        The method adds text annotations with the total amount of
        energy per stacked bar.

        Parameters
        ----------
        xaxis
            The subplot xaxis (a dictionary-like object).
        """
        idx = xaxis["anchor"].lstrip("y")
        sector = xaxis["title"]["text"].lstrip("<b>")
        values = self.df.query(f"{self.cfg.facet_column} == '{sector}'").copy()

        values["pos"] = values[self.col_values].where(values[self.col_values].gt(0))
        values["neg"] = values[self.col_values].where(values[self.col_values].le(0))

        totals = values.groupby(self.cfg.plot_xaxis).sum(numeric_only=True)
        totals["pos_display"] = totals["pos"].apply(prettify_number)
        totals["neg_display"] = totals["neg"].apply(prettify_number)

        if totals["pos"].sum() > 0:
            scatter = go.Scatter(
                x=totals.index,
                y=totals["pos"] + totals["pos"].abs().max() / 100,
                text=totals["pos_display"],
                texttemplate="<b>%{text}</b>",
                mode="text",
                textposition="top center",
                showlegend=False,
                name="Sum",
                textfont={"size": 18},
                hoverinfo="skip",
            )

            self.fig.add_trace(scatter, col=int(idx) if idx else 1, row=1)

        if totals["neg"].sum() < 0:
            scatter = go.Scatter(
                x=totals.index,
                y=totals["neg"] - totals["neg"].abs().max() / 100,
                text=totals["neg_display"],
                texttemplate="<b>%{text}</b>",
                mode="text",
                textposition="bottom center",
                showlegend=False,
                name="Sum",
                textfont={"size": 18},
                hoverinfo="skip",
            )
            self.fig.add_trace(scatter, col=int(idx) if idx else 1, row=1)

        # totals = values.groupby(self.cfg.plot_xaxis).sum(numeric_only=True)
        # totals["display_value"] = totals[self.col_values].apply(prettify_number)
        # y_offset = totals[self.col_values].abs().max() / 100
        #
        # scatter = go.Scatter(
        #     x=totals.index,
        #     y=totals[self.col_values] + y_offset,
        #     text=totals["display_value"],
        #     texttemplate="<b>%{text}</b>",
        #     mode="text",
        #     textposition="top center",
        #     showlegend=False,
        #     name="Sum",
        #     textfont={"size": 18},
        #     hoverinfo="skip",
        # )
        #
        # self.fig.add_trace(scatter, col=int(idx) if idx else 1, row=1)
