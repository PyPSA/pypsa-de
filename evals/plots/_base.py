# -*- coding: utf-8 -*-
"""Common graph bases and emtpy figures."""

import pathlib
import typing

import pandas as pd
from jinja2 import Template
from plotly import express as px
from plotly import graph_objects as go
from plotly.offline.offline import get_plotlyjs

from evals.configs import PlotConfig
from evals.constants import ALIAS_LOCATION_REV, RUN_META_DATA


class ESMChart:
    """A base class for Energy System Modeling graphs using Plotly.

    Parameters
    ----------
    df
        The data frame with the plot data. The class expects a data
        frame that complies with the metric data model, i.e. has
        the expected column and index labels.

    cfg
        The plotly configuration.
    """

    # todo: avoid inheritance and use MixIn instead

    def __init__(self, df: pd.DataFrame, cfg: PlotConfig) -> None:
        self._df = df
        self.cfg = cfg
        self.fig = go.Figure()
        self.unit = self.cfg.unit or df.attrs["unit"]
        self.metric_name = df.attrs["name"]
        self.location = ""
        self.col_values = ""

    @property
    def empty_input(self) -> bool:
        """Determine if the input DataFrame is empty or all NaN.

        Returns
        -------
        :
            True if the input DataFrame is empty or contains only NaN
            values, False otherwise.
        """
        return self._df.empty or self._df.isna().all().all()

    def to_html(
        self, output_path: pathlib.Path, groupby: list[str], idx: typing.Hashable
    ) -> pathlib.Path:
        """Serialize the Plotly figure to an HTML file.

        Parameters
        ----------
        output_path
            The folder to save the HTML file under.
        groupby
            List of groupby keys needed to fill the file name template.
        idx
            The data frame index from the gropuby clause needed to fill
            the file name template.

        Returns
        -------
        :
            The path of the file written.
        """
        file_name = f"{self.construct_file_name(groupby, idx)}.html"
        file_path = output_path / "HTML" / file_name
        template_html = """\
<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8" />
<meta name="viewport" content="width=device-width, initial-scale=1.0" />
<meta name="esmtools" content="{{ esmtools }}" />
</head>
<body>
    {{ fig }}
</body>
</html>"""

        div = self.fig.to_html(include_plotlyjs="directory", full_html=False)
        with file_path.open("w", encoding="utf-8") as fh:
            fh.write(Template(template_html).render(fig=div, **RUN_META_DATA))

        # need to write the plotly.js too, because to_html does not
        bundle_path = file_path.parent / "plotly.min.js"
        if not bundle_path.exists():
            bundle_path.write_text(get_plotlyjs(), encoding="utf-8")

        return file_path

    def to_json(
        self, output_path: pathlib.Path, groupby: list[str], idx: typing.Hashable
    ) -> pathlib.Path:
        """Serialize the Plotly figure to a JSON file.

        Parameters
        ----------
        output_path
            The folder to save the JSON file under.
        groupby
            List of groupby keys needed to fill the file name template.
        idx
            The data frame index from the gropuby clause needed to fill
            the file name template.

        Returns
        -------
        :
            The path of the file written.
        """
        file_name = f"{self.construct_file_name(groupby, idx)}.json"
        file_path = output_path / "JSON" / file_name
        self.fig.write_json(file_path, engine="auto")

        return file_path

    def construct_file_name(self, groupby: list[str], idx: typing.Hashable) -> str:
        """Construct the file name based on the provided template.

        Parameters
        ----------
        groupby
            List of groupby values.
        idx
            The index used for constructing the file name.

        Returns
        -------
        :
            The constructed file name based on the template and
            provided values.
        """
        idx = [idx] if isinstance(idx, str) else idx
        parts = {"metric": self.metric_name} | {
            g: ALIAS_LOCATION_REV.get(i, i) for g, i in zip(groupby, idx, strict=True)
        }
        return self.cfg.file_name_template.format(**parts)

    @staticmethod
    def custom_sort(
        df: pd.DataFrame, by: str, values: tuple, ascending: bool = False
    ) -> pd.DataFrame:
        """Sort a data frame by first appearance in values.

        Sort a data frame by the given column and first appearance
        in a given iterable.

        Parameters
        ----------
        df
            The dataframe to sort.
        by
            The column name to find values in.
        values
            The values to sort by. The order in this collection defines
            the sort result.
        ascending
            Whether, or not to reverse the result (Plotly inserts
            legend items from top down).

        Returns
        -------
        :
            The sorted data frame.
        """
        if not values:
            return df

        def _custom_order(ser: pd.Series) -> pd.Series:
            """Sort by first appearance in an iterable.

            First, construct a dictionary from the input values with the
            series value as key and the position as value.

            Second, use the dictionary in an anonymous function to get
            the position, or 1000 (to put it last) if a value is not
            found in the data series.

            Parameters
            ----------
            ser
                The pandas Series that should become sorted.

            Returns
            -------
            The sorted pandas Series.
            """
            order = {s: i for i, s in enumerate(values)}
            return ser.apply(lambda x: order.get(x, 1000))

        return df.sort_values(by=by, key=_custom_order, ascending=ascending)

    def _set_base_layout(self) -> None:
        """Set various figure properties."""
        self.fig.update_layout(
            height=800,
            font_family="Calibri",
            plot_bgcolor="#ffffff",
            legend_title_text=self.cfg.legend_header,
        )
        # update axes
        self.fig.update_yaxes(
            showgrid=self.cfg.yaxes_showgrid, visible=self.cfg.yaxes_visible
        )
        self.fig.update_xaxes(
            showgrid=False,
            tickprefix="<b>",
            ticksuffix="</b>",
            tickfont_size=20,
            title_font={"size": 20},
        )
        self.fig.update_layout(
            xaxis={"categoryorder": "category ascending"},
            hovermode="x",  # all categories are shown by mouse-over
        )
        # trace order always needs to be reversed to show correct order
        # of legend entries for relative bar charts
        self.fig.update_layout(legend={"traceorder": "reversed"})

        # export the metadata directly in the Layout property for JSON
        self.fig.update_layout(meta=[RUN_META_DATA])

    def _append_footnotes(self) -> None:
        """Append the footnote(s) at the bottom of the figure."""
        self._append_footnote(self.cfg.footnotes[0], align="left")
        self._append_footnote(self.cfg.footnotes[1], y=-0.2)
        if lines := self._count_footnote_lines():
            self.fig.update_layout(margin={"b": 125 + 50 * lines})

    def _append_footnote(
        self, footnote_text: str, y: float = -0.17, align: str = None
    ) -> None:
        """Append a footnote at the bottom of the Figure.

        Parameters
        ----------
        footnote_text
            The text displayed at the bottom of figures.
        y
            The vertical position of the footnote. Negative values
            move the footnote down.
        align
            The text alignment mode.
        """
        if footnote_text:
            self.fig.add_annotation(
                text=footnote_text,
                xref="paper",
                yref="paper",
                xanchor="left",
                yanchor="top",
                x=0,
                y=y,
                showarrow=False,
                font={"size": 15},
                align=align,
            )

    def _style_title_and_legend_and_xaxis_label(self) -> None:
        """Update figure title and legend."""
        self.fig.update_layout(
            title_font_size=self.cfg.title_font_size,
            font_size=self.cfg.font_size,
            legend={
                "x": 1,
                "y": 1,
                "font": {"size": self.cfg.legend_font_size},
            },
        )
        if self.cfg.xaxis_title:  # allow skipping via empty string
            self.fig.update_layout(xaxis_title=self.cfg.xaxis_title)

    def _count_footnote_lines(self) -> int:
        """Count the number of lines in footnote texts.

        Returns
        -------
        :
            The number of text lines required to write the
            footnote text.
        """
        return "".join(self.cfg.footnotes).count("<br>")


def empty_figure(title: str) -> go.Figure:
    """Return an empty graph with explanation text.

    Parameters
    ----------
    title
        The figure title displayed at the top of the graph.

    Returns
    -------
    :
        The plotly figure with a text that explains that there is no
        data available for this view.
    """
    fig = px.bar(pd.DataFrame(), title=title)
    fig.add_annotation(
        text="No Values to be displayed",
        xref="paper",
        yref="paper",
        xanchor="center",
        yanchor="middle",
        x=0.5,
        y=0.5,
        showarrow=False,
        font={"size": 20},
    )
    fig.update_xaxes(showgrid=False, showticklabels=False)
    fig.update_yaxes(showgrid=False, showticklabels=False)
    fig.update_layout(xaxis_title="", yaxis_title="", plot_bgcolor="white")
    fig.update_layout(meta=[RUN_META_DATA])

    return fig
