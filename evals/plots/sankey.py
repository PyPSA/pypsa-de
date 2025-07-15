"""Module for Sankey diagram."""

import numpy as np
import pandas as pd
import plotly
import pyam
from plotly.graph_objs import Figure, Sankey

from evals.constants import DataModel
from evals.utils import filter_by, insert_index_level, rename_aggregate

pd.set_option("display.width", 250)
pd.set_option("display.max_columns", 20)


class ProtoSankey:
    """A Sankey tracer bullet to enhance understanding."""

    def __init__(self) -> None:
        idx = pd.MultiIndex.from_product(
            (["a", "b", "c"], ["AC", "gas"]),
            names=["carrier", "bus_carrier"],
        )
        df_start = pd.DataFrame(index=idx, columns=["value"], data=10)
        idx = pd.MultiIndex.from_product(
            (["a", "d", "e"], ["AC", "gas"]),
            names=["carrier", "bus_carrier"],
        )
        df_end = pd.DataFrame(index=idx, columns=["value"], data=10)

        idx = pd.MultiIndex.from_product(
            (["b", "e"], ["AC", "gas"]),
            names=["carrier", "bus_carrier"],
        )
        df_transform = pd.DataFrame(index=idx, columns=["value"], data=10)

        self.slices = [df_start, df_transform, df_end]
        self.nodes = {"label": []}
        self.links = {"source": [], "target": [], "value": []}

    def plot(self):
        raise NotImplementedError
        # nodes and links are to be maintained separately

        # a node levels are a data frame with the following columns:
        # index: (level, carrier, bus_carrier)
        # columns: [
        #     "value",  # the energy in TWh
        #     "node_id",  # the integer unique id of the node
        #     "link_label",  # the display name of the link
        #     "direction"  # if the "value" is oriented to the left
        #                  # side, or to the right side
        # ]

        # prepare start nodes (generation)
        # prepare first summary nodes
        # prepare transformation & storage nodes left side (input)
        # prepare transformation & storage nodes right (output)
        # prepare passthrough nodes (remainder from transformation)
        # prepare second summary nodes
        # prepare end nodes (final demand)

        # link start nodes to first summary

    def prepend_slice_level_to_multiindex(self):
        """Prepend the slice level to differentiate duplicated nodes."""
        self.slices = [
            insert_index_level(df, str(i), "level") for i, df in enumerate(self.slices)
        ]

    def node_map(self):
        """Return the node map"""
        _df = pd.concat(self.slices)
        _df["node_id"] = [*range(len(_df))]
        return _df.drop("value", axis=1)

    def make_links(self):
        raise NotImplementedError


class OverviewSankey:
    """
    Overview Sankey diagram for one country and year.

    Follows the (MVC) model view component architecture. The data model
    is a collection of data frames with "carrier" and bus_carrier index
    levels. See the input arguments description for the 4 expected data
    frames.
    The view is the transformed input data stored in the "nodes" and
    "links" instance attributes. This data structure is close to the
    component: the plotly.Sankey graph object.

    Parameters
    ----------
    df
    year
    location
    """

    def __init__(self, df: pd.DataFrame, year: str, location: str) -> None:
        self._df = df  # keep a reference to the original data set
        self.location = location
        self.year = year
        self.model = df.droplevel(DataModel.LOCATION).droplevel(DataModel.YEAR)
        self.view = {
            "link": {
                "source": np.array([], dtype=np.int32),
                "target": np.array([], dtype=np.int32),
                "value": np.array([], dtype=np.float32),
                "label": np.array([], dtype=str),
            },
            "node": {
                "label": np.array([], dtype=str),
            },
        }
        self.fig: Figure = Figure()  # component in MVC

        # store references to left and right side link pairs
        self.pairs = []

    def plot(self) -> None:
        """Create the plotly Figure object from the view."""
        self.link_generation_to_first_summary_level()
        self.link_first_summary_level_to_transformation_input()

        sankey_trace = Sankey(
            link=self.view["link"],
            node=self.view["node"],
        )
        self.fig = Figure(sankey_trace)

    def check_input(self) -> None:
        """
        Check the slice inputs.

        This should not be the responsibility of the plotter class.

        Raises
        ------
        AssertionError
            If validity checks fail.
        """
        # locations = pd.Index()
        # years = pd.Index()
        # for df in self.raw_input:
        #     locations.append(df.index.unique(DataModel.LOCATION))
        #     years.append(df.index.unique(DataModel.YEAR))
        # assert len(locations) == 1, "Multiple locations are not supported."
        # assert len(years) == 1, "Multiple years are not supported."

        # todo: assert transformation sums are equal
        # todo: assert generation and demand sums are equal

    def link_generation_to_first_summary_level(self) -> None:
        """Prepare first links."""
        start = self.get_start_node_id()
        left = self.prepare_link_side(
            self.model["Generation"].to_frame("value"),
            nodes_idx_name=DataModel.CARRIER,
            links_idx_name=DataModel.BUS_CARRIER,
            start=start,
        )

        # # we must add left side node ids before preparing the right
        # # side to update the max node id
        # self.extend_view("link", "source", left["node_id"].to_numpy())

        # right side is a summary node:
        right = self.prepare_link_side(
            self.model["Generation"]
            .to_frame("value")
            .groupby(DataModel.BUS_CARRIER)
            .sum(),
            nodes_idx_name=DataModel.BUS_CARRIER,
            links_idx_name=DataModel.BUS_CARRIER,
            start=left["node_id"].max() + 1,
        )

        # many-to-one relation between left and right
        link = left.merge(
            right,
            on=[DataModel.BUS_CARRIER],  # DataModel.CARRIER,
            how="inner",
            suffixes=("_left", "_right"),
        )

        self.pairs.append((left, right))

        self.extend_view("link", "source", link["node_id_left"].to_numpy())
        self.extend_view("link", "target", link["node_id_right"].to_numpy())
        self.extend_view("link", "value", link["value_left"].to_numpy())
        self.extend_view("link", "label", link["link_label_left"].to_numpy())

        left_labels = left.set_index("node_id")["node_label"].drop_duplicates()
        self.extend_view("node", "label", left_labels.to_numpy())

        right_labels = right.set_index("node_id")["node_label"].drop_duplicates()
        self.extend_view("node", "label", right_labels.to_numpy())

    def link_first_summary_level_to_transformation_input(self) -> None:
        """Prepare second links."""
        start = self.get_start_node_id()
        left = self.pairs[0][1]
        # left = self.prepare_link_side(
        #     self.model["Transformation Input"].to_frame("value"),
        #     nodes_idx_name=DataModel.BUS_CARRIER,
        #     links_idx_name=DataModel.BUS_CARRIER,
        #     start=start,
        # )
        # self.extend_view("links", "source", left["node_id"].to_numpy())

        right = self.prepare_link_side(
            self.model["Transformation Input"]
            .to_frame("value")
            .groupby(DataModel.BUS_CARRIER)
            .sum(),
            nodes_idx_name="Transformation & Storage",
            links_idx_name=DataModel.BUS_CARRIER,
            start=start,
        )
        # right = right.groupby(DataModel.BUS_CARRIER).sum()

        link = left.merge(
            right,
            on=DataModel.BUS_CARRIER,
            how="inner",
            suffixes=("_left", "_right"),
        )

        self.extend_view("link", "source", link["node_id_left"].to_numpy())
        self.extend_view("link", "target", link["node_id_right"].to_numpy())
        self.extend_view("link", "value", link["value_right"].to_numpy())
        self.extend_view("link", "label", link["link_label_left"].to_numpy())

        # left side already exists in nodes
        right_labels = right.set_index("node_id")["node_label"].drop_duplicates()
        self.extend_view("node", "label", right_labels.to_numpy())

        # todo: must map remainder from left to additional
        #  nodes on the right
        remainder = left["value"].sub(right["value"], fill_value=0)
        remainder.clip(lower=0)  # todo: remove me once data is correct

        passthrough = self.prepare_link_side(
            remainder.to_frame("value"),
            nodes_idx_name=DataModel.BUS_CARRIER,
            links_idx_name=DataModel.BUS_CARRIER,
            start=self.get_start_node_id(),
        )

        link = left.merge(
            passthrough,
            on=DataModel.BUS_CARRIER,
            how="inner",
            suffixes=("_left", "_right"),
        )

        link = link[["node_id_left", "node_id_right", "value_right", "link_label_left"]]
        self.extend_view("link", "source", link["node_id_left"].to_numpy())
        self.extend_view("link", "target", link["node_id_right"].to_numpy())
        self.extend_view("link", "value", link["value_right"].to_numpy())
        self.extend_view("link", "label", link["link_label_left"].to_numpy())

        # left side already exists in nodes
        right_labels = passthrough.set_index("node_id")["node_label"].drop_duplicates()
        self.extend_view("node", "label", right_labels.to_numpy())

    def prepare_link_side(
        self,
        metric_slice: pd.DataFrame,
        nodes_idx_name: str,
        links_idx_name: str,
        start: int,
    ) -> pd.DataFrame:
        """
        Add ids and labels to an input dataframe.

        Set

        Parameters
        ----------
        metric_slice
        nodes_idx_name
        links_idx_name
        start

        Returns
        -------
        :
            A copy of the input dara frame with additional
            columns for node and link identifier and display
            names.
        """
        # df = self.model[metric].copy()

        # set labels to values from a given index level,
        # or to an arbitrary string
        if nodes_idx_name in metric_slice.index.names:
            enumerated_labels = enumerate(
                metric_slice.index.unique(nodes_idx_name), start=start
            )
            node_ids = {carrier: i for i, carrier in enumerated_labels}
            metric_slice["node_id"] = metric_slice.index.get_level_values(
                nodes_idx_name
            ).map(node_ids)
            metric_slice["node_label"] = metric_slice.index.get_level_values(
                nodes_idx_name
            )
        else:
            metric_slice["node_id"] = start
            metric_slice["node_label"] = nodes_idx_name

        if links_idx_name in metric_slice.index.names:
            metric_slice["link_label"] = metric_slice.index.get_level_values(
                links_idx_name
            )
        else:
            metric_slice["link_label"] = links_idx_name

        return metric_slice

    def get_start_node_id(self) -> int:
        """
        Return the largest node id from the view.

        The returned value is used as a first value for node ids.
        The very first id is zero, subsequent ids are the largest
        existing id plus one.

        Returns
        -------
        :
            The largest node id value from the source and target link
            entries.
        """
        ids_source = self.view["link"]["source"]
        ids_target = self.view["link"]["target"]

        # edge case: empty arrays return 0
        max_source = int(ids_source.max()) + 1 if len(ids_source) > 0 else 0
        max_target = int(ids_target.max()) + 1 if len(ids_target) > 0 else 0

        return max(max_source, max_target)

    def extend_view(self, link_or_node: str, field: str, items: np.array) -> None:
        """
        Add given sources to the view.

        Parameters
        ----------
        link_or_node
            The name of the instance attribute to update.
        field
            The link dictionary key to update.
        items
            Integer source IDs to add to the view.

        Returns
        -------
        :
            Inplace updates view attributes.
        """
        value = self.view[link_or_node][field]
        self.view[link_or_node][field] = np.concatenate([value, items])


def read_iamc_data_frame():
    xls = pd.read_excel(
        "/IdeaProjects/pypsa-at/results/v2025.02/KN2045_Mix/evaluation/exported_iamc_variables.xlsx",
        index_col=[0, 1, 2, 3, 4],
    )
    xls.columns.name = "Year"
    return xls.stack()


def main():
    year = "2050"
    region = "GB0"

    df = read_iamc_data_frame()
    df = rename_aggregate(df, "TWh", level="Unit").div(1e6)
    df = filter_by(df, Year=year, Region=region)

    # ['AC',
    #  'Biomass',
    #  'Coal',
    #  'Gas',
    #  'H2',
    #  'Heat',
    #  'Methanol',
    #  'NH3',
    #  'Oil',
    #  'Uranium',
    #  'Waste']

    sankey_mapping = {
        # AC
        "Primary Energy|AC|Import Domestic": ("Domestic Import AC", "AC Primary"),
        "Primary Energy|AC|Import Foreign": ("Foreign Import AC", "AC Primary"),
        "Primary Energy|AC|Run-of-River": ("Run-of-River", "AC Primary"),
        "Primary Energy|AC|Solar Rooftop": ("Solar (Rooftop)", "AC Primary"),
        "Primary Energy|AC|Solar Utility": ("Solar (Utility)", "AC Primary"),
        "Primary Energy|AC|Wind Offshore": ("Wind (Offshore)", "AC Primary"),
        "Primary Energy|AC|Wind Onshore": ("Wind (Onshore)", "AC Primary"),
        "Secondary Energy|Demand AC": ("AC Primary", "Transformation Input"),
        "Secondary Energy|Bypass AC": (
            "AC Primary",
            "AC Bypass",
        ),  # bypass = final - lossy secondary supply
        # todo: not everything goes to transformation, some parts go to export and final demand
        #  need to connect secondary AC demand to Transformation Input! Not sum(primary).
        #  the rest is final AC demand (the bypass amount)
        # Primary Energy|Biomass|Solid
        # Primary Energy|Biogas
        # Primary Energy|Coal|Hard
        # Primary Energy|Coal|Lignite
        # Primary Energy|Gas|Biogas|w CC
        # Primary Energy|Gas|Biogas|w/o CC
        # Primary Energy|Gas|Import Domestic
        # Primary Energy|Gas|Import Foreign
        # Primary Energy|Gas|Production
        "Primary Energy|H2|Import Global": ("Global Import H2", "H2 Primary"),
        "Primary Energy|H2|Import Domestic": ("Domestic Import H2", "H2 Primary"),
        "Primary Energy|H2|Import Foreign": ("Foreign Import H2", "H2 Primary"),
        # "Secondary Energy|Demand H2": ("H2 Primary", "Transformation Input"),
        # "Secondary Energy|Bypass H2": ("H2 Primary", "H2 Bypass"),
        # Secondary Energy|Ambient Heat|AC|Air Heat Pump
        # Secondary Energy|Ambient Heat|AC|Ground Heat Pump
        # Secondary Energy|Ambient Heat|Biomass|CHP
        # Secondary Energy|Ambient Heat|Coal|CHP
        # Secondary Energy|Ambient Heat|Gas|Boiler
        # Primary Energy|Methanol
        # Primary Energy|NH3  # skip?
        # Primary Energy|Oil|Import
        # Primary Energy|Uranium
        # Primary Energy|Waste|HVC from naphtha
        # Primary Energy|Waste|Import Foreign
        # Primary Energy|Waste|Solid
    }
    variables = df.index.unique("Variable")
    for k in sankey_mapping:
        if k in variables:
            continue
        print(f"Skipping '{k}' because it does not exist in {year} and {region}.")
        sankey_mapping.pop(k)

    sankey_mapping = {
        k: v for k, v in sankey_mapping.items() if k in df.index.unique("Variable")
    }

    df = pyam.IamDataFrame(df)
    fig = df.plot.sankey(mapping=sankey_mapping)
    plotly.io.show(fig)


if __name__ == "__main__":
    FILEPATH = "/IdeaProjects/pypsa-at/results/v2025.02/KN2045_Mix/evaluation/exported_iamc_variables.xlsx"
    main()
