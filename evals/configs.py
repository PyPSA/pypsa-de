"""Module to collect configuration items and their default values."""

from dataclasses import dataclass, field

from evals.constants import COLOUR_SCHEME_BMK, DataModel, Group


@dataclass()
class PlotConfig:
    """Holds configuration items for Plotly figures."""

    title: str = None
    chart = None  # ESMBarChart | ESMGroupedBarChart | ESMTimeSeriesChart
    file_name_template: str = "{metric}_{year}_{location}"
    unit: str = ""  # default is metric.df.attrs["unit"]

    # the metric data frame is grouped this index level before plotting.
    # One html figure is created per resulting group.
    plotby: list = field(default_factory=lambda: [DataModel.LOCATION])

    # Used to pivot the data frame before sending it to the plotter. The
    # specified index/column levels will be in the plot data frame. The
    # rest is aggregated (summed up).
    pivot_index: list = field(default_factory=lambda: DataModel.YEAR_IDX_NAMES)
    pivot_columns: list = field(default_factory=lambda: [])

    plot_category: str = DataModel.CARRIER
    plot_xaxis: str = DataModel.YEAR

    # defines the subplots in GroupedBarChart
    facet_column: str = DataModel.BUS_CARRIER

    category_orders: tuple = ()
    colors: dict = field(default_factory=lambda: COLOUR_SCHEME_BMK)
    pattern: dict = field(
        default_factory=lambda: dict.fromkeys(
            [
                Group.import_foreign,
                Group.export_foreign,
                Group.import_domestic,
                Group.export_domestic,
                Group.import_net,
                Group.export_net,
                Group.import_global,
            ],
            "/",
        )
    )
    fill: dict = field(default_factory=dict)
    stacked: bool = True
    line_dash: dict = field(default_factory=dict)
    line_width: dict = field(default_factory=dict)
    line_shape: str = "hv"
    legend_header: str = "Categories"
    xaxis_title: str = "<b>Years</b>"
    yaxis_color: str = "DarkSlateGrey"
    footnotes: tuple = ("", "")
    cutoff: float = 0.0001  # needs update depending on unit
    cutoff_drop: bool = True  # only effective in BarCharts

    legend_font_size: int = 20
    title_font_size: int = 30
    font_size: int = 20
    xaxis_font_size: int = 20
    yaxes_showgrid: bool = False
    yaxes_visible: bool = False


@dataclass()
class ExcelConfig:
    """Holds configuration items for Excel file."""

    axis_labels: list = None
    chart: str = "stacked"  # 'stacked', 'clustered', 'standard', 'percentStacked', None
    chart_title: str = None
    chart_width: int = 20  # cm
    chart_switch_axis: bool = False  # switch categories with x-axis
    chart_colors: dict = field(
        default_factory=lambda: {k: v.lstrip("#") for k, v in COLOUR_SCHEME_BMK.items()}
    )
    # pivot tables to use the following labels as index or column
    pivot_index: str | list = field(
        default_factory=lambda: [DataModel.LOCATION, DataModel.CARRIER]
    )
    pivot_columns: str | list = DataModel.YEAR


@dataclass()
class ViewDefaults:
    """
    Holds all configuration items needed to export Metrics.

    The 'excel' and 'plotly' fields are processed by the export_excel
    and export_plotly methods, respectively. Both configuration spaces
    are kept separate to keep the variable space small during export.
    """

    excel: ExcelConfig = field(default_factory=lambda: ExcelConfig())
    plotly: PlotConfig = field(default_factory=lambda: PlotConfig())
