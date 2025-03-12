# -*- coding: utf-8 -*-
import logging
import os
import sys

sys.path.append(os.path.abspath(os.path.dirname(__file__)))  # Adds 'scripts/' to path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))  # Adds repo root

import pypsa
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from datetime import datetime
import pickle
from pypsa.descriptors import get_switchable_as_dense as as_dense
from _helpers import configure_logging, mock_snakemake
from pricing_analysis import timestep_before, get_supply_demand, supply_price_link, demand_price_link, get_compressed_demand

logger = logging.getLogger(__name__)

year_colors = ["dimgrey", "darkorange", "seagreen", "cadetblue", "hotpink", "darkviolet"]
markers =["v", "^", "<", ">", "1", "2", "3", "4", "*", "+", "d", "o", "|", "s", "P", "p", "h"]
date_format = "%Y-%m-%d %H:%M:%S"

carrier_renaming = {
    'urban central solid biomass CHP CC': 'biomass CHP CC',
    'urban central solid biomass CHP': 'biomass CHP',
    'urban central gas CHP': 'gas CHP',
    'urban central gas CHP CC': 'gas CHP CC',
    'urban central air heat pump': 'air heat pump',
    'urban central resistive heater': 'resistive heater',
    'urban central lignite CHP': 'lignite CHP'
}

carrier_renaming_reverse = {
    'biomass CHP CC': 'urban central solid biomass CHP CC',
    'biomass CHP' :'urban central solid biomass CHP' ,
    'gas CHP': 'urban central gas CHP' ,
    'gas CHP CC' : 'urban central gas CHP CC',
    'air heat pump' : 'urban central air heat pump',
    'resistive heater': 'urban central resistive heater',
    'lignite CHP': 'urban central lignite CHP'
}

carrier_groups = {
    'CCGT': 'gas turbines',
    'OCGT': 'gas turbines',
    'coal': 'coal incl. CHP',
    'urban central coal CHP': 'coal incl. CHP',
    'urban central lignite CHP': 'lignite incl. CHP',
    'urban central gas CHP': 'gas CHP',
    'urban central gas CHP CC': 'gas CHP',
    'onwind': 'solar+wind',
    'offwind-ac': 'solar+wind',
    'offwind-dc': 'solar+wind',
    'solar': 'solar+wind',
    'solar-hsat': 'solar+wind',
    'urban central solid biomass CHP': 'biomass CHP',
    'biogas': 'gas turbines',
    'H2 OCGT': 'H2 turbines',
    'H2 retrofit OCGT': 'H2 turbines',
    'urban central oil CHP': 'oil inlc. CHP',
    'lignite': 'lignite incl. CHP',
    'waste CHP': 'waste CHP',
    'waste CHP CC': 'waste CHP',

}

group_colors = {
    'coal incl. CHP': '#545454',
    'lignite incl. CHP': '#826837',
    'gas turbines': '#e05b09',
    'gas CHP': 'darkred',
    'solar+wind': 'lightgreen',
    'biomass CHP': 'forestgreen',
    'PHS': '#51dbcc',
    'battery discharger': 'darkviolet',
    'hydro': '#298c81',
    'H2 turbines': 'darkblue',
    'nuclear': '#ff8c00',
    'oil inlc. CHP': '#c9c9c9',
    'waste CHP': '#e3d37d',
    'resistive heater': 'indianred',
    'air heat pump': 'salmon',
    'ground heat pump': 'firebrick'
 }

carrier_groups_d = {
    'urban central resistive heater': 'resistive heater',
    'urban decentral resistive heater': 'resistive heater',
    'rural resistive heater': 'resistive heater',
    'rural air heat pump': 'air heat pump',
    'urban central air heat pump': 'air heat pump',
    'urban decentral air heat pump': 'air heat pump',
    'rural ground heat pump': 'ground heat pump',
}


group_colors_demand = {
    'resistive heater': 'indianred',
    'air heat pump': 'salmon',
    'ground heat pump': 'firebrick',
}

def plot_supply_demand(n, 
                    supply, 
                    demand, 
                    buses, 
                    timestep,
                    tech_colors, 
                    ylim=None, 
                    p="p_nom_opt", 
                    d="p_nom_opt", 
                    mc="mc", 
                    savepath=None,
                    only_carriers=False, 
                    whole_system=False, 
                    demand_plot=True, 
                    demand_text=True, 
                    compress_demand=False, 
                    year=9999):
    # Filter out technologies with negative supply
    drop_c = ["electricity distribution grid", "BEV charger"]
    supply = supply[(supply[p] >= 1) & ~(supply.carrier.isin(drop_c))]

    # Convert to GW
    supply.loc[:, p] = supply[p] / 1e3
    
    if whole_system:
        supply = supply[~supply.carrier.isin(["AC", "DC"])]
        demand = demand[~demand.carrier.isin(["AC", "DC"])]

    # Sort the supply DataFrame by the marginal cost column
    supply = supply.sort_values(by=mc)

    # Process the demand data
    demand = demand[demand[d] > 1]  # Drop all technologies with demand less than 1e-1 MW
    demand.loc[:,d] = demand[d] / 1e3  # Convert to GW
    demand.loc[:,"bidding_price"] = demand.bidding_price.clip(lower=0) # Clip negative bidding prices to 0 
    # group demand technologies together, where the bidding price differene is lower than a threshold
    if compress_demand:
        demand = get_compressed_demand(demand, th=0.1)
        only_carriers = True # no more index available as data is grouped
    demand = demand.sort_values(by='bidding_price', ascending=False)
    
    # Add start and end of demand curve for plotting
    start_row = pd.DataFrame([{"bidding_price": np.inf, "p": 0, "p_nom_opt": 0, "volume_demand": 0, "carrier": "start"}], index=["start"])
    end_row = pd.DataFrame([{"bidding_price": 0, "p": 0, "p_nom_opt": 0,"volume_demand": 0, "carrier": "end"}], index=["end"])
    demand = pd.concat([start_row, demand, end_row])
    
    # Replace infinite bidding price with 10% more than the maximum finite bidding price
    max_finite_bidding_price = demand.bidding_price[demand.bidding_price.apply(np.isfinite)].max()
    demand["bidding_price"] = demand["bidding_price"].replace(np.inf, max_finite_bidding_price * 1.1)
    
    # Cumulative sum of demand
    demand.loc[:,d] = demand[d].cumsum()

    if compress_demand:
        demand.reset_index(drop=True, inplace=True)

    # Create the plot
    fig, ax = plt.subplots(figsize=(8, 6))

    x_position = 0

    supply_labels = set()  # Set to track labels to avoid duplicates

    # Plot supply bars
    for index, row in supply.iterrows():
        height = row[mc]
        width = row[p]
        carrier_label = str(row['carrier']) if only_carriers else index
        ax.add_patch(plt.Rectangle((x_position, 0), width, height, color=tech_colors[str(row['carrier'])], 
                                   alpha=0.5, label=carrier_label if carrier_label not in supply_labels else "_nolegend_"))
        supply_labels.add(carrier_label)
        x_position += width

    ax.set_xlim(0, max(x_position, demand[d].max()) * 1.11)
    ax.set_ylim(0, max(demand.bidding_price.max(), supply[mc].max()) * 1.11)
    if ylim is not None:
        ax.set_ylim(0, ylim)

    plt.xlabel('bid and ask volume [GW]')
    plt.ylabel('bid and ask price [€/MWh]')
    plt.title(f'Market clearing for electricity at {timestep}')

    demand_legend = {}  # Dictionary to store the number-label pairs for the legend

    if demand_plot:
        # Plot demand curve
        for i in range(len(demand) - 1):
            x1 = demand.iloc[i][d]
            x2 = demand.iloc[i + 1][d]
            y1 = demand.iloc[i + 1]['bidding_price']
            y2 = demand.iloc[i + 1]['bidding_price']
            ax.plot([x1, x2], [y1, y2], color='black', linestyle='--')
            if (demand_text & (i < len(demand) - 2)):
                number_label = str(i + 1)  # Number to display in the plot
                ax.text((x1 + x2) / 2, ((y1 + y2) / 2)+ 7, number_label, fontsize=12, ha='center', rotation=0)
                demand_legend[number_label] = (demand.index[i + 2] if not only_carriers else demand.carrier[i + 1])
            ax.plot([x1, x1], [demand.iloc[i]['bidding_price'], y2], color='grey', linestyle='--')


    # Create supply legend handles and labels
    supply_handles, supply_labels = ax.get_legend_handles_labels()
    supply_handles_labels = dict(zip(supply_labels, supply_handles))
    supply_handles_labels.pop("_nolegend_", None)  # Remove the dummy label used to avoid duplicates

    # Create demand legend handles and labels
    demand_legend_handles = [plt.Line2D([0], [0], color='white', label=f"{num}: {label}") for num, label in demand_legend.items()]

    # Plot market clearing point if there is only one bus
    if len(buses) == 1:
        ax.plot(n.statistics.supply(bus_carrier="AC", aggregate_time=False)[timestep].sum()/1e3, n.buses_t.marginal_price.loc[timestep, buses], 
                 marker='x', markersize=7, color="red", label='market clearing')


    legend1 = plt.legend(handles=supply_handles_labels.values(), labels=supply_handles_labels.keys(), title="Supply",
                        loc='upper left', bbox_to_anchor=(1, 1.02))
    ax.add_artist(legend1)

    legend2 = plt.legend(handles=demand_legend_handles, title="Demand", loc='upper left', bbox_to_anchor=(0, -0.12), ncol=2)
    ax.add_artist(legend2)
    # plt.tight_layout()

    plt.grid(True)
    fig.savefig(savepath, bbox_extra_artists=(legend1,legend2), bbox_inches='tight')

def get_compressed_demand(demand, th):
    demand.loc[:, "carrier"] = demand["carrier"].replace({
        "electricity": "electricity load", 
        "industry electricity": "electricity load", 
        "agriculture electricity": "electricity load"
    })
    df = demand.sort_values(by=['carrier', 'bidding_price']).reset_index(drop=True)
    group = (df['bidding_price'].diff().abs() > th).cumsum()
    df.loc[:,'group'] = group

    grouped_df = df.groupby(['carrier', 'group']).agg({
        'bidding_price': 'mean',  
        'p': 'sum',               
        'p_nom_opt': 'sum',       
        'volume_demand': 'sum',    
        'carrier': 'first',  
    }).reset_index(drop=True)

    return grouped_df

    # side by side merit order plot
def plot_supply_demand_s(n, supply, demand, buses, timestep, ylim=None, p="p_nom_opt", d="p_nom_opt", mc="mc", 
                       only_carriers=False, whole_system=False, demand_plot=True, demand_text=True, 
                       compress_demand=False, ax=None, year=9999):
    # Filter out technologies with negative supply
    drop_c = ["electricity distribution grid", "BEV charger"]
    supply = supply[(supply[p] >= 5) & ~(supply.carrier.isin(drop_c))]
    supply.loc[:,p] = supply[p] / 1e3 # Convert to GW
    
    if whole_system:
        supply = supply[~supply.carrier.isin(["AC", "DC"])]
        demand = demand[~demand.carrier.isin(["AC", "DC"])]
    
    supply = supply.sort_values(by=mc)
    demand = demand[demand[d] > 1]
    demand.loc[:,d] = demand[d] / 1e3
    demand.loc["bidding_price"] = demand.bidding_price.clip(lower=0) # Clip negative bidding prices to 0
    
    if compress_demand:
        demand = get_compressed_demand(demand, th=0.1)
        only_carriers = True

    demand = demand.sort_values(by='bidding_price', ascending=False)
    
    start_row = pd.DataFrame([{"bidding_price": np.inf, "p": 0, "p_nom_opt": 0, "volume_demand": 0, "carrier": "start"}], index=["start"])
    end_row = pd.DataFrame([{"bidding_price": 0, "p": 0, "p_nom_opt": 0, "volume_demand": 0, "carrier": "end"}], index=["end"])
    demand = pd.concat([start_row, demand, end_row])
    
    max_finite_bidding_price = demand.bidding_price[demand.bidding_price.apply(np.isfinite)].max()
    demand["bidding_price"] = demand["bidding_price"].replace(np.inf, max_finite_bidding_price * 1.1)
    demand[d] = demand[d].cumsum()

    if compress_demand:
        demand.reset_index(drop=True, inplace=True)
    
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 6))

    x_position = 0
    supply_labels = set()
    
    for index, row in supply.iterrows():
        height = row[mc]
        width = row[p]
        carrier_label = str(row['carrier']) if only_carriers else index
        ax.add_patch(plt.Rectangle((x_position, 0), width, height, color=tech_colors[str(row['carrier'])], 
                                   alpha=0.5, label=carrier_label if carrier_label not in supply_labels else "_nolegend_"))
        supply_labels.add(carrier_label)
        x_position += width

    ax.set_xlim(0, max(x_position, demand[d].max()) * 1.11)
    ax.set_ylim(0, max(demand.bidding_price.max(), supply[mc].max()) * 1.11)
    if ylim is not None:
        ax.set_ylim(0, ylim)
        
    if len(buses) == 1:
        ax.plot(n.statistics.supply(bus_carrier="AC", aggregate_time=False)[timestep].sum(), 
                n.buses_t.marginal_price.loc[timestep, buses], marker='x', markersize=7, color="red", label='market clearing')
        
    ax.set_xlabel('bid and ask volume [GW]')
    ax.set_ylabel('bid and ask price [€/MWh]')
    ax.set_title(f'Market clearing for electricity at {timestep} in year {year}')
    
    demand_legend = {}

    if demand_plot:
        for i in range(len(demand) - 1):
            x1 = demand.iloc[i][d]
            x2 = demand.iloc[i + 1][d]
            y1 = demand.iloc[i + 1]['bidding_price']
            y2 = demand.iloc[i + 1]['bidding_price']
            ax.plot([x1, x2], [y1, y2], color='black', linestyle='--')
            if (demand_text & (i < len(demand) - 2)):
                number_label = str(i + 1)
                ax.text((x1 + x2) / 2, ((y1 + y2) / 2) + 7, number_label, fontsize=12, ha='center', rotation=0)
                demand_legend[number_label] = (demand.index[i + 2] if not only_carriers else demand.carrier[i + 1])
            ax.plot([x1, x1], [demand.iloc[i]['bidding_price'], y2], color='grey', linestyle='--')

    supply_handles, supply_labels = ax.get_legend_handles_labels()
    supply_handles_labels = dict(zip(supply_labels, supply_handles))
    supply_handles_labels.pop("_nolegend_", None)

    demand_legend_handles = [plt.Line2D([0], [0], color='white', label=f"{num}: {label}") for num, label in demand_legend.items()]

    # Plot market clearing point if there is only one bus
    if len(buses) == 1:
        ax.plot(n.statistics.supply(bus_carrier="AC", aggregate_time=False)[timestep].sum()/1e3, n.buses_t.marginal_price.loc[timestep, buses], 
                 marker='x', markersize=7, color="red", label='market clearing')

    
    supply_handles, supply_labels = ax.get_legend_handles_labels()
    supply_handles_labels = dict(zip(supply_labels, supply_handles))
    supply_handles_labels.pop("_nolegend_", None)

    legend2 = ax.legend(handles=demand_legend_handles, title="Demand", loc='upper right', ncol=2, fancybox=True, framealpha=0.3)
    ax.add_artist(legend2)

    ax.grid(True)

    # Return the legend handles and labels
    return supply_handles_labels

def get_all_supply_prices(n, bus, period=None, carriers=None):

    index = n.snapshots if period is None else period

    # Initialize a dictionary to store results temporarily
    res_dict = {ts: {} for ts in index}

    for gen in n.generators.index[(n.generators.bus == bus) & ((n.generators.carrier.isin(carriers)) if carriers is not None else True)]:
        marginal_cost = n.generators.loc[gen].marginal_cost
        for ts in index:
            res_dict[ts][gen] = marginal_cost

    for su in n.storage_units.index[(n.storage_units.bus == bus) & ((n.storage_units.carrier.isin(carriers)) if carriers is not None else True)]:
        for ts in index:
            res_dict[ts][su] = (
                n.storage_units.loc[su].marginal_cost +
                n.storage_units_t.mu_energy_balance.loc[ts, su] *
                1 / n.storage_units.efficiency_dispatch.loc[su]
            )

    for st in n.stores.index[(n.stores.bus == bus) & ((n.stores.carrier.isin(carriers)) if carriers is not None else True)]:
        for ts in index:
            res_dict[ts][st] = (
                n.stores.loc[st].marginal_cost +
                n.stores_t.mu_energy_balance.loc[ts, st]
            )

    loc_buses = ["bus" + str(i) for i in np.arange(0, 5)]
    for link in n.links.index[(n.links.bus0 != bus) & (n.links[loc_buses].isin([bus]).any(axis=1)) & ((n.links.carrier.isin(carriers)) if carriers is not None else True)]:
        for ts in index:
            res_dict[ts][link] = supply_price_link(n, link, ts, bus)

    # Convert the dictionary to a DataFrame
    res = pd.DataFrame.from_dict(res_dict, orient='index')
    
    return res


if __name__ == "__main__":
    if "snakemake" not in globals():
        import os
        import sys
        from _helpers import mock_snakemake

        snakemake = mock_snakemake(
            "pricing_plots",
            simpl="",
            clusters=1,
            opts="",
            ll="vopt",
            sector_opts="None",
            run="KN2045_Bal_v4_24H",
        )
    
    # ensure output directory exist
    for dir in snakemake.output[2:]:
        if not os.path.exists(dir):
            os.makedirs(dir)

    configure_logging(snakemake)
    config = snakemake.config
    planning_horizons = snakemake.params.planning_horizons
    nhours = int(snakemake.params.hours[:-1])
    nyears = nhours / 8760
    tech_colors = snakemake.params.plotting["tech_colors"]

    # Load networks
    networks = [pypsa.Network(fn) for fn in snakemake.input.networks]
    modelyears = [fn[-7:-3] for fn in snakemake.input.networks]

    # save as dict
    n_dict = {}
    years = np.arange(2020,2050,5)
    for i, year in enumerate(years):
        n_dict[year] = networks[i]
    networks = n_dict

    # Load price setting data
    with open(snakemake.input.pricing_s, 'rb') as file:
        results_s = pickle.load(file)
    with open(snakemake.input.pricing_d, 'rb') as file:
        results_d = pickle.load(file)

    # convert to datetime
    for year in years:
        results_s[year].timestep = pd.to_datetime(results_s[year].timestep, format=date_format)
        results_d[year].timestep = pd.to_datetime(results_d[year].timestep, format=date_format)
     

    # update tech_colors
    colors_update = networks[2020].carriers.color.rename(networks[2020].carriers.nice_name).to_dict()
    colors_update = {k: v for k, v in colors_update.items() if v != ""}
    tech_colors.update(colors_update)

    # carrier names manual
    tech_colors["urban central oil CHP"] = tech_colors["oil"]
    tech_colors["Solar"] = tech_colors["solar"]
    tech_colors["Electricity load"] = tech_colors["electricity"]
    tech_colors["Electricity trade"] = tech_colors["AC"]
    tech_colors["Offshore Wind"] = tech_colors["offwind-ac"]
    tech_colors["urban decentral heat"] = tech_colors["urban central heat"]
    tech_colors["urban decentral biomass boiler"] = tech_colors["biomass boiler"]
    tech_colors["rural biomass boiler"] = tech_colors["biomass boiler"]
    tech_colors["urban decentral oil boiler"] = tech_colors["oil boiler"]
    tech_colors["rural oil boiler"] = tech_colors["oil boiler"]
    tech_colors["rural ground heat pump"] = tech_colors["ground heat pump"]
    tech_colors['gas CHP'] = 'darkorange'
    tech_colors['urban decentral resistive heater'] = 'indianred' 
    tech_colors['urban decentral air heat pump'] = 'salmon'
    tech_colors['rural resistive heater'] = "indianred"
    tech_colors['rural air heat pump'] = "salmon"
    tech_colors["H2 OCGT"] = tech_colors["H2"]
    tech_colors["H2 retrofit OCGT"] = tech_colors["H2"]
    tech_colors["urban central H2 retrofit CHP"] = "turquoise"
    tech_colors["battery discharger"] = "darkgoldenrod"
    tech_colors["urban central H2 retrofit OCGT"] = "seagreen"
    tech_colors["load-shedding"] = "slategrey"

    # plotting - merit order 3 cases
    n = networks[2020]
    if "2019-01-11 15:00:00" in  n.snapshots:
        ts = ["2019-01-11 15:00:00", "2019-08-05 18:00:00",  "2019-06-02 12:00:00"] 
    elif "2013-02-18 15:00:00" in  n.snapshots:
        ts = ["2013-02-18 15:00:00", "2013-12-31 15:00:00", "2013-07-07 12:00:00"]
    else: 
        ts = n.snapshots[[0, 1, -1]]

    for year in planning_horizons:

        all_supply_handles_labels = {}
        num_subplots = 3
        fig, axes = plt.subplots(3, 1, figsize=(8, 3*6))
        axes = axes.flatten()

        # Plot each subplot and collect legend handles and labels
        for i in range(num_subplots):
            n = networks[year]
            buses = ["DE0 0"]
            timestep = str(ts[i])
            supply, demand = get_supply_demand(n, buses, timestep)
            supply_handles_labels = \
                plot_supply_demand_s(n, supply, demand, buses, timestep, p="volume_bid", d="volume_demand", mc="mc_final",
                                only_carriers=True, demand_text=True, compress_demand=True, ax=axes[i], year=year)

            # Merge the current subplot's supply legend with the global collection
            all_supply_handles_labels.update(supply_handles_labels)
        del all_supply_handles_labels["market clearing"]    

        # Create the combined supply legend
        fig.legend(handles=all_supply_handles_labels.values(),
                labels=all_supply_handles_labels.keys(),
                title="Supply", loc='lower center', ncol=3, bbox_to_anchor=(0.5, -0.1))

        plt.tight_layout(pad=2)
        plt.savefig(f"{snakemake.output.merit_order_3cases}/{year}.png", bbox_inches='tight')
        plt.close(fig) 

    # plotting - merit orders
    for year in planning_horizons:
        n = networks[year]
        # plot only every 100th timestep
        for timestep in n.snapshots[::5]:
            supply, demand = get_supply_demand(n, buses, str(timestep))
            plot_supply_demand(n, 
                    supply, 
                    demand, 
                    buses, 
                    str(timestep), 
                    tech_colors,
                    ylim=None, 
                    p="volume_bid", 
                    d="volume_demand", 
                    mc="mc_final",
                    savepath=f"{snakemake.output.merit_order_all}/{year}-{timestep}.png",
                    only_carriers=True, 
                    whole_system=False, 
                    demand_plot=True, 
                    demand_text=True, 
                    compress_demand=True, 
                    year=year)

    # plotting - price setter temporal (supply)
    data = results_s
    markers_here = [6 , 7]

    for year, ylim in zip(planning_horizons , [(-10, 100) , (-10, 150) , (-10, 200), (-10 , 300) , (-10 , 400) , (-10 , 500)]):
        df = data[year][data[year].valid].copy()
        df.set_index("timestep", inplace=True)

        fig, ax = plt.subplots(figsize=(8, 6))
        for carrier in df.carrier.unique():
            ax.plot(
                df[df.carrier == carrier].index, 
                df["marginal price @ bus"][df.carrier == carrier],
                marker=markers_here[len(carrier) % 2],
                markersize=7, 
                linestyle="", 
                label=carrier, 
                color=tech_colors[carrier]
                )
        plt.ylim(ylim)
        plt.title(f"Price setter for electricity in {year} (supply)")
        plt.legend(bbox_to_anchor=(1, 1))
        plt.savefig(f"{snakemake.output.price_setter}/{year}-price-setter.png", bbox_inches='tight')

    
    # plotting - price taker temporal (demand)
    data = results_d

    for year, ylim in zip(planning_horizons , [(-10, 100) , (-10, 150) , (-10, 200), (-10 , 300) , (-10 , 400) , (-10 , 500)]):
        df = data[year][data[year].valid].copy()
        df.set_index("timestep", inplace=True)

        fig, ax = plt.subplots(figsize=(8, 6))
        for carrier in df.carrier.unique():
            ax.plot(
                df[df.carrier == carrier].index, 
                df["marginal price @ bus"][df.carrier == carrier],
                marker=markers_here[len(carrier) % 2],
                markersize=7, 
                linestyle="", 
                label=carrier, 
                color=tech_colors[carrier]
                )
        plt.ylim(ylim)
        plt.title(f"Least price taker for electricity in {year} (demand)")
        plt.legend(bbox_to_anchor=(1, 1))
        plt.savefig(f"{snakemake.output.price_taker}/{year}-price-taker.png", bbox_inches='tight')

    
    # plotting - market clearing price duration curve (price setter)
    data = results_s

    for year in planning_horizons:
        df = data[year][data[year].valid].copy()
        # select only every 5th row
        df = df.iloc[::5, :]
        df.sort_values(by="marginal price @ bus", ascending=False, inplace=True)
        df.reset_index(inplace=True)
        
        fig, ax = plt.subplots(figsize=(10, 6))
        for carrier in df.carrier.unique():
            plt.plot(df[df.carrier == carrier]["marginal price @ bus"], 
                    marker='|', 
                    linestyle="", 
                    label=carrier, 
                    color=tech_colors[carrier])

        plt.ylim(-10, 250)
        plt.title(f"Price duration curve for electricity by price setter in {year}")
        plt.xlabel("Hours of the year in 3h steps (sorted)")
        plt.ylabel("Market clearing price in EUR/MWh)")
        plt.legend(bbox_to_anchor=(1, 1),fancybox=True, shadow=True, ncol=1)
        plt.savefig(f"{snakemake.output.pdc_price_setter}/{year}-pdc-price-setter.png", bbox_inches='tight')

    # plotting - market clearing price duration curve (price taker)
    data = results_d

    for year in years:
        df = data[year][data[year].valid].copy()
        df = df[df.bidding_price.notna()]
        df.sort_values(by="marginal price @ bus", ascending=False, inplace=True)
        df.reset_index(inplace=True)
        
        fig, ax = plt.subplots(figsize=(10, 6))

        for carrier in df.carrier.unique():
            plt.plot(df[df.carrier == carrier]["marginal price @ bus"], "*", label=carrier, color=tech_colors[carrier])

        plt.ylim(-10, 200)
        plt.title(f"Price duration curve for electricity by highest price taker in {year}")
        plt.xlabel("Hours of the year in 3h steps (sorted)")
        plt.ylabel("Market clearing price in EUR/MWh)")
        plt.legend(bbox_to_anchor=(1, 1),fancybox=True, shadow=True, ncol=1)
        plt.savefig(f"{snakemake.output.pdc_price_taker}/{year}-pdc-price-taker.png", bbox_inches='tight')

    
    # plotting - price duration curves
    # Fraction of time [%]
    bus = "DE0 0" 
    fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(6, 5))

    for i, n in enumerate(networks.values()):
        lmps = pd.DataFrame(n.buses_t.marginal_price[bus])
        lmps.sort_values(by=bus, ascending=False, inplace=True)
        lmps["percentage"] = np.arange(len(lmps)) / len(lmps) * 100
        ax.plot(lmps["percentage"], lmps[bus], label=years[i], color=year_colors[i])

        ax.set_ylim([-50, 300])
        ax.set_ylabel("Electricity Price [$€/MWh_{el}$")
        ax.set_xlabel("Fraction of time [%]")
        ax.set_title(f"Electricity price duration curves", fontsize=16)
        ax.legend()
        ax.grid(True)

    fig.tight_layout()
    plt.savefig(snakemake.output.elec_pdc, bbox_inches='tight')

    
    # plotting - price setting development

    # supply
    data = results_s
    carriers = []
    for year in planning_horizons:
        df = data[year][data[year]["valid"] == True].copy()
        carriers.extend(df.carrier.unique().tolist())
    carriers = list(set(carriers))

    df_res = pd.DataFrame(index = planning_horizons, columns = carriers)
    df_res_absolut = pd.DataFrame(index = planning_horizons, columns = carriers)

    for year in planning_horizons:
        df = data[year][data[year]["valid"] == True].copy()
        df = df[df.supply_price.notna()]
        res = df.carrier.value_counts() / df.carrier.value_counts().sum()
        res_absolut = df.carrier.value_counts()
        df_res.loc[year] = res[res > 0.03]
        df_res_absolut.loc[year] = res_absolut[res_absolut > 100]

    df_s_grouped = df_res_absolut.T.groupby(carrier_groups).sum().T
    df_s_grouped = pd.concat([df_s_grouped, df_res_absolut.loc[:, ~df_res_absolut.columns.isin(carrier_groups.keys())]], axis=1)
    df_s_grouped = df_s_grouped[df_s_grouped.sum().sort_values(ascending=False).index]

    # demand
    data = results_d
    carriers = []
    for year in planning_horizons:
        df = data[year][data[year]["valid"]].copy()
        carriers.extend(df.carrier.unique().tolist())
    carriers = list(set(carriers))

    df_res = pd.DataFrame(index = planning_horizons, columns = carriers)
    df_res_absolut = pd.DataFrame(index = planning_horizons, columns = carriers)

    for year in planning_horizons:
        df = data[year][data[year]["valid"]].copy()
        res = df.carrier.value_counts() / df.carrier.value_counts().sum()
        res_absolut = df.carrier.value_counts()
        df_res.loc[year] = res[res > 0]
        df_res_absolut.loc[year] = res_absolut[res_absolut > 10]
    
    group_colors.update(group_colors_demand)

    df_d_grouped = df_res_absolut.T.groupby(carrier_groups_d).sum().T
    df_d_grouped = pd.concat([df_d_grouped, df_res_absolut.loc[:, ~df_res_absolut.columns.isin(carrier_groups_d.keys())]], axis=1)
    df_d_grouped = df_d_grouped[df_d_grouped.sum().sort_values(ascending=False).index]

    # combined plot
    df_all_grouped = pd.concat([df_s_grouped, df_d_grouped], axis=1)
    df_all_grouped = df_all_grouped.T.groupby(df_all_grouped.columns).sum().T
    df_all_grouped = df_all_grouped.loc[:, ~(df_all_grouped == 0).all()]

    res = pd.DataFrame(index=planning_horizons, columns=["Only supply price setter", "Only demand price setter", "Supply and demand price setter"])
    for year in planning_horizons:
        df_s = results_s[year].copy().set_index("timestep")
        df_d = results_d[year].copy().set_index("timestep")
        df_d = df_s.reset_index().drop_duplicates(subset=['timestep'], keep='first').set_index('timestep')
        df_s = df_d.reset_index().drop_duplicates(subset=['timestep'], keep='first').set_index('timestep')
        both_i = df_s[(df_s["valid"] & df_d["valid"])].index
        only_s_i = df_s[(df_s["valid"] & ~df_d["valid"])].index
        only_d_i = df_d[(df_d["valid"] & ~df_s["valid"])].index
        res.loc[year] = [len(only_s_i), len(only_d_i), len(both_i)]
    
    # Plot the data
    fig, ax = plt.subplots(figsize=(8, 6))

    color_list = [group_colors.get(carrier, tech_colors.get(carrier)) for carrier in df_all_grouped.columns]
    df_all_grouped.plot(kind="area", stacked=True, figsize=(8, 6), color=color_list, ax=ax)
    res.plot(ax=ax, color=["black", "grey", "brown"], linestyle="--", marker="o")
    plt.title("Development of Price setting Technology for Electricity (total)")
    plt.xlabel("Year")
    plt.ylabel("Price Setter [%]")
    plt.legend(bbox_to_anchor=(1, -0.1), fancybox=True, shadow=True, ncol=4)
    plt.savefig(snakemake.output.price_setting_dev, bbox_inches='tight')
