# PyPSA-DE - Hochaufgelöstes, sektorengekoppeltes Modell des deutschen Energiesystems

PyPSA-DE ist ein sektorengekoppeltes Energiesystem-Modell auf Basis der Toolbox [PyPSA](https://github.com/PyPSA/pypsa) und des europäischen Modells [PyPSA-Eur](https://github.com/PyPSA/pypsa-eur). Der PyPSA-DE Workflow modelliert das deutsche Energiesystem mit deutschlandspezifischen Datensätzen (MaStR, Netzentwicklungsplan,...) im Verbund mit den direkten Stromnachbarn sowie Spanien und Italien. Der Ausbau und der Betrieb von Kraftwerken, des  Strom- und Wasserstoffübertragunsnetzes und die Energieversorgung aller Sektoren werden dann in einem linearen Optimierungsproblem gelöst, mit hoher zeitlicher und räumlicher Auflösung. PyPSA-DE wurde im Rahmen des Kopernikus-Projekts [Ariadne](https://ariadneprojekt.de/) entwickelt in dem Szenarien für ein klimaneutrales Deutschland untersucht werden, und spielt eine zentrale Rolle im [Ariadne Szenarienreport](https://ariadneprojekt.de/publikation/report-szenarien-zur-klimaneutralitat-2045/), als Leitmodell für den [Sektor Energiewirtschaft und Infrastruktur](https://ariadneprojekt.de/publikation/report-szenarien-zur-klimaneutralitat-2045/#6-sektorale-perspektive-energiewirtschaft) und als eines von drei Gesamtsystemmodellen. Die Ergebnisse aus der Modellierung mit PyPSA-DE werden auch im [Ariadne-Webinar zu den Kernaussagen des Berichts](https://youtu.be/UL3KAH7e0zs) ([Folien](https://ariadneprojekt.de/media/2025/03/Ariadne_Szen2025_Webinar_Folien_Kernaussagen.pdf)) und im [Ariadne-Webinar zur Energiewirtschaft](https://youtu.be/FcmHBL1MKQA) ([Folien](https://ariadneprojekt.de/media/2025/03/Ariadne_Szen2025_Webinar_Folien_Energiewirtschaft.pdf)) vorgestellt

# PyPSA-DE - High resolution, sector-coupled model of the German Energy System

PyPSA-DE is a sector-coupled energy system model based on the toolbox [PyPSA](https://github.com/PyPSA/pypsa) and the European model [PyPSA-Eur](https://github.com/PyPSA/pypsa-eur). It solves a linear optimization problem to simulate the electricty and hydrogen transmission networks, as well as supply, demand and storage in all sectors of the energy system in Germany and its neighboring countries, as well as Italy and Spain, with high spatial and temporal resolution. PyPSA-DE was developed in the context of the Kopernikus-Projekt [Ariadne](https://ariadneprojekt.de/en/), which studies scenarios of a carbon-neutral German economcy, and plays a decisive role in the [Ariadne Szenarienreport](https://ariadneprojekt.de/publikation/report-szenarien-zur-klimaneutralitat-2045/), as reference model for the [energy and infrastructure sectors](https://ariadneprojekt.de/publikation/report-szenarien-zur-klimaneutralitat-2045/#6-sektorale-perspektive-energiewirtschaft). The results of modeling with PyPSA-DE are also presented in the [Ariadne-Webinar on the core messages of the report](https://youtu.be/UL3KAH7e0zs) ([slides](https://ariadneprojekt.de/media/2025/03/Ariadne_Szen2025_Webinar_Folien_Kernaussagen.pdf)) and in the [Ariadne-Webinar on the energy sector](https://youtu.be/FcmHBL1MKQA) ([slides](https://ariadneprojekt.de/media/2025/03/Ariadne_Szen2025_Webinar_Folien_Energiewirtschaft.pdf)).

This repository contains the entire scientific project, including data sources and code. The philosophy behind this repository is that no intermediary results are included, but all results are computed from raw data and code.

[<img src="https://github.com/PyPSA/pypsa-de/blob/main/doc/img/INFRA_Stromnetzausbau.png?raw=true" width="400"/>](https://github.com/PyPSA/pypsa-de/blob/main/doc/img/INFRA_Stromnetzausbau.png?raw=true)

## Getting ready

You need `conda` or `mamba` to run the analysis. Using conda, you can create an environment from within which you can run the analysis:

```
conda env create -f envs/{os}-pinned.yaml
```

Where `{os}` should be replaced with your operating system, e.g. for linux the command would be:

```
conda env create -f envs/linux-pinned.yaml
```

## Connecting to the Ariadne-Database

### For external users: Use config.public.yaml

The default workflow configured for this repository assumes access to the internal Ariadne2 database. The database will soon be publicly available. Until then, users that do not have the required login details can run the analysis based on the data published during the [first phase of the Ariadne project](https://data.ece.iiasa.ac.at/ariadne/).

This is possible by providing an additional config to the snakemake workflow. For every `snakemake COMMAND` specified in the instructions below, public users should use:

```
snakemake COMMAND --configfile=config/config.public.yaml
```

The additional config file specifies the required database, model, and scenario names for Ariadne1. If public users wish to edit the default scenario specifications, they can do so by changing `scenarios.public.yaml` to `scenarios.manual.yaml`. More details on using scenarios are given below.

### For internal users: Provide login details

The snakemake rule `retrieve_ariadne_database` logs into the interal Ariadne IIASA Database via the [`pyam`](https://pyam-iamc.readthedocs.io/en/stable/tutorials/iiasa.html) package. The credentials for logging into this database have to be stored locally on your machine with `ixmp4`. To do this activate the project environment and run

```
ixmp4 login <username>
```

You will be prompted to enter your `<password>`.

Caveat: These credentials are stored on your machine in plain text.

To switch between internal and public use, the command `ixmp4 logout` may be necessary.

## Run the analysis

Before running any analysis with scenarios, the rule `build_scenarios` must be executed. This will create the file `config/scenarios.automated.yaml` which includes input data and CO2 targets from the IIASA Ariadne database as well as the specifications from the manual scenario file. [This file is specified in the default config.yaml via they key `run:scenarios:manual_file` (by default located at `config/scenarios.manual.yaml`)].

    snakemake build_scenarios -f
or in case of using the public database

    snakemake build_scenarios --configfile=config/config.public.yaml -f

Note that the hierarchy of scenario files is the following: `scenarios.automated.yaml` > (any `explicitly specified --configfiles`) > `config.yaml `> `config.default.yaml `Changes in the file `scenarios.manual.yaml `are only taken into account if the rule `build_scenarios` is executed.

To run the analysis use

    snakemake ariadne_all

This will run all analysis steps to reproduce results. If computational resources on your local machine are limited you may decrease the number of cores by adding, e.g. `-c4` to the call to get only 4 cores. For more option please refer to the [snakemake](https://snakemake.readthedocs.io/en/stable/) documentation.

## Repo structure

* `config`: configuration files
* `ariadne-data`: Germany specific data from the Ariadne project
* `scripts`: contains the Python scripts for the workflow, the Germany specific code needed to run this repo is contained in `scripts/pypsa-de`
* `cutouts`: very large weather data cutouts supplied by atlite library (does not exist initially)
* `data`: place for raw data (does not exist initially)
* `resources`: place for intermediate/processing data for the workflow (does not exist initially)
* `results`: will contain all results (does not exist initially)
* `logs` and `benchmarks`
* The `Snakefile` contains the PyPSA-DE specific snakemake workflow

## Differences to PyPSA-EUR

PyPSA-DE is a softfork of PyPSA-EUR. As such, large parts of the functionality are similar, and the [documentation](https://pypsa-eur.readthedocs.io/en/latest/) of PyPSA-Eur is a good starting point to get acquainted with the model. On top of that, PyPSA-DE adds several data sources and workflow steps that improve the representation of the German Energy System. Below is a non-conclusive list of the most important changes.

- Default resolution of 16 regions in Germany and 13 region for neighboring countries
- 10 pre-defined scenarios (1 Current Policies, 3 Net-Zero Scenarios (Balanced, Focus H2, Focus Electricity), 2 Demand Variations based on the Balanced Scenario, 4 Demand Variations Based on the Current Policies Scenario)
- Specific cost assumption for Germany:

  - Gas, Oil, Coal prices
  - electrolysis and heat-pump costs
  - Infrastructure costs [according to the Netzentwicklungsplan](https://github.com/PyPSA/pypsa-ariadne/pull/193) 2021 and 2023
  - option for pessimstic, mean and optimistic cost development
- Transport and Industry demands as well as heating stock imported from the sectoral models in the Ariadne consortium ([Aladin](https://ariadneprojekt.de/modell-dokumentation-aladin/), [REMOD](https://ariadneprojekt.de/modell-dokumentation-remod/), [FORECAST](https://ariadneprojekt.de/modell-dokumentation-forecast/) and [REMIND](https://ariadneprojekt.de/modell-dokumentation-remind/))
- More detailed data on CHPs in Germany
- The model has been validated against 2020 electricity data for Germany
- National CO2-Targets according to the Klimaschutzgesetz
- Additional constraints that limit maximum capacity of specific technologies
- Import constraints on Efuels, hydrogen and electricity
- Renewable build out according to the Wind-an-Land, Wind-auf-See and Solarstrategie laws
- A comprehensive reporting  module that exports Capacity Expansion, Primary/Secondary/Final Energy, CO2 Emissions per Sector, Trade, Investments, and more.
- Plotting functionality to compare different scenarios
- Electricity Network development until 2030 (and for AC beyond) according to the NEP23
- Offshore development until 2030 according to the Offshore NEP23
- Hydrogen network development until 2028 according to the Wasserstoffkernnetz. PCI / IPCEI projects for later years are included as well.
- `costs:horizon` - specify if technology costs are expected to follow an `optimistic, mean` or `pessimistic` trajectory

## New Config Options

- `iiasa_database` - interaction with IIASA database. Specify a database, and `leitmodelle` for demand and co2 emissions data in specific sectors
- `wasserstoff_kernnetz` - configure which parts of the Wasserstoff Kernnetz should be included in the model
- `new_decentral_fossil_boiler_ban` - specify in which country and which years to ban fossil boilers
- `coal_generation_ban` - specify in which country and which years to ban electricity generation from coal
- `nuclear_generation_ban` - specify in which country and which years to ban electricity generation from nuclear
- `first_technology_occurrence` - specify the year form which on specific technologies are available
- `solving:constraints` - specify PyPSA-DE specific limits, e.g. on capacity, trade and generation
- `co2_budget_DE_source` specify the carbon trajectory for Germany: Following the projections of the Umweltbundestamt (`UBA`) or targeting net zero with the Klimaschutzgesetz(`KSG`)
- `costs:NEP` and `costs:transmission` - specify which year of the Netzentwicklungsplan should be used as basis for the transmission line costs (`2021,2023`) and if new HVDC links should be built with `overhead` or `underground` cables

## License

The code in this repo is MIT licensed, see `./LICENSE.md`.
