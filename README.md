# Kopernikus-Projekt Ariadne - Gesamtsystemmodell PyPSA-DE

Dieses Repository enthält das Gesamtsystemmodell PyPSA-DE für das Kopernikus-Projekt Ariadne, basierend auf der Toolbox PyPSA und dem Datensatz PyPSA-Eur. Das Modell bildet Deutschland mit hoher geographischer Auflösung, mit voller Sektorenkopplung und mit Integration in das europäische Energiesystem ab.

This repository contains the entire scientific project, including data sources and code. The philosophy behind this repository is that no intermediary results are included, but all results are computed from raw data and code.

## Getting ready

You need conda or [mamba](https://mamba.readthedocs.io/en/latest/) to run the analysis. Using mamba, you can create an environment from within you can run it:

    mamba env create -f environment.yaml

## For external users: Use config.public.yaml

The default workflow configured for this repository assumes access to the internal Ariadne2 database. Users that do not have the required login details can run the analysis based on the data published during the [first phase of the Ariadne project](https://data.ece.iiasa.ac.at/ariadne/).

This is possible by providing an additional config to the snakemake workflow. For every `snakemake COMMAND` specified in the instructions below, public users should use:

```
snakemake --configfile=config/config.public.yaml COMMAND
```

The additional config file specifies the required database, model, and scenario names for Ariadne1. If public users wish to edit the default scenario specifications, they should change `scenarios.public.yaml` instead of `scenarios.manual.yaml`. More details on using scenarios are given below.

## For internal users: Provide login details

The snakemake rule `retrieve_ariadne_database` logs into the interal Ariadne IIASA Database via the [`pyam`](https://pyam-iamc.readthedocs.io/en/stable/tutorials/iiasa.html) package. The credentials for logging into this database have to be stored locally on your machine with `ixmp4`. To do this activate the project environment and run

```
ixmp4 login <username>
```

You will be prompted to enter your `<password>`.

Caveat: These credentials are stored on your machine in plain text.

To switch between internal and public use, the command `ixmp4 logout` may be necessary.

## Run the analysis

Before running any analysis with scenarios, the rule `build_scenarios` must be executed. This will create the file `config/scenarios.automated.yaml` which includes input data and CO2 targets from the IIASA Ariadne database as well as the specifications from the manual scenario file. [This file is specified in the default config.yaml via they key `run:scenarios:manual_file` (by default located at `config/scenarios.manual.yaml`)].

    snakemake -call build_scenarios -f

Note that the hierarchy of scenario files is the following: `scenarios.automated.yaml` > (any `explicitly specified --configfiles`) > `config.yaml `> `config.default.yaml `Changes in the file `scenarios.manual.yaml `are only taken into account if the rule `build_scenarios` is executed.

For the first run, open config.yaml and set

    enable:
        retrieve: true # set to false once initial data is retrieved
        retrieve_cutout: true # set to false once initial data is retrieved

and then run from main repository

    snakemake -call

This will run all analysis steps to reproduce results.

To generate a PDF of the dependency graph of all steps `build/dag.pdf` run:

    snakemake -c1 --use-conda -f dag

## Repo structure

* `config`: configuration files
* `ariadne-data`: Germany specific data from the Ariadne project
* `workflow`: contains the Snakemake workflow, including the submodule PyPSA-Eur and specific scripts for Germany
* `cutouts`: very large weather data cutouts supplied by atlite library (does not exist initially)
* `data`: place for raw data (does not exist initially)
* `resources`: place for intermediate/processing data for the workflow (does not exist initially)
* `results`: will contain all results (does not exist initially)

## Differences to PyPSA-EUR

- Specific cost assumption for Germany:
  - Gas, Oil, Coal prices
  - electrolysis and heat-pump costs
  - Infrastructure costs according to the Netzentwicklungsplan 23 (NEP23)
  - option for pessimstic, mean and optimistic cost development
- Transport and Industry demands as well as heating stock imported from the sectoral models in the Ariadne consortium
- More detailed data on CHPs in Germany
- Option for building the German Wasserstoffkernnetz
- The model has been validated against 2020 electricity data for Germany
- National CO2-Targets according to the Klimaschutzgesetz
- Additional constraints that limit maximum capacity of specific technologies
- Import constraints
- Renewable build out according to the Wind-an-Land, Wind-auf-See and Solarstrategie laws
- A comprehensive reporting  module that exports Capacity Expansion, Primary/Secondary/Final Energy, CO2 Emissions per Sector, Trade, Investments, ...
- Plotting functionality to compare different scenarios

## License

The code in this repo is MIT licensed, see `./LICENSE.md`.
