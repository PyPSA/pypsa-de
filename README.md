
# PyPSA-AT: A Sector-Coupled Open Optimisation Model of the Austrian Energy System

**PyPSA-AT** is an Austrian adaptation of the open European energy system model [PyPSA-Eur](https://github.com/pypsa/pypsa-eur). It provides a detailed model of the Austrian energy system, including electricity, heating, and transport sectors.

The project builds upon the methodologies developed in [PyPSA-DE](https://github.com/pypsa/pypsa-de) - the German adaptation of PyPSA-Eur - while incorporating Austrian-specific network topology, energy system characteristics, and regulatory frameworks. It leverages established modeling approaches for electricity system calibration and infrastructure planning, adapted for the Austrian context.

## Features

- High-resolution model of the Austrian transmission system
- Integration of Austrian-specific energy data sources
- Detailed representation of district heating networks
- Consideration of Austrian energy policies and targets
- Enhanced spatial resolution for Austria while maintaining compatibility with neighboring countries

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/AGGM-AG/pypsa-at.git
   cd pypsa-at
   ```

2. Installation using pixi (recommended):
   ```bash
   pixi install
   ```

## Data Requirements

The model requires various data sources that are either downloaded automatically or need to be obtained manually due to license restrictions. Key data sources include:

- Austrian transmission grid data
- District heating network information
- Renewable energy potential maps
- Demand data and load profiles

Refer to the documentation for detailed information about data sources and preprocessing steps.

## Usage

1. Configure your model by adjusting the base scenario in `config/config.at.yaml`
2. Include scenario settings that differ from the base scenario in `config/scenarios.manual.yaml`  
3. Generate the scenarios file picked up by the snakemake workflow:
   ```bash
   snakemake build_scenarios -f --cores 'all'
   ```
   This will populate `config/scenarios.autoamted.yaml`. Do not forget to enable sc

4. Run the model using the default rule `all`:
   ```bash
   snakemake -call all --cores 'all'
   ```
   or simply
   ```bash
   snakemake
   ```
## Documentation

Detailed documentation is available at [docs-at folder](./docs-at).

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License - see the [LICENSE.txt](LICENSE.txt) file for details. 

Parts of the code that originate from [PyPSA-DE](https://github.com/pypsa/pypsa-de) or [PyPSA-Eur](https://github.com/pypsa/pypsa-eur) remain under their original MIT licenses. The copyright and attribution notices from these original projects are preserved in the respective source files.

## Acknowledgments

PyPSA-AT builds upon [PyPSA-Eur](https://github.com/pypsa/pypsa-eur) and [PyPSA-DE](https://github.com/pypsa/pypsa-de), developed by the PyPSA team at TU Berlin and other contributors.

## Citation

If you use PyPSA-AT in your research, please cite it as:

[Add your preferred citation format]