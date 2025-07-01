# FGSIM Application: Influenza Surveillance in Germany

This repository contains the application code for the Fine-Grid Spatial Interaction Matrix (FGSIM) methodology applied to German influenza surveillance data, as described in the paper "Fine-Grid Spatial Interaction Matrices for Surveillance Models with Application to Influenza in Germany" by Manuel Stapper and Sebastian Funk.

## Overview

FGSIM is a novel spatial modeling approach for infectious disease surveillance that defines transmission risk at the individual level and aggregates it to the district level. Unlike traditional methods that rely on administrative boundaries or simple proximity measures, FGSIM uses high-resolution population density data to model contact patterns between individuals across different geographical areas.

### Key Features

- **Individual-level modeling**: Defines transmission risk between randomly selected individuals from different districts
- **Multiple distance measures**: Supports 5 different distance metrics including beeline distance, travel time, gravity model, circle distance, and radiation model
- **High-resolution risk mapping**: Generates risk maps at 100m resolution, independent of administrative boundaries
- **Endemic-epidemic framework**: Built on the well-established surveillance modeling framework

## Methodology

FGSIM extends traditional spatial epidemiological models by:

1. **Sampling individuals** from districts according to population density (using WorldPop 100m grid data)
2. **Computing distances** between randomly selected individual pairs using various metrics
3. **Fitting distributions** to distance samples (mixture of truncated log-normal distributions)
4. **Calculating transmission weights** using power-law decay functions with closed-form expectations
5. **Aggregating** individual-level risks to district-level transmission matrices

### Distance Measures

The implementation includes five distance measures:

- **Beeline Distance**: Geodesic distance using Haversine formula
- **Travel Time**: Car travel time using OpenStreetMap routing data
- **Gravity Model**: Incorporates population densities at origin and destination
- **Circle Distance**: Accounts for population within circular catchment areas
- **Radiation Model**: Based on intervening opportunities theory

## Data

The application uses:
- **Influenza surveillance data**: Weekly case counts from Robert Koch Institute (2001-2020)
- **Population data**: WorldPop 100m resolution grid data for Germany
- **Administrative boundaries**: 401 German districts grouped into 4 regions (North, East, South, West)
- **Geographic data**: OpenStreetMap for routing calculations

### Data Exclusions
- H1N1 pandemic period (weeks 17/2009 - 16/2010) excluded due to different transmission mechanisms
- COVID-19 overlap period handled as endemic-level influenza


## Dependencies

- **FGSIM Toolbox**: General tool, not applciation-specific at [GitHub Repo](https://github.com/ManuelStapper/FGSIM)
- **surveillance**: Endemic-epidemic modeling framework
- **geosphere**: Geodesic distance calculations  
- **osrm**: OpenStreetMap routing interface
- **WorldPop**: High-resolution population data
- **mixtools**: Expectation-maximization algorithms
- **sp/sf**: Spatial data handling

## File overview

- Beeline.R, Gravity.R, Radiation.R, RouteTime.R: Computing the distance measures
- Radiation.cpp: C++ helper functions for circle distance and radiation distance without correction
- Radiation.cpp: C++ helper functions for circle distance and radiation distance with correction
- Fitting.R: Functions to fit mixture of truncated LogNormal distribution
- PrepInput.R: Translating parameter estimates to a better format
- IDlookup.cdv: Mapping shapefile to district IDs
- Dates.csv: Overview of dates included in study
- Application.R: Main application file


## License

This project is licensed under the MIT License

## Acknowledgments

- **Funding**: Wellcome Trust (210758/Z/18/Z)
- **Data**: Robert Koch Institute, WorldPop, OpenStreetMap
- **Institution**: Department of Infectious Disease Epidemiology and Dynamics, London School of Hygiene and Tropical Medicine

## Contact

- **Manuel Stapper**: [manuel.stapper@lshtm.ac.uk]
- **Sebastian Funk**: [sebastian.funk@lshtm.ac.uk]
