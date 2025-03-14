# Functional Factor Models for Yield Curve Forecasting

## Project Overview
This repository contains the code for my master's thesis *Functional Factor Models for Yield Curve Forecasting*. The study extends the approximate functional factor model by Otto and Salish (2024) to improve yield curve forecasting by incorporating macroeconomic principal component factors and modeling conditional heteroskedasticity in the factor dynamics.

The model is compared to the widely used Dynamic Nelson-Siegel (DNS) model in an empirical application to the U.S. Treasury yields data from [Liu & Wu (2021)](https://doi.org/10.1016/j.jfineco.2021.05.059). Forecasting performance is evaluated for different model specifications and forecast horizons, focusing on point and interval forecasts across the yield curve.

## Installation

- Macro data has to be obtained from [the FRED-MD database](https://www.stlouisfed.org/research/economists/mccracken/fred-databases) and the current .csv including all of 2023 has to be inserted into `/data/` as `FRED-MD_current.csv`.

## Usage

- The script of the whole analysis is contained in `main.R`.

## Future Extensions
- Use other types of yield curve data, such as the novel data set of 
- Improve selection of macro factors by using partial least squares or ML-based index creation.
- Use SSM models estimated by the extended Kalman Filter to incorporate conditional heteroskedasticity into the point forecasts.

