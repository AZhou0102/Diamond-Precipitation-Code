# Planetary Modeling

This repository contains Python code to model the structure and properties of planets, using different equations of state (EoS) and planetary compositions.

## Key Components

**`looped_solver.ipynb`**/**`looped_solver.py`**: This module iteratively calculates density and pressure profiles inside a planet given a list of radii and a dictionary describing its layered composition. See the linked notebook below for a walkthrough of the code and its usage.
- https://github.com/AZhou0102/Diamond-Precipitation-Code/blob/main/looped_solver.ipynb

**`solve_adams_williamson.py`**: Returns a list of gravities and a list of pressures corresponding to each radius within a planet given a list of radii and densities at each radius.
- https://github.com/AZhou0102/Diamond-Precipitation-Code/blob/main/solve_adams_williamson.py

**`EoS_Bits.py`**: Contains functions and constants related to different equations of state (EoS) for modeling material properties under various pressures and densities.
- https://github.com/AZhou0102/Diamond-Precipitation-Code/blob/main/EoS_Bits.py

**`planetary_dictionary.py`**: Generates a dictionary representing a planet's internal composition, including core, mantle, and other layers. It assigns densities to each region based on given parameters like the planet's radius and material composition.
- https://github.com/AZhou0102/Diamond-Precipitation-Code/blob/main/planetary_dictionary.py

**`monte_carlo_planets.ipynb`**: Monte Carlo simulation of hypothetical exoplanets with randomized composition. Evaluates hypothetical planets for diamond precipitation candidacy. Outputs plots and csv dataframes.
- https://github.com/AZhou0102/Diamond-Precipitation-Code/blob/main/monte_carlo_planets.ipynb

**`prem.ipynb`**: Preliminary Reference Earth Model
- https://github.com/AZhou0102/Diamond-Precipitation-Code/blob/main/prem.ipynb

## Installation

This project requires the following Python packages, in addition to the custom modules available in this repository. Note that other modules may be needed for some of the plotting notebooks and scripts that aren't listed above in the "Key Components" section:
- `numpy`
- `matplotlib`
- `math (part of Python's standard library)`
- `random (part of Python's standard library)`
- `csv (part of Python's standard library)`
- `os (part of Python's standard library)`
