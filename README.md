# Numerical Simulation and Data Analysis Pipeline

Python pipeline for numerical simulation, parameter sweeps, data analysis, and visualization.

This project applies scientific computing techniques to a hidden vector dark matter model. The code solves a Boltzmann differential equation, runs simulations over different model parameters, analyzes the numerical outputs, and generates visualizations of the results.

## What this project demonstrates

* Numerical solution of differential equations
* Parameter sweeps and simulation analysis
* Data processing and filtering
* Scientific visualization
* Modular Python programming
* Use of NumPy, SciPy, Matplotlib, Pandas, and multiprocessing

## Repository structure

```text
evolver/      Core code for the model, equations, integrator, and utilities
scripts/      Scripts for running simulations, sweeps, and plots
notebooks/    Exploratory analysis notebooks
figures/      Representative visualizations
data/         Output folder for generated simulation files
```

## Example use

Run a single simulation:

```bash
python scripts/run_single.py
```

Run a parameter sweep:

```bash
python scripts/run_sweep.py
```

Analyze sweep results:

```bash
python scripts/analyze_sweep.py
```

## Example output

![Direct Detection Constraints](figures/direct_detection_1_pretty.png)

## About

Although the case study comes from theoretical physics, the workflow demonstrates transferable skills in data analysis, numerical modeling, simulation, and visualization.

