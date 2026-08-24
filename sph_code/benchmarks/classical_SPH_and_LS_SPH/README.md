# Benchmarks with Classical SPH and Least Squares-SPH (LS-SPH)

## 🚩Implemented Examples

This code implements:  
1. Taylor-Green vortex (TG)
2. Lid-driven cavity flow (CF)
3. Boussinesq convection (BC)

## 📑 Fundamental Files

|File name|Description|
|:---|:---|
|`input.f90`|Defines all simulation parameters (see `input.f90` for details).|
|`main.f90`|The main program for the SPH simulation, which calls Fortran subroutines.|
|`*_main.py`|Used for visualizing simulation results.|
|`Makefile`|Compiles the Fortran files.|
|`initialize_local.py`|Resets this local directory to its initial, distributed state.|

## 📁 Simulation Results Storage

- **Simulation results** are saved as binary files in the `output` directory, which is created when the simulation begins.
- **Simulation figures** are saved in the `fig` directory, which is created after `*_main.py` is run.

## 📗 Manual

Manual for this code is available at the link below.
[https://doi.org/10.1016/j.rinam.2025.100594](https://doi.org/10.1016/j.rinam.2025.100594)  