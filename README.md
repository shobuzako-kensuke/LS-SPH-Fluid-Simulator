# LS-SPH-Fluid-Simulator

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15709255.svg)](https://doi.org/10.5281/zenodo.15709255)

[日本語版はこちら](./README_ja.md)  

This open-source code implements various simulations based on the **Least Squares Smoothed Particle Hydrodynamics (LS-SPH) method** [[1](#ref1)] and its advanced version, the **Enhanced and Explicit LS-SPH (ELS-SPH) method**.  

This code includes the following examples:  

1. Fluid benchmark problems
   - Taylor-Green vortex
   - Lid-driven cavity flow
   - Boussinesq convection
   - ~~Oscillating drop~~
   - ~~Dam break~~
2. Applications to geophysical and engineering problems (currently under development)

|Taylor-Green vortex | Lid-driven cavity flow | Boussinesq convection |
|:---:|:---:|:---:|
|<img src="https://github.com/user-attachments/assets/cb0b81a5-61bb-4860-b0f8-f94014b2cc68" alt="granular_column_collapse" width=300>|<img src="https://github.com/user-attachments/assets/fd3fed95-f3f8-4f3f-941c-7b6204d75d1b" alt="cylinder_lift" width=300>|<img src="https://github.com/user-attachments/assets/fd3fed95-f3f8-4f3f-941c-7b6204d75d1b" alt="cylinder_lift" width=300>|


## ⚙️ Requirements

This code is written in **Fortran** (for calculations) and **Python** (for visualization).

| Category | Requirement | Notes |
|:---|:---|:---|
|Operating System |Unix-like environment | Tested on Windows Subsystem for Linux (WSL2)|
|Compiler | Intel Fortran | Tested with `ifx` |
|Build Tool | `make` | For building Fortran files|
|Visualization | `Python` | Tested with `Python 3.12.0` and requires basic libraries (e.g., `matplotlib`)|
|Movie Generation| `ffmpeg` | Required for Python to generate movies|

> [!TIP]
> If `ffmpeg` is not installed, run `sudo apt install ffmpeg` .


## 🖥️ Usage

1. Navigate to the simulation directory (e.g., `sph_code/benchmarks/ELS_SPH/Taylor_Green_vortex/source_code`)
2. Run `make` to build the Fortran programs.
3. Run `./start_calculation` to start the calculation.
4. After the simulation finishes, run `python main.py` to generate figures.

> [!NOTE]
> Each problem setting is described in the **README.md** in its directory.


## 🧑‍💻 Citation

Please cite the following ***two*** references:

<a id="ref1">[1]</a>  
Shobuzako, K., Yoshida, S., Kawada, Y., Nakashima, R., Fujioka, S., & Asai, M. (2025).  
A generalized smoothed particle hydrodynamics method based on the moving least squares method and its discretization error estimation.  
*Results in Applied Mathematics*, 26, 100594. [https://doi.org/10.1016/j.rinam.2025.100594](https://doi.org/10.1016/j.rinam.2025.100594)

[2]  
Shobuzako, K. (2025). *LS-SPH-Fluid-Simulator* (Version 1.1.0) [Computer software]. Zenodo.  
[https://doi.org/10.5281/zenodo.15709255](https://doi.org/10.5281/zenodo.15709255)


## 🤝 Contributing
If you would like to improve the code, report a bug, or add a new feature, feel free to submit a pull request.


## 🪪 License

This project is licensed under [the MIT License](./LICENSE) .
