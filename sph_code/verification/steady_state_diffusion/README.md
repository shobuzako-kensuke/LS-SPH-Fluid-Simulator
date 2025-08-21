# Steady-State Diffusion

## 🚩Overview

### Problem Statement
This code solves for a steady-state solution of the following two-dimensional diffusion equation:
```math
\begin{align}
\frac{\partial f}{\partial t} = \frac{\partial^{2} f}{\partial x^{2}} + \frac{\partial^{2} f}{\partial y^{2}} + 2\pi^{2} \sin(\pi x) \cos(\pi y)~~,~~x,y \in [0,1]~.
\end{align}
```
The steady-state solution is analytically given by
```math
\begin{align}
f(x,y) = \sin(\pi x) \cos(\pi y)~.
\end{align}
```
The boundary conditions are described below:  
| Wall |Dirichlet | Neumann |
|:---|:---:|:---:|
|Bottom wall ($y=0$)|$f=\sin(\pi x)$|$f_{y}=0$|
|Top wall ($y=1$)|$f=-\sin(\pi x)$|$f_{y}=0$|
|Left wall ($x=0$)|$f=0$|$f_{x}=\pi \cos(\pi y)$|
|Right wall ($x=1$)|$f=0$|$f_{x}=-\pi \cos(\pi y)$|

> [!NOTE]
> The values of $f$ at the four corners are zero.

|Granular Column Collapse (2D) | Initial setting (an example) |
|:---:|:---:|
|<img src="https://github.com/user-attachments/assets/cb0b81a5-61bb-4860-b0f8-f94014b2cc68" alt="granular_column_collapse" width=300>|<img src="https://github.com/user-attachments/assets/c3709de8-9f0e-4d21-8b0a-f8f83a1f303c" alt="initial_setting_granular_column_collapse_2D" width=300>|


### Numerical Models

This code supports the following spatial discretization models: classical SPH, LS-SPH, and ELS-SPH.
Please see [Implemented SPH Models](#implemented-sph-models) for more details.
Eq. (1) is solved by second-order Runge-Kutta method with the initial condition $f(x,y)=0$.
The convergence threshold is set by `threshold`.

### Calculation Options

In Makefile, you can select the following options.
| Options Name| Description | Available |
|:---|:---|:---|
|`SPH_model` | discretization model | `classical_SPH_Laplacian`, `LS_SPH_2ND`, `LS_SPH_3RD`|
|`kernel_function`| kernel function | `CUBIC_SPLINE`, `QUINTIC_SPLINE`, `WENDLAND_C2`, `WENDLAND_C4`, `WENDLAND_C6`|
|`variable_density`| variable mass density or not | `VARIABLE_DENSITY`, `CONSTANT_DENSITY`|
|`virtual_marker`| virtual markers for wall particles |`VM_ON`, `VM_OFF`|
|`wall_accuracy`| discretization accuracy for virtual markers |`WL_1ST`, `WL_2ND`, `WL_3RD`, `WL_4TH`|
|`wall_*`| boundary conditions | `DIRICHLET_*`, `NEUMANN_*`|


### Implemented SPH models

<!-- `classical_SPH_0th`  
```math
\begin{align}
f _{i}=\sum_{j\in \Omega _{i}}f _{j}W_{ij}V_{j} + \mathcal{O}\left((\Delta x)^{0} \right)~,
\end{align}
```
where $f _{i} \equiv f(\bm{x} _{i})$, $\Omega _{i}$ is the index set of the particles located within the support of the kernel function centered at $\bm{x} _{i}$, $W _{ij} \equiv W(r _{ij};h)$ is the kernel function ($r _{ij} \equiv |\bm{x} _{i}-\bm{x} _{j}|$ and $h$ is the smoothing length), $V _{j}$ is the volume, and $\Delta x$ is the typical particle distance.


`classical_SPH_1st`  
```math
\begin{align}
\nabla f|_{\bm{x}=\bm{x} _{i}} = \sum_{j\in \Omega _{i}}f _{j}\nabla _{i} W_{ij} V_{j} + \mathcal{O}\left((\Delta x)^{-1} \right)~.
\end{align}
```

`classical_SPH_1st_sum`  
```math
\begin{align}
\nabla f|_{\bm{x}=\bm{x} _{i}} = \sum_{j\in \Omega _{i}}\left(f _{j}+f_{i}\right)\nabla _{i} W_{ij} V_{j} + \mathcal{O}\left((\Delta x)^{-1} \right)~.
\end{align}
```

`classical_SPH_1st_dif`  
```math
\begin{align}
\nabla f|_{\bm{x}=\bm{x} _{i}} = \sum_{j\in \Omega _{i}}\left(f _{j}-f_{i}\right) \nabla _{i} W_{ij} V_{j} + \mathcal{O}\left((\Delta x)^{0} \right)~.
\end{align}
``` -->

`classical_SPH_Laplacian`  
```math
\begin{align}
\nabla^{2} f|_{\bm{x}=\bm{x} _{i}} = 2\sum_{j\in \Omega _{i}}\left(f _{j}-f_{i}\right) \frac{\bm{x} _{ji} \cdot \nabla _{i} W_{ij}}{r _{ij}^{2}} V_{j} + \mathcal{O}\left((\Delta x)^{-1} \right)~,
\end{align}
```
where $f _{i} \equiv f(\bm{x} _{i})$, $\Omega _{i}$ is the index set of the particles located within the support of the kernel function centered at $\bm{x} _{i}$, $\bm{x}_{ji} \equiv \bm{x}_{j}-\bm{x}_{i}$, $W _{ij} \equiv W(r _{ij};h)$ is the kernel function ($r _{ij} \equiv |\bm{x} _{i}-\bm{x} _{j}|$ and $h$ is the smoothing length), $V _{j}$ is the volume of particle $j$, and $\Delta x$ is the typical particle distance.

<!-- `corrected_SPH_0th`  
```math
\begin{align}
f _{i}=C _{i}^{-1}\sum_{j\in \Omega _{i}}f _{j}W_{ij}V_{j} + \mathcal{O}\left((\Delta x)^{1} \right)~,
\end{align}
```
where
```math
\begin{align}
C _{i} \equiv \sum_{j\in \Omega _{i}}W_{ij}V_{j}~.
\end{align}
```

`corrected_SPH_1st`  
```math
\begin{align}
\nabla f|_{\bm{x}=\bm{x} _{i}}=\bm{B} _{i}^{-1}\sum_{j\in \Omega _{i}}\left(f _{j} - f _{i}\right) \nabla _{i} W_{ij}V_{j} + \mathcal{O}\left((\Delta x)^{1} \right)~,
\end{align}
```
where
```math
\begin{align}
\bm{B} _{i} \equiv \sum_{j\in \Omega _{i}} (\nabla _{i}W_{ij}) \bm{x} _{ji}^{T} V_{j}~.
\end{align}
``` -->

`LS_SPH_2ND` and `LS_SPH_3RD` consider up to the second- and third-order terms in the 2D Taylor expansion, respectively.
For `LS_SPH_2ND`, 
```math
\begin{align}
\bm{M}\tilde{\bm{d}}=\bm{b}~,
\end{align}
```
where
```math
\begin{align}
\bm{M} &\equiv \sum _{j\in \Omega _{i}} \bm{a} \bm{a}^{T} W_{ij} V_{j}~, \\
\bm{a} &\equiv 
\begin{bmatrix}
\dfrac{x_{ji}}{h} & \dfrac{y_{ji}}{h} & \dfrac{1}{2!} \dfrac{x_{ji}^{2}}{h^{2}} & \dfrac{x_{ji} y_{ji}}{h^{2}} & \dfrac{1}{2!} \dfrac{y_{ji}^{2}}{h^{2}}
\end{bmatrix}^{T}~, \\
\tilde{\bm{d}} &\equiv 
\begin{bmatrix}
h~\widetilde{\dfrac{\partial f}{\partial x}\Bigg|_{\bm{x}_{i}}} &h~\widetilde{\dfrac{\partial f}{\partial y}\Bigg|_{\bm{x}_{i}}} & h^{2}~\widetilde{\dfrac{\partial^{2} f}{\partial x^{2}}\Bigg|_{\bm{x}_{i}}} & h^{2}~\widetilde{\dfrac{\partial^{2} f}{\partial x \partial y}\Bigg|_{\bm{x}_{i}}} & h^{2}~\widetilde{\dfrac{\partial^{2} f}{\partial y^{2}}\Bigg|_{\bm{x}_{i}}}
\end{bmatrix}^{T}~, \\
\bm{b} &\equiv \sum_{j \in \Omega _{i}} \bm{a} f_{j} W_{ij} V_{j}~,
\end{align}
```
where a tilde denotes an approximation, and the discretization error of their derivatives are of second-order and first-order accuracies, respectively.  

<!-- `LS_SPH_type_B` is produced by eliminating the first term in the vectors $\bm{a}$ and  $\tilde{\bm{d}}$. -->

> [!TIP]
> Further details of the LS-SPH model are described in [Shobuzako et al. (2025)](https://www.sciencedirect.com/science/article/pii/S2590037425000585).


## 📑 Fundamental Files

|File name|Description|
|:---|:---|
|`input.f90`|Defines all simulation parameters (see `input.f90` for details).|
|`main.f90`|The main program for the simulation, which calls Fortran subroutines.|
|`main.py`|Used for visualizing simulation results.|
|`Makefile`|Compiles the Fortran files and Sets calculation options|
|`initialize_local.py`|Resets this local directory to its initial, distributed state.|

## 📁 Simulation Results Storage

- **Simulation results** are saved as binary files in the `output` directory, which is created when the simulation begins.
- **Simulation figures** are saved in the `fig` directory, which is created after `main.py` is run.
