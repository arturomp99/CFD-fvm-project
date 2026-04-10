# CFD FVM Project

Finite Volume Method (FVM) solver for the 1D Euler equations (inviscid compressible flow).  
Given a mesh and initial conditions, it advances the conserved variables — density, momentum and total energy — in time using an implicit or explicit Euler integrator.

---

# Installation

## Option A — Clone with Git (recommended)

Requires [Git](https://git-scm.com/) installed.

```bash
git clone https://github.com/arturomp99/CFD-fvm-project.git
cd CFD-fvm-project
```

To pull future updates:

```bash
git pull
```

## Option B — Download as ZIP

1. Go to the repository page: https://github.com/arturomp99/CFD-fvm-project
2. Click the green **Code** button → **Download ZIP**
3. Unzip the archive and open the resulting folder in MATLAB

---

# Requirements

- MATLAB R2019b or later (uses `polyshape`, `readlines`, `linsolve`)

---

# How to run

### 1. Provide mesh files

Create a `Meshes/` directory at the project root and populate it with `.dat` files:

| File                     | Content                                                                        |
| ------------------------ | ------------------------------------------------------------------------------ |
| Nodes file               | Tab-separated node coordinates `(x  y)`, one node per line                     |
| Cells file               | Tab-separated node indices forming each cell, one cell per line                |
| Boundary condition files | One file per boundary patch; each file contains the node indices on that patch |

### 2. Set the mesh paths — `constants/FilePaths.m`

Edit `FilePaths.m` to point to your mesh files:

```matlab
NODES = 'Meshes/coords_1D_10nodes.dat';
CELLS = 'Meshes/cells_1D_10nodes_allsquares_shuffled.dat';
BOUNDARY_CONDITIONS = { ...
    'Meshes/bc_left_1D_10nodes.dat', ...
    'Meshes/bc_right_1D_10nodes.dat', ...
};
```

### 3. Configure the simulation — `constants/Config.m`

`Config.m` is the main configuration file. It exposes the following constants:

#### `INITIAL_CONDITIONS`

A function handle `@(cells_centroid_x)` that returns the initial state vector  
`[density (N×1); momentum (N×1); total_energy (N×1)]` given the cell centroid x-coordinates.

Two built-in options are provided in `initial_conditions/`:

- **`sod`** — Sod shock tube: assigns different primitive variables on each side of a discontinuity.
  ```matlab
  INITIAL_CONDITIONS = @(pos) sod( ...
      struct('left', 1.,    'right', 0.1),   ... % pressure [Pa]
      struct('left', 1.,    'right', 0.125), ... % density [kg/m^3]
      struct('left', 0.,    'right', 0.0),   ... % velocity [m/s]
      0.5, pos ...                               % shock position [m]
  );
  ```
- **`uniform`** — Uniform initial state across all cells.
  ```matlab
  INITIAL_CONDITIONS = @(pos) uniform(101325., 1.225, 0., pos);
  ```

#### `CONVECTIVE_FLUX_INTERPOLATOR`

A function handle `@(state, cells)` that selects the numerical scheme used to compute convective fluxes at cell faces.

Two schemes are available in `convective_flux/interpolators/`:

- **`upwind_interpolator`** _(default)_ — First-order upwind scheme. Stable but diffusive.
  ```matlab
  CONVECTIVE_FLUX_INTERPOLATOR = @(state, cells) upwind_interpolator(state, cells);
  ```
- **`linear_interpolator`** — Central (second-order) scheme. Less diffusive but may be unstable without additional stabilisation.
  ```matlab
  CONVECTIVE_FLUX_INTERPOLATOR = @(state, cells) linear_interpolator(state, cells);
  ```

### 4. Run

Open MATLAB, set the working directory to the project root, and run:

```matlab
main
```

Results are displayed as three space-time contour plots (density, momentum, total energy).

---

# Project structure

```
main.m                      Entry point
constants/
  Config.m                  Simulation configuration (initial conditions, scheme)
  FilePaths.m               Paths to mesh .dat files
  Air.m                     Physical constants for air (γ, Cp, Cv, R)
mesh_processing/            Mesh loading and geometry computation
initial_conditions/         Built-in initial condition functions (sod, uniform)
convective_flux/            Convective flux assembly
  interpolators/            Upwind and linear interpolation schemes
propagators/                Time integrators (forward Euler, backward Euler)
timestep_control/           Timestep calculators (constant, CFL-based)
stopping_criteria/          Simulation stopping conditions
results_manager/            Results sampling and storage
utils/                      Shared utilities (thermodynamics, mesh helpers)
visualizer.m                Post-processing: space-time contour plots
```
