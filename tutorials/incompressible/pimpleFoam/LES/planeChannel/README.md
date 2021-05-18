<!------------------------------------------------------------------------- -->

# Overview

#### Aim

- Validation of the following synthetic turbulence inflow generators:
  - `turbulentDFSEMInlet`

#### Characteristics

- Newtonian
- Incompressible
- Unsteady



# Methodology

- The test case is based on a direct-numerical simulation by
[(Moser et al., 1999)](#ref2). The test case is a smooth-wall turbulent plane
channel flow wherein an internal flow statistically develops downstream (until
is fully developed) through parallel smooth walls that are two characteristic-
length apart. The flow is statistically stationary and fully-developed.
The friction Reynolds number of the flow is $`\mathrm{Re}_\tau = 395`$.

#### Metrics

- The metrics to quantify the validations are:
  - $`\overline{u}`$ vs $`y/\delta`$
  - $`\overline{u^\prime v^\prime}`$ vs $`y/\delta`$
  - $`\overline{u^\prime u^\prime}`$ vs $`y/\delta`$
  - $`\overline{v^\prime v^\prime}`$ vs $`y/\delta`$
  - $`\overline{w^\prime w^\prime}`$ vs $`y/\delta`$
  - $`x/\delta`$ vs $`\overline{C_f}`$

#### Benchmarks

- The numerical results reported in [(Moser et al., 1999)](#ref2)
were obtained from the official website of the data:
[external link](https://turbulence.oden.utexas.edu/data/MKM/chan395/).
  - [Means](https://turbulence.oden.utexas.edu/data/MKM/chan395/profiles/chan395.means)
  - [Reynolds stresses](https://turbulence.oden.utexas.edu/data/MKM/chan395/profiles/chan395.reystress)
  - [Kinetic energy balance](https://turbulence.oden.utexas.edu/data/MKM/chan395/balances/chan395.kbal)
  - [Probability density functions](https://turbulence.oden.utexas.edu/data/MKM/chan395/pdfs/)
  - [Spectra](https://turbulence.oden.utexas.edu/data/MKM/chan395/spectra/)
  - [Correlations](https://turbulence.oden.utexas.edu/data/MKM/chan395/correlations/)
  - [Other profiles](https://turbulence.oden.utexas.edu/data/MKM/chan395/profiles/)

#### Definitions

Various useful expressions are defined as follows:

```math
u_\tau = \sqrt{\frac{\tau_w}{\rho}}
```

```math
\mathrm{Re}_\tau = \frac{u_\tau \, \delta}{\nu}
```

```math
y^+ = \frac{u_\tau \, y}{\nu}
```

```math
u^+ = \frac{u}{u_\tau}
```

```math
C_f = \frac{ |\tau_w| }{0.5 u_b^2}
```

where

| Parameter            | Explanation     |
| ---                  | ---             |
| $`u_\tau`$           | Friction velocity $`\text{[} \mathrm{m}/\mathrm{s} \text{]}`$ |
| $`u_b`$              | Bulk speed $`\text{[} \mathrm{m}/\mathrm{s} \text{]}`$ |
| $`C_f`$              | Skin friction coefficient (based on only x-direction shear stress) |
| $`\mathrm{Re}_\tau`$ | Reynolds number based on friction velocity |
| $`\delta`$           | Characteristic length (i.e. channel half height) $`\text{[} \mathrm{m} \text{]}`$ |
| $`\nu`$              | Fluid kinematic viscosity $`\text{[} \mathrm{m}^2/\mathrm{s} \text{]}`$ |


# Numerics

## Physical domain

- A rectangular-prism channel with infinite spanwise
length and walls two-characteristic length apart.

## Physical modelling

#### Reference quantities

| Metric                | Value     |
| ---                   | ---       |
| $`\mathrm{Re}_\tau`$  | $`395`$   |
| $`\mathrm{Re}_{\tau, \mathrm{DNS}}`$ | $`392.24`$  |
| $`\nu`$               | $`1/392.24 = 0.0025494595145829`$ |
| $`u_\tau`$            | $`1 \text{[} \mathrm{m}/\mathrm{s} \text{]}`$ |
| Bulk speed            | $`| \mathbf{u}_b | \approx 17.55 \text{[} \mathrm{m}/\mathrm{s} \text{]}`$ |

#### Turbulence modelling

- [Large eddy simulation](https://en.wikipedia.org/wiki/Large_eddy_simulation)
with the Smagorinsky sub-filter scale model utilising the van Driest
wall-damping function is used. The sub-filter scale model constants:
  - $`C_k \approx 0.0265463553`$
  - $`C_e = 1.048`$
  - $`C_s = (C_k (C_k/C_e)^{0.5})^{0.5} \approx 0.065`$

## Numerical domain modelling

- A constant-volume rectangular-prism domain.
- Orientation:

| Direction             | Value     |
| ---                   | ---       |
| $`\mathrm{x}`$        | Longitudinal direction (streamwise flow)       |
| $`\mathrm{y}`$        | Wall-normal direction                          |
| $`\mathrm{z}`$        | Spanwise direction (statistically homogeneous) |

- Dimensions:

| Dimension             | Value     |
| ---                   | ---       |
| $`\mathrm{x}`$        | $`60 \text{[} \mathrm{m} \text{]}`$ |
| $`\mathrm{y}`$        | $`2 \text{[} \mathrm{m} \text{]}`$ |
| $`\mathrm{z}`$        | $`\pi \text{[} \mathrm{m} \text{]}`$ |

- Sketch:

N/A

## Numerical domain discretisation

#### Spatial-domain resolution

- Resolution:

| Direction             | Value     |
| ---                   | ---       |
| $`\mathrm{x}`$        | 600 [cells]  |
| $`\mathrm{y}`$        | 64 [cells]   |
| $`\mathrm{z}`$        | 70 [cells]   |

- Resolution distribution:
  - Uniform in $`(x,z)`$-directions
    - $`\Delta_x^+ = (\Delta_x \lvert\mathbf{u}_\tau\rvert)/\nu \approx 39.5`$
    - $`\Delta_z^+ \approx 17.7`$
  - Stretched in $`y`$-direction, progressively finer towards walls
    - Wall-normal cell-to-cell expansion ratio: $`10.7028`$
    - First wall-normal cell-centre height in wall units: $`y^+ \approx 3.8\text{e-}3`$

- Sketch:

N/A

#### Temporal-domain resolution

- Resolution: Constant time-step size of $`\Delta_t = 2\text{e-}3 \text{[} \mathrm{s} \text{]}`$.
  - (Estimated) Courant-Friedrichs-Lewy (CFL) number based on
  $`\{\overline{u_x}\}_{y^+ = 392} = 20.133 \text{[} \mathrm{m}/\mathrm{s} \text{]}`$:
  $`\mathrm{CFL} \approx 0.4`$

## Equation discretisation

- Spatial derivatives and variables: refer to `system/fvSchemes`.
- Temporal derivatives and variables: refer to `system/fvSchemes`.

## Numerical boundary/initial conditions

- Refer to `0.orig`.

## Pressure-velocity coupling algorithm

- Refer to `system/fvSolution`.

## Linear solvers

- Refer to `system/fvSolution`.

## Initialisation and sampling

- Initialisation:
  - A domain pass-through (based on
  $`\{\overline{u_x}\}_{y^+ = 392} = 20.133 \text{[} \mathrm{m}/\mathrm{s} \text{]}`$):
  $`T \approx 3 \text{[} \mathrm{s} \text{]}`$.
  - Number of pass-throughs for initialisation: $`20 \text{[-]}`$.
  - Initialisation duration: $`60 \text{[} \mathrm{s} \text{]}`$.
- Sampling:
  - Number of pass-throughs for sampling: $`40 \text{[-]}`$
  - Sampling duration: $`120 \text{[} \mathrm{s} \text{]}`$.



# Execution and configuration

- Execution:
  - For single-processor or multi-processor (default) runs:
  ```
  ./Allrun
  ./plot
  ```
  - For inspection of results and plots,
  see `results` and `plots` directories after case runs.
  - To clean/reset the cases:
  ```
  ./Allclean
  ```
- Configuration:
  - Input settings of `Allrun` and `plot` can be
  configured from their `settings` section.
  - New setups can be added under `setups.orig` directory.



# Results

- Last test:
  - api      = 2102
  - patch    = 210414
  - HEAD     = 739c1c1d61
  - compiler = Clang (system) clang version 9.0.1
  - mpi      = SYSTEMOPENMPI mpirun (Open MPI) 1.10.7.0.5e373bf1fd
  - OS       = openSUSE Leap 15.1
  - opts     = linux64ClangDPInt32Opt

#### Illustrations

N/A



# Discussion

- For DFSEM, "although the above results show the benefit of taking large
values for `d`, the parameter is also significant from a computational
standpoint...somewhat logically, that a sensible compromise is to take an eddy
density of `d` close to unity, which appears to give a good representation of
the velocity PDF with a relatively fast simulation time."
[(Poletto et al., (2013))](#ref1).



# References

Warning [^1].

[1]<a name="ref1"></a>:
Poletto, R., Craft, T., & Revell, A. (2013).
A new divergence free synthetic eddy method for
the reproduction of inlet flow conditions for LES.
Flow, turbulence and combustion, 91(3), 519-539.
DOI:[10.1007/s10494-013-9488-2](https://doi.org/10.1007/s10494-013-9488-2)

[2]<a name="ref2"></a>:
Moser, R. D., Kim, J., & Mansour, N. N. (1999).
Direct numerical simulation of turbulent
channel flow up to $`\mathrm{Re}_\tau = 590`$.
Physics of fluids, 11(4), 943-945.
DOI:[10.1063/1.869966](https://doi.org/10.1063/1.869966)


[^1]: This document must not be used for academic publications and purposes.

<!------------------------------------------------------------------------- -->
