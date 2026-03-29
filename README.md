# ECARM — Energy-Consistent Advanced Regional Model

[![DOI](https://zenodo.org/badge/1195137030.svg)](https://doi.org/10.5281/zenodo.19311595)

ECARM (Energy-Consistent Advanced Regional Model) is a modified dynamical
core integrated into the WRF (Weather Research and Forecasting) modeling
framework (version 4.2). It replaces the default ARW (Advanced Research WRF)
dynamical core in the `dyn_em` directory with an energy-consistent
discretization scheme designed to reduce spurious energy sources and sinks
in numerical simulations.

The model's compilation, execution, and output are fully compatible with
standard WRF v4.2 — the same executables (`wrf.exe`, `real.exe`,
`ideal.exe`), the same output files (`wrfout`, `wrfinput`, `wrfbdy`), and
the same `namelist.input` configuration are used.

## Design Philosophy

The design objective of ECARM's dynamical core discretization scheme is to
reduce spurious energy sources and sinks, thereby minimizing numerical
dissipation and dispersion, and enhancing the model's capability to reproduce
the spatiotemporal evolution of simulation targets.

The discretization design principle involves constructing spatial
discretization schemes constrained by the energy properties of the model's
dynamical core, where energy properties refer to spatial energy transport and
interconversion among different energy forms. Specifically, this is
accomplished by deriving discrete expressions of energy properties based on
discrete grid cells as fundamental units.

## Governing Equations

The dynamical core adopts the **Lamb form of the Navier–Stokes equations**.
The basic prognostic variables are:

- *u*, *v*, *w* — wind components
- *ρ* — density
- *ρT* — density-weighted temperature
- *q* — moisture variables

Terrain-following coordinates employ:

- *m* = map-scale factor in *x*-direction (grid spacing / true spacing)
- *n* = map-scale factor in *y*-direction (grid spacing / true spacing)

To alleviate roundoff errors associated with pressure gradient force
calculations under hydrostatic equilibrium, a **hydrostatic balance
decomposition** is introduced: *ρ* and *p* are decomposed into the sum of
basic-state quantities satisfying hydrostatic equilibrium and perturbation
quantities, with perturbation quantities replacing total quantities as
prognostic variables.

## Key Features

- **Energy-consistent spatial discretization**: Discrete schemes are
  constrained by energy transport and interconversion properties, reducing
  spurious numerical energy sources and sinks
- **Reduced numerical dissipation and dispersion**: Improved fidelity in
  reproducing spatiotemporal evolution of atmospheric phenomena
- **Lamb-form Navier–Stokes equations**: With hydrostatic balance
  decomposition for robust pressure gradient force computation over steep
  terrain
- **Full WRF compatibility**: All WRF executables, output file names, physics
  parameterizations, I/O infrastructure, and WPS preprocessing remain
  unchanged — ECARM runs identically to standard WRF from a user perspective

## Relationship to WRF

This code is based on **WRF version 4.2**
([official release](https://github.com/wrf-model/WRF/releases/tag/v4.2)).
The WRF model is in the public domain and its use is unrestricted.

**What is modified:** the dynamical core source code, primarily in `dyn_em/`
and related interface files.

**What is NOT modified:** executable names (`wrf.exe`, `real.exe`,
`ideal.exe`), output file names (`wrfout`, `wrfinput`, `wrfbdy`, `wrfrst`),
namelist structure, physics parameterizations, WPS preprocessing, and all
other WRF components.

| Directory     | Scope of Changes                                             |
|:--------------|:-------------------------------------------------------------|
| `dyn_em/`     | Replaced dynamical core with ECARM energy-consistent solver  |
| `Registry/`   | Updated variable registrations for ECARM prognostic fields   |
| `share/`      | Interface adaptations for the new dynamical core             |
| All others    | Unchanged from WRF v4.2                                      |

## Building

### Prerequisites

Same as WRF 4.2:

- Fortran and C compilers (Intel, GNU, or PGI)
- NetCDF library (with Fortran bindings)
- MPI library (for parallel runs)

### Compilation

```bash
cd ECARM
./configure          # Select your compiler and parallelism option
./compile em_b_wave 2>&1 | tee compile.log
```

After successful compilation, the standard WRF executables appear in `main/`:

- `wrf.exe` — Main model integration (with ECARM dynamical core)
- `real.exe` — Real-data case preprocessor
- `ideal.exe` — Idealized case preprocessor

### Running

ECARM is run in exactly the same way as standard WRF:

```bash
cd test/em_b_wave
ln -sf ../../main/ideal.exe .
ln -sf ../../main/wrf.exe .
./ideal.exe
mpirun -np 4 ./wrf.exe
```

Output files follow the standard WRF naming convention:
`wrfout_d01_YYYY-MM-DD_HH:MM:SS`.

## Citation

If you use ECARM in your research, please cite:

> Ren, X.-Y.: ECARM: Energy-Consistent Advanced Regional Model
> (Version 1.0.0) [code], Zenodo,
> https://doi.org/10.5281/zenodo.19311595, 2026.

Also cite the underlying WRF framework:

> Skamarock, W. C., Klemp, J. B., Dudhia, J., Gill, D. O., Liu, Z.,
> Berner, J., Wang, W., Powers, J. G., Duda, M. G., Barker, D. M., and
> Huang, X.-Y.: A Description of the Advanced Research WRF Model
> Version 4, NCAR Technical Note NCAR/TN-556+STR,
> https://doi.org/10.5065/1dfh-6p97, 2019.

## License

This code is released as public domain, consistent with the WRF model
license. See [LICENSE](LICENSE) for details.

## Contact

- Ren X.-Y. — renxingyu22@mails.ucas.ac.cn
- Institute of Atmospheric Physics， Chinese Academy of Sciences
