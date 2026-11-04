# detonationFoam-Flash: a detonation solver with Lagrangian droplet spray

`detonationFoam-Flash` is a fork of [detonationFoam](https://github.com/JieSun-pku/detonationFoam) (V2.0, based on OpenFOAM 8).
The three solver types of detonationFoam_V2.0 (`solverTypeEuler`, `solverTypeNS_Sutherland`,
`solverTypeNS_mixtureAverage`) are merged into a **single executable** whose functionality is
switched entirely through a dictionary, and an optional **Lagrangian droplet spray cloud**
(OpenFOAM `basicSprayCloud`) is added with two-way coupling to the density-based gas phase.

## What's new?

1. **Unified solver.** The three per-type solver loop trees are replaced by one solver,
   `detonationFoam-Flash`, with a single time loop. Everything is switched through
   `constant/solverTypeProperties` of the case:

   ```C++
   solverType        Euler;        // Euler | NS_Sutherland | NS_mixtureAverage
   viscous           false;        // optional explicit override: on | off
   transportModel    none;         // none | Sutherland | mixtureAverage
   droplets          off;          // Lagrangian droplet cloud on | off
   cloudName         sprayCloud;   // cloud name (constant/<cloudName>Properties)
   SW_track          on;           // leading shock position tracking
   SW_position_limit 1.0;          // normal exit once the shock reaches this x [m]
   ```

   The equation treatment of each mode is identical to the original V2.0 solver-type loops.

2. **Lagrangian droplet spray.** With `droplets on`, a `basicSprayCloud` (injection,
   breakup, evaporation and wall interaction submodels, as in OpenFOAM-8 sprayFoam) is
   coupled two-way to the gas phase through the standard source terms:

   | Gas equation | Source | Physics |
   |--------------|--------|---------|
   | continuity   | `Srho(rho)` | evaporation mass |
   | momentum     | `SU(U)` | vapor-carrying momentum + semi-implicit drag |
   | species      | `SYi(i, Yi)` | vapor composition mapped to gas species |
   | energy       | `Sh(e)` | convective heat transfer + vapor enthalpy |

   The solver requires `constant/g` and a `liquids { ... }` entry in
   `constant/thermophysicalProperties` when droplets are enabled, and the vaporizing
   species must not be the `inertSpecie`.

3. **Tutorial.** An example case with evaporating H2O droplets in a 1D detonation is
   provided in [tutorials/1D_NH3_O2_droplets](tutorials/1D_NH3_O2_droplets).

4. The legacy `solverTypeEuler`, `solverTypeNS_Sutherland` and `solverTypeNS_mixtureAverage`
   directories were removed from the release; they remain retrievable from the git history
   of this fork and of the upstream repository.


## Compiling

1. Install OpenFOAM version 8.

2. Compile:
   ```
   cd detonationFoam-Flash
   ./Allwmake
   ```
   The executable is `detonationFoam-Flash` (installed to `$FOAM_USER_APPBIN`).

## Credit

- Base solver: [detonationFoam](https://github.com/JieSun-pku/detonationFoam) by Jie Sun
  (JieSun-pku); solver documentation of previous versions is available in the upstream
  repository (version 1.0 / V2.0 branches).
- Built on [OpenFOAM](https://openfoam.org) (version 8, GNU GPL), Foundation version.
- Dynamic load-balancing chemistry by the [DLBFoam](https://github.com/blttkgl/DLBFoam)
  authors (Aalto University).
- Low-dissipation schemes from [ROUNDSchemes](https://github.com/ROUNDschemes/libROUNDSchemes).
- The dynamicMesh2D adaptive refinement library originates from detonationFoam_V2.0.

```
