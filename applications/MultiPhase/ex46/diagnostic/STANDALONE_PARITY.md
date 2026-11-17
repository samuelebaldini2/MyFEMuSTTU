# ex46 standalone-parity diagnostic controls

This experimental branch starts from `125cfb95d6afaab20340344a3c84df0fb038912a` on `diagnostic/surface-pressure-range`. It is intended to isolate discrepancies in the surface-pressure force balance without altering the standard ex46 algorithm when the diagnostic is disabled.

## Usage

In `run_ex46.sh`, set:

```bash
SURFACE_PRESSURE_DIAGNOSTIC=1
SURFACE_REFERENCE_PROFILE=1
SURFACE_WEIGHTS=3   # 3 for normal moments (q=10), 1 for scalar (q=17)
NPROC=1             # first validation; then compare with NPROC=2,8
```

The standalone parity profile sets the moment order to 10 for three weights, or 17 for one weight; tolerance 1e-12; quad-precision moment-to-weight fitting but double output weights; all four outer GK rules to 15/31; a trilinear Q1 projection of the loaded HEX27 geometry before refinement; three-normal-moment first call order; extended endpoint budget; strict warning rejection; unconstrained pressure-space postprocessing; and diagnostic-only termination after the initial NS solve.

These defaults are individually overridable via `SURFACE_MOMENT_ORDER`, `SURFACE_MOMENT_TOL`, `SURFACE_FIT_PRECISION`, `SURFACE_GK`, `SURFACE_PROJECT_Q1`, `SURFACE_NORMAL_FIRST`, `SURFACE_UNCONSTRAINED`, `SURFACE_STRICT_WARNINGS`, `SURFACE_EXTENDED_BUDGET`, and `SURFACE_PRESSURE_ONLY`.

`SURFACE_REFERENCE_PROFILE=1` requires the diagnostic flag to be explicitly enabled. The comparison controls are disabled by default.

## What is and is not matched

- Matching: input HEX27 mesh, normalized sphere shape (up to the application length scale), moment order, GK rules, 3D endpoint budget, moment tolerances, double output/quad fitting, optional Q1 projection, positive-volume integration call order, pressure columns retained in least squares, and surface/pressure balance definitions.
- Deliberately unchanged: application physical scale and surface-tension coefficient; existing Navier–Stokes matrix, ghost penalties and pressure constraints; physical solver formulation.
- The pressure columns remain constrained in the **actual production NS solve**. Only the postprocessed range diagnostic can disregard their fixed status. This compares a standalone-like pressure space without altering FEMuS.
- The Q1 geometry control changes a diagnostic run's mesh coordinates before refinement. Never use this to characterize the original stored HEX27 geometry; repeat with `SURFACE_PROJECT_Q1=0`.
- Normal-first mode computes the three normal moments first, preserves their positive phase weights, and then computes the tangential scalar weights used by other production terms.
- `FEMUS_CUTFEM_UNLIMITED=1` is an increased finite budget (not literally unlimited). Global warnings across ranks are reported; strict mode rejects warned snapshots.
- `SURFACE_PRESSURE_ONLY=1` still executes the first NS solve but skips all subsequent time iterations. It is not a solver-free assembly mode.

## Output / interpretation

The postprocessor reports `prescribed_relative` (denominator surface norm), `prescribed_relative_maxnorm` (standalone denominator max(surface,pressure)), and `range_relative` (surface norm); the pressure-action check independently verifies matrix extraction. It also enforces an optimality tolerance on the scaled normal residual. Output metadata contains the actual moment order/tolerance, GK rule orders, fitting precision selection, projection flag, pressure-space choice, cut-element count and warning count.

Reproducibility should be established through a one-rank run, then a two-rank run with the same configuration; inspect matching global matrix/force entries before interpreting scalar balance differences. No production MPI/PETSc build or runtime comparison was possible in the editing environment; these are not claimed as validated numerical results.

## Known risks / remaining differences

1. Projection averages corner-defined Q1 positions across owned HEX27 elements using distributed vectors; this is intended to be rank invariant, but still requires an MPI run and an audit of the shared-node consistency.
2. The FEMuS application retains its physical scaling (cube side 0.01, sigma 25); only dimensionless/relative residuals are directly comparable to the standalone normalized-cube calculation.
3. Long-double moment integration still returns moments in long double; selecting quad improves only the moment-to-weight fitting arithmetic, as in the controlled standalone double-output comparison.
4. Standalone and FEMuS use different element assembly implementations and accumulation orders; machine-identical results are not guaranteed.
5. The Python unconstrained-pressure option skips FEMuS fixed pressure columns in the range calculation; it does not change velocity Dirichlet row filtering.
6. The mesh must be HEX27 for the Q1 projection mode. Other geometries and axisymmetric/symmetric configurations are intentionally rejected for the parity control.
7. Very high moment orders plus quad fitting can require substantial memory, especially with multiple MPI ranks.

## Suggested controlled experiment

Run three-weight `q=10` at 1e-12 tolerance, GK 15/31, quad fitting, Q1-projected, unconstrained pressure. Verify zero warnings and compare `range_relative` to the standalone 5.776985586474213e-15. Next turn **one knob at a time**: first Q1 projection, then quad fitting, then call order, then GK, then tolerance and pressure constraints. Do not infer the cause of a numerical floor without these paired runs.

The original source comparison remains in the user-supplied `standalone_vs_ex46_differences.md`; this README describes the implemented experimental controls.
