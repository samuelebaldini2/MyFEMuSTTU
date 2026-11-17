# Static surface-force range diagnostic

This standalone diagnostic uses the branch's native CutFEM and FEM basis sources to assemble the two Q1 pressure-to-Q2-velocity blocks and the surface-tension force on a normalized spherical drop. It is not the full PETSc application residual. It removes normal velocity components on cube faces, retains both Q1 phase fields with no pressure reclassification, and omits stabilization. `CUTFEM_PROJECT_Q1=1` projects the mesh's high-order geometry nodes to the corner-defined Q1 geometry. Both moment methods use the identical mesh and zero-velocity pressure-jump setup.

`FEMUS_CUTFEM_UNLIMITED=1` raises the evaluation budget to one billion and depth to 48; it does not literally remove limits. The probe rejects moment warnings. The three-normal calculation obtains the volume weights directly with its first normal moment, instead of ex46's current scalar-first cache path. Area is not measured in this mode.

Build with Boost, Eigen, CMake and a C++14 compiler; the solver requires Python, NumPy and SciPy:

```sh
cmake -S unittests/testCutFemWeightLSNative -B build-range -DFEMUS_NATIVE_BOOST=/path/to/boost -DFEMUS_NATIVE_EIGEN=/path/to/eigen -DCMAKE_BUILD_TYPE=Release
cmake --build build-range --target cutfem_surface_range -j 2
FEMUS_CUTFEM_UNLIMITED=1 CUTFEM_PROJECT_Q1=1 CUTFEM_GK=31 build-range/cutfem_surface_range run B applications/MultiPhase/ex46/input/unstructured_cube_hex27.neu 17 1e-12 0 /tmp/range-one
FEMUS_CUTFEM_UNLIMITED=1 CUTFEM_PROJECT_Q1=1 CUTFEM_GK=31 CUTFEM_THREE=1 build-range/cutfem_surface_range run B applications/MultiPhase/ex46/input/unstructured_cube_hex27.neu 10 1e-12 0 /tmp/range-three
OPENBLAS_NUM_THREADS=1 python3 validation/surface-range/solve.py /tmp/range-one --max-relative 1e-12
OPENBLAS_NUM_THREADS=1 python3 validation/surface-range/solve.py /tmp/range-three --max-relative 1e-12
```

The solver assembles the sparse operator, verifies its action on the prescribed Laplace pressure and the constant-pressure gauge, and uses column-scaled LSMR from zero to minimize `||R_sigma + B^T p||`. It recomputes the residual using long-double accumulation and divides by `||R_sigma||`. A seeded random-force negative control checks that arbitrary forces are not declared balanced. The pressure is nonunique due to the gauge and weakly supported DOFs. The output is a numerical distance, not an exact-arithmetic certificate.

The three-moment q=10 check on the included mesh with the commands above returned a relative range distance of `5.776985586474213e-15` (absolute `2.4844336522281788e-15`), using 18,405 free velocity DOFs and 1,068 active pressure DOFs. The random-force negative-control distance was `0.9715449212191232`.

The one-moment q=17 check returned `9.189790727230102e-15` relative (absolute `3.952134724573427e-15`) with the same number of free velocity and active pressure DOFs. Its random-force distance was `0.9715449212188259`. Both results use the enlarged budget and quad moment-to-weight arithmetic; integration warnings must be treated as failures.
