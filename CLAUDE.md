# KEnRef — project conventions

C++ (Eigen) port of the R `ke` package for kinetic ensemble refinement. Core numerical code lives
in `src/core/KEnRef.cpp` with declarations in `include_core/core/` (`KEnRef.h`, `IoUtils.h`,
`NamedMatrix.h`). The GROMACS/PLUMED interface is built separately (`BUILD_KENREF_GMX`).

## Project layout (three projects)
KEnRef is consumed by two MD engines, so changes to the core's public API ripple into both:
- **kenref_core** — the numerical library (`src/core/`, `include_core/core/`). Holds each energy model's
  compute/backprop: PLATEAUS (`coord_array_to_g_energy`), SIGMA (`coord_array_to_sigma_energy`), RELAX
  (`coord_array_to_relax_energy`). Each model is one self-registering `EnergyModel` subclass
  (`SigmaModel`/`PlateausModel`/`RelaxModel`) in the `ModelRegistry`.
- **kenref_gmx** — GROMACS integration: `KEnRefForceProvider` (`src/gmxinterface/`, `BUILD_KENREF_GMX`,
  build dir `cmake-build-*-gmx-2025.4`).
- **plumed_kenref** — the PLUMED fork at **`/home/amr/git/plumed2`** (separate git repo) with `KEnRefBias`.

The model-abstraction restructure is DONE: all three consumers (the offline tools, the GROMACS force
provider, and PLUMED `KEnRefBias`) select the model by NAME via `ModelRegistry` and run it through the
shared `KEnRefDriver`/`buildModelIndexing` — no per-model enum or switch. **SIGMA + PLATEAUS + RELAX are
all wired** in every consumer; adding a model = one `EnergyModel` subclass + one CMake list entry, no
consumer edit.

## Build & test
- Build dir: `cmake-build-debug` (ninja). Core lib target `kenref_core`; test exe target
  `Google_tests_kenref_core_exe` (gtests live in `google_tests/`, fixtures in `res/google_tests/`).
- Build the tests: `ninja -C cmake-build-debug Google_tests_kenref_core_exe`
- Run one test: `cd cmake-build-debug/google_tests && ./Google_tests_kenref_core_exe --gtest_filter="..."`
- The `TestCoordArrayToSigmaEnergyFD` finite-difference test is slow (~200s); exclude it with
  `--gtest_filter="-KEnRefTestSuite.TestCoordArrayToSigmaEnergyFD"` for quick runs, but run it before
  claiming a numerical kernel is correct.
- Scalar type is `KEnRef_Real_t`, set by the `KENREF_DOUBLE` compile option (double by default).
- SIMD width comes from the `ACCEL` cmake option's `-march` (default `AVX2_256`; also `AVX_512`,
  `AVX_256`).

## Numerical correctness: validate against R
Every numerical function is a faithful port of an R `ke` reference and must be validated against
R-generated ground truth (energy/gradient/intermediate values), not just self-consistency. The R
package is at `/home/amr/PycharmProjects/ke` (load with `pkgload::load_all(...)`, not `library(ke)`).
The `port-ke-function` skill documents the full porting + fixture-generation + test workflow. Do not
optimize or refactor a kernel in a way that breaks the bit-for-bit match with R.

## Hot-kernel / OpenMP conventions
KEnRef kernels run inside MD refinement (per step, over many atom pairs × models), so inner kernels
in `KEnRef.cpp` are hot. Optimize them — but only after a correct, R-validated baseline exists.

- **`numOmpThreads` parameter**: every kernel takes `int numOmpThreads`, forwarded to
  `num_threads(...)`. `0` means **"use all available threads"** — NEVER treat `0` as serial.
- **Determinism**: prefer deterministic paths (BLAS-style GEMM/GEMV reformulations of accumulator
  loops; OpenMP parallelism over *free* axes with disjoint output writes). Atomic/reduction
  accumulation is ULP-nondeterministic — if a kernel uses it, document that bitwise-reproducible
  output requires `numOmpThreads=1` / `OMP_NUM_THREADS=1`.
- **Performance idioms** (reference kernel: `a_matrix_to_relax`): hoist loop-invariants; cache-block
  on the ColMajor-contiguous (pairs/rows) axis (e.g. 128-row blocks); size-gate parallelism with an
  env-configurable threshold (e.g. `KENREF_RELAX_PARALLEL_THRESHOLD`, default 256) so small inputs
  stay serial; write Eigen array expressions rather than scalar loops or hardcoded intrinsics.
- **OpenMP race-safety**: no `#pragma omp atomic` on an Eigen row (use per-scalar atomics); no
  structured-binding capture inside omp regions (use `std::get<>`); never
  `reduction(+:<dynamic-size Eigen>)` (the private copy default-constructs empty → crash; use
  thread-local + `critical`, sized with `Matrix::Zero(r,c)`); parallelize free axes, not the
  reduction axis.
- After optimizing, re-run the R-ground-truth test at **both** full threads and `OMP_NUM_THREADS=1`;
  passing serial + failing threaded ⇒ a race.

## Git
- **`master` is the development trunk** — new work lands there. Release branches are cut from it and stay
  quiet: `release/1.x` (tags `v1.0.0`, `v1.1.0`) and `release/2.x` (tag `v2.0.0`). The branches are
  named for the MAJOR line, not one version, because each accumulates minor releases.
- **Fixes flow forward, never backward** (PLUMED's model, which this mirrors): a fix goes on the *oldest*
  release branch it applies to, then is merged up the ladder, so nothing on an older branch is ever missing
  from a newer one. New features go straight to `master`.
  ```
  release/1.x  ──►  release/2.x  ──►  master
  ```
  Invariant, checkable any time: `git merge-base --is-ancestor release/1.x release/2.x` and
  `... release/2.x master` must both hold.
- Tags live on the release branches. Install directory names in the deployment tree come from
  `git describe --tags --always`, so an untagged commit deploys as `2.0.0-7-gabc123456` rather than `2.0.0`
  — which is how you tell a release apart from a snapshot.
- Commit/push only when asked.
