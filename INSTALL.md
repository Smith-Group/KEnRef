# Installing KEnRef (from the KEnRef side)

KEnRef is a **library** (`kenref_core`) plus an optional **kenref-gmx** force provider (executables). It can
also build **PLUMED** with the `kenref` module — by **delegating to the PLUMED side** (kenref builds kenref;
each MD engine builds itself). This page covers starting **from the KEnRef repo**. To start from PLUMED
instead, see `src/kenref/install.md` in the PLUMED repo.

> Design: the build managers own the logic. KEnRef's CMake builds kenref_core / kenref-gmx (and find-or-fetch+
> builds the stock GROMACS the force provider needs, delegating to gromacs's own cmake). Anything PLUMED is
> delegated to the plumed scripts. `build.sh` is a shallow shell: it maps flags to `cmake -D…` and hands off.

## Requirements
CMake ≥ 3.24, a build tool (**Make or Ninja** — CMake's default, or set `CMAKE_GENERATOR`), a C++17 compiler,
and git. Eigen: if not found it is **fetched, built and shipped INTO the prefix** (so the install is
self-contained) — unless you point at an external one (`-DEigen3_DIR=…`) or disable shipping
(`-DKENREF_INSTALL_EIGEN=OFF`). GROMACS / PLUMED are only fetched when you explicitly request a component that
needs them (see *Components* below).

**MPI is required only for `--with-gmx`.** `kenref_core` contains no MPI code, and the PLUMED bias reaches
replicas through PLUMED's own `Communicator` (PLUMED supplies its MPI). So a core-only or core+PLUMED build
needs no MPI at all. Only the kenref-gmx force provider and the `energycalc`/`s2calc` tools use MPI directly.

## Components are tri-state — and AUTO never downloads a heavy engine
The three component switches are **`ON` / `OFF` / `AUTO`**, default **`AUTO`**:

- **`ON`** — always build it; a missing sub-dependency is **fetched** (when downloads are on) or the configure
  **FATALs** (when off).
- **`OFF`** — never build it.
- **`AUTO`** — build it **only if its heavy engine (GROMACS / PLUMED) is already present** (found or provided).
  **AUTO never downloads+compiles GROMACS or PLUMED on its own** — no surprise 20-minute builds. To pull in a
  heavy engine you must ask explicitly with **`=ON`**.

Downloads of *sub-dependencies* (Eigen, and — for a committed `=ON` component — GROMACS/PLUMED) are governed by
one switch: **`--download ON|OFF`** (`-DKENREF_FETCH_MISSING`, default `ON`). With `--download OFF` a missing
piece is a clean error telling you what to provide.

## Quick start (a few commands)

```bash
./build.sh -y                        # kenref_core (library)                       -> /usr/local/kenref
./build.sh --with-gmx -y             # + kenref-gmx executables (force GMX on; fetches a stock GROMACS if none)
./build.sh --with-plumed -y          # + PLUMED (kenref module)   — delegated to the plumed side
./build.sh --with-plumed-gromacs -y  # + PLUMED + a batched GROMACS 2025.x         — delegated to build-and-batch.sh
```

`--with-gmx` / `--with-plumed` take an optional value: `--with-gmx=auto` (build gmx *iff* a GROMACS is already
present), `--with-gmx=off`, or bare `--with-gmx` (= `on`, forces it). **Everything installs under ONE prefix**
(`--prefix`, default `/usr/local/kenref`).

## Pointing `--with-gmx` at your GROMACS — one directory

kenref-gmx is a GROMACS **plugin** (a ForceProvider/MDModule), so it needs three things from GROMACS, each
found in exactly one place: the **generated `config.h`** (only in the *build* tree), GROMACS's **internal**
MDModules headers (only in the *source* tree — GROMACS does not install them), and **`libgromacs`** + its CMake
package (only in the *install* tree). An installed GROMACS **alone can never** satisfy kenref-gmx.

You normally pass just **one** option:

```bash
./build.sh --with-gmx -y --gmx-gromacs-build ~/gromacs/cmake-build-release
```

- **`--gmx-gromacs-build DIR`** (`-DGROMACS_BUILD_DIR`) — **GROMACS's own CMake *build* directory**: the one you
  ran `cmake -B` into when you built GROMACS. It contains `CMakeCache.txt` and `src/include/config.h`.
  It is **not** the install prefix (the tree with `bin/gmx` and `share/cmake/gromacs*`) and **not** the source
  checkout. KEnRef reads that cache to derive the matching **source** (`CMAKE_HOME_DIRECTORY`) and **install**
  (`CMAKE_INSTALL_PREFIX`) and the MPI flavor — so the three trees stay consistent and cannot silently drift.
- **`--gmx-gromacs-src DIR`** (`-DGROMACS_SRC_DIR`) — only if you have **no build yet** and want KEnRef to build
  *your* GROMACS checkout for you.
- Pass neither → KEnRef fetches + builds a stock GROMACS (needs `--download ON`).

There is **no `GROMACS_INSTALL_DIR` option**: it is always derived. Nothing is guessed either — no
`../gromacs` sibling assumption and no `/usr/local/gromacs` default, because a wrong guess can silently pick a
different GROMACS **version**.

### Your GROMACS must be built with these two options

```bash
cmake -S gromacs -B gromacs/build -DGMX_MPI=ON -DGMX_INSTALL_LEGACY_API=ON …
```

- **`-DGMX_MPI=ON`** — kenref-gmx runs **multi-simulation** (replica) refinement, and GROMACS supports
  multi-simulations *only* with a real external MPI library (it throws otherwise). A **thread-MPI or serial**
  GROMACS therefore cannot run kenref-gmx; KEnRef detects this from `config.h` (`GMX_LIB_MPI`) and tells you up
  front — `--with-gmx=AUTO` simply skips gmx, `=ON` fails with the reason.
- **`-DGMX_INSTALL_LEGACY_API=ON`** — GROMACS defaults this **OFF**, and then does not install the public
  headers kenref-gmx compiles against. KEnRef warns if your build dir has it off.

Run `./build.sh` with no args for interactive mode. Toolchain comes from the environment
(`CXX=mpicxx CXXFLAGS="-stdlib=libc++" ./build.sh …`). `--build-type` defaults to Release.

## Acceleration (`--accel`, matched to GROMACS's `GMX_SIMD`)
`--accel` selects the SIMD tier; the names mirror GROMACS's `GMX_SIMD` so you can pick the **same** tier for
kenref and gromacs. Left unset, the common x86 tier is auto-detected.

| `--accel`                        | flags                       | Eigen align |
|----------------------------------|-----------------------------|-------------|
| `REFERENCE` / `SSE2` / `SSE4.1`  | scalar / `-msse2` / `-msse4.1` | 16 |
| `AVX_128_FMA`, `AVX_256`, `AVX2_128`, `AVX2_256` | `-mavx…` / `-march=…` | 32 |
| `AVX_512` (`AVX_512_KNL` deprecated) | `-march=skylake-avx512` / … | 64 |
| `ARM_NEON_ASIMD`, `ARM_SVE`, `IBM_VSX` | AArch64 / SVE / POWER | varies |

**Ambiguity to know:** a GROMACS tier name encodes GROMACS's internal *kernel datapath* width, but what
matters for KEnRef is the width the **compiler** enables for Eigen (which sets the alignment/ABI). So
`AVX_128_FMA` and `AVX2_128` — despite "128" — emit 256-bit AVX ⇒ Eigen alignment **32**, the same class as
the `*_256` tiers. `AVX_128_FMA` FMA4-vs-FMA3 and `ARM_SVE` vector length are CPU-dependent.

## SIMD / Eigen ABI safety (important)
kenref_core stores Eigen objects **inside its own containers** (`std::vector<CoordsMatrixType>`, struct
members), so **Eigen's alignment is baked into `libkenref_core`'s ABI**. Building a consumer (kenref-gmx,
PLUMED, or your own code) against an Eigen with a **different alignment** — from a different `--accel`/`-march`,
a different Eigen version, or an Eigen config macro — corrupts memory. KEnRef guards this two ways:

- **Compile time:** including `core/KEnRef.h` pulls in `core/EigenAbiCheck.h`, which `static_assert`s your
  Eigen's `EIGEN_MAX_ALIGN_BYTES` / version against the values kenref_core was built with — a mismatch is a
  **build error**. (Bypass with `-DKENREF_NO_EIGEN_ABI_CHECK` at your own risk.)
- **Runtime:** `cmake --install` runs a check exe that FATALs if the just-built consumer's Eigen alignment
  differs from the linked core's (covers the cross-`--accel` installed-core case). The same check is a gtest.

**Match `--accel` (and the Eigen) between kenref_core and every consumer.**

## What installs where — ONE prefix

Everything lands under a single prefix (default `/usr/local/kenref`), like GROMACS's `/usr/local/gromacs`:

```
/usr/local/kenref/
  bin/       KEnRef, energycalc, s2calc            (only with --with-gmx)
  include/   core/  gmxinterface/  eigen3/
  lib/       libkenref_core.{a,so.*}  libkenref_and_eigen3.a  libkenref_gmxinterface.a  pkgconfig/  cmake/KEnRef/
  share/     pkgconfig/  licenses/  (shipped eigen3/)
  modulefiles/kenref/<version>       (the TCL modulefile — see below)
  env.sh     kenref-build-manifest.txt
```

Override the prefix with `--prefix` (or `-DCMAKE_INSTALL_PREFIX`), or `-DKENREF_INSTALL_BASE=$HOME/.local` for a
no-sudo build. Installing under `/usr/local` uses `sudo` **only** if the prefix isn't writable.

**PLUMED delegation.** `--with-plumed[-gromacs]` sets CMake options; **CMake itself**, at `cmake --install` time
(after kenref_core is installed), clones PLUMED into the build dir and runs the plumed side's `build-only.sh` /
`build-and-batch.sh` (`cmake/InvokePlumed.cmake`), which **reuse** the just-installed kenref — no rebuild, no
loop. Point at a local plumed with `-DKENREF_PLUMED_SRC_DIR` (with `--download OFF` a missing plumed is an error).

**Forwarding arbitrary flags.** Anything after `--` goes straight to CMake, e.g.
`./build.sh -- -G Ninja -DKENREF_EIGEN_INSTALL_DIR=/usr/local`.

## Using an install

```bash
source /usr/local/kenref/env.sh                       # PATH / LD_LIBRARY_PATH / PKG_CONFIG_PATH / CMAKE_PREFIX_PATH
# or environment-modules (modulefile lives UNDER the prefix):
module use /usr/local/kenref/modulefiles && module load kenref/<version>
```

The modulefile is generated at `<prefix>/modulefiles/kenref/<version>` (`-DKENREF_MODULEFILE_DIR` /
`--modulefile-dir`; set empty to skip). To also expose it from a shared system dir, pass
`--link-modulefiles[=DIR]` (default `/usr/local/modulefiles`) — build.sh offers this interactively **only when
the prefix isn't `/usr/local`** (there the modulefile already lives under `/usr/local/modulefiles`). Then
`module use /usr/local/modulefiles && module load kenref/<version>` works too.

## Power user (raw cmake, no script)

```bash
# library + (optionally) gmx executables in ONE configure, ONE prefix. gmx=AUTO builds only if GROMACS is found.
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=/usr/local/kenref \
      -DBUILD_KENREF_CORE=ON -DBUILD_KENREF_GMX=AUTO \
      -DGROMACS_BUILD_DIR=~/gromacs/cmake-build-release    # gromacs's BUILD dir; src+install+MPI derived from it
cmake --build build -j && cmake --install build     # the install runs the SIMD/Eigen ABI gate
```

`CMAKE_BUILD_TYPE` defaults to **Release** when you don't pass one (as GROMACS does), so a raw `cmake` build
is optimised by default; pass `-DCMAKE_BUILD_TYPE=Debug` etc. to change it.

`BUILD_KENREF_CORE` / `BUILD_KENREF_GMX` / `BUILD_KENREF_PLUMED` and `KENREF_WITH_PLUMED` are all tri-state
(`ON`/`OFF`/`AUTO`); `KENREF_FETCH_MISSING` is the download switch. KEnRef's CMake does **not** build PLUMED
directly — `KENREF_WITH_PLUMED=ON` delegates to the plumed side at install time.

## Starting from PLUMED instead
See **`src/kenref/install.md`** in the PLUMED repo — `build-only.sh` (PLUMED only) and `build-and-batch.sh`
(PLUMED + a batched GROMACS), which delegate the kenref build back to this repo.
