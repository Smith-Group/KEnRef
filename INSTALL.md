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
MPI, and git. Eigen: if not found, it is **fetched, built (with the same build type as KEnRef) and shipped
INTO the kenref prefix** (so the install is self-contained) — unless you point at an external one
(`-DEigen3_DIR=…`) or disable shipping (`-DKENREF_INSTALL_EIGEN=OFF`). For `--with-gmx`, either
provide a stock **GROMACS** (source+build+
install) or let CMake fetch+build one (heavy, one-time). For `--with-plumed[-gromacs]`, PLUMED is cloned into
the build dir (and, for the gromacs variant, a GROMACS 2025.x is fetched by the plumed side).

## Quick start (a few commands)

```bash
./build.sh -y                       # kenref_core (library)                 -> /usr/local/kenref
./build.sh --with-gmx -y            # + kenref-gmx executables              -> /usr/local/kenref-gmx
                                    #   (CMake fetches+builds a stock GROMACS if you don't provide one)
./build.sh --with-plumed -y         # + PLUMED (kenref module)  — delegated to the plumed side
./build.sh --with-plumed-gromacs -y # + PLUMED + a batched GROMACS 2025.x   — delegated to build-and-batch.sh
```

Provide your own stock GROMACS to `--with-gmx` (skips the fetch+build):

```bash
./build.sh --with-gmx -y \
     --gmx-gromacs-src ~/gromacs --gmx-gromacs-build ~/gromacs/build --gmx-gromacs-install /usr/local/gromacs
```

Run `./build.sh` with no args for interactive mode. Toolchain comes from the environment
(`CXX=mpicxx CXXFLAGS="-stdlib=libc++" ./build.sh …`). `--build-type` defaults to Release; `--accel`
**auto-detects** the machine's SIMD.

## What installs where

| Component | Contents | Default location |
|---|---|---|
| **kenref** (library) | `libkenref_core.{a,so}`, `libkenref_and_eigen3.a`, headers, `.pc`, cmake config | `/usr/local/kenref` |
| **kenref-gmx** (executables) | `KEnRef`, `energycalc`, `s2calc`, `libkenref_gmxinterface.a` | `/usr/local/kenref-gmx` |
| **PLUMED** (+ GROMACS) | via delegation | the plumed side's prefixes (see its `install.md`) |

**One set of core libraries.** `--with-gmx` builds in two phases — core → `/usr/local/kenref`, then the
executables → `/usr/local/kenref-gmx` **consuming** that core (never a duplicate `libkenref_*`). Override the
roots with `--prefix` (library) and `--gmx-prefix` (executables), or `-DKENREF_INSTALL_BASE=$HOME/.local` for a
no-sudo build. Installing under `/usr/local` uses `sudo` **only** if the prefix isn't writable.

**PLUMED delegation.** `--with-plumed[-gromacs]` just sets CMake options (`-DKENREF_WITH_PLUMED[_GROMACS]=ON`).
**CMake itself** then, at `cmake --install` time (after kenref_core is installed), clones PLUMED into the build
dir and runs the plumed side's `build-only.sh` / `build-and-batch.sh` (`cmake/InvokePlumed.cmake`). Those
**reuse** the just-installed kenref (found via the installed `env.sh` the hook sources) — no rebuild, no loop —
and do the plumed/gromacs work in their own build dir. Point at a local plumed with `-DKENREF_PLUMED_SRC_DIR`.

**Forwarding arbitrary flags.** Anything after `--` goes straight to CMake, e.g.
`./build.sh -- -G Ninja -DKENREF_EIGEN_INSTALL_DIR=/usr/local`.

## Using an install

Every install prints a summary and writes `<prefix>/env.sh`, and **always generates a TCL modulefile**:

```bash
source /usr/local/kenref/env.sh                       # PATH / LD_LIBRARY_PATH / PKG_CONFIG_PATH / CMAKE_PREFIX_PATH
# or environment-modules:
module use /usr/local/modulefiles/kenref && module load kenref/<version>
```

Modulefile location is `-DKENREF_MODULEFILE_DIR` / `--modulefile-dir` (default
`<install-base>/modulefiles/<component>/<version>`; set empty to skip). The **kenref-gmx** env also references
the **kenref** lib prefix so the executables resolve `libkenref_core` at runtime; the core env includes the
external `eigen3.pc` so a downstream PLUMED build prefers `kenref_core + external eigen3`.

Paths to add manually: `PATH` += `…/kenref-gmx/bin`; `LD_LIBRARY_PATH` += `…/kenref/lib`;
`PKG_CONFIG_PATH`/`CMAKE_PREFIX_PATH` += `…/kenref` (+ `…/kenref-gmx`).

## Power user (raw cmake, no script)

```bash
# library (phase 1):
cmake -S . -B build-core -DBUILD_KENREF_CORE=ON -DBUILD_KENREF_GMX=OFF -DCMAKE_INSTALL_PREFIX=/usr/local/kenref
cmake --build build-core -j && cmake --install build-core

# kenref-gmx executables (phase 2, consumes the installed core). Provide the 3 GROMACS dirs, or omit them and
# CMake find-or-fetch+builds a stock gromacs (kenref_provide_gromacs) into the build dir:
cmake -S . -B build-gmx -DBUILD_KENREF_CORE=OFF -DBUILD_KENREF_GMX=ON \
      -DCMAKE_PREFIX_PATH=/usr/local/kenref -DCMAKE_INSTALL_PREFIX=/usr/local/kenref-gmx \
      -DGROMACS_SRC_DIR=… -DGROMACS_BUILD_DIR=… -DGROMACS_INSTALL_DIR=…    # (optional; else auto)
cmake --build build-gmx -j && cmake --install build-gmx
```

Building CORE+GMX in **one** build co-locates them under a single prefix (CMake has one install prefix) and
warns — use the two phases above for one canonical copy of the libraries. KEnRef's CMake does **not** build
PLUMED; use `build.sh --with-plumed[-gromacs]` (or the plumed side directly) for that.

## Starting from PLUMED instead
See **`src/kenref/install.md`** in the PLUMED repo — `build-only.sh` (PLUMED only) and `build-and-batch.sh`
(PLUMED + a batched GROMACS), which delegate the kenref build back to this repo.
