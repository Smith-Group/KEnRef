# Where KEnRef installs — the three tiers

A KEnRef build can land in one of three places, and they exist for different audiences. Mixing them up is
the single easiest way to test the wrong binary, so each has a distinct top-level name.

| tier | who it is for | location | layout |
|---|---|---|---|
| **development** | you, right now | the working repos / CLion build dirs | whatever your IDE does |
| **staging** | internal testing | `/smithlab/opt/<pkg>-dev/…` | hierarchical: one directory level per dimension |
| **deployment** | users | `/smithlab/opt/<pkg>/…` | flat: the whole config encoded in one directory NAME |

The `-dev` suffix is not decoration. It is what stops the testing matrix from occupying the names users
load, so `module load plumed/…` reaches a release and `module load plumed-dev/…` reaches the matrix.

> **A tier is never copied into another. It is COMPILED again.**
> KEnRef's `.pc` files bake `prefix=<configure-time path>` into their **content**, and neither `cp -a` nor
> `cmake --install --prefix <other>` rewrites file content. A copied install would advertise the ORIGINAL
> location's `-I`/`-L`, so both "installs" would silently compile and link against the same files — a tree
> that looks right and tests the wrong thing, and breaks outright if the first is removed. (GROMACS is
> different: its rpath is `$ORIGIN`-relative, so it genuinely relocates. KEnRef is not.)

---

## The user's own install (neither tier)

Someone who just builds KEnRef gets **one prefix with everything in it** — see `INSTALL.md`. That is the
default and needs none of this document:

```
/usr/local/kenref/          bin/ include/ lib/ share/ modulefiles/ env.sh
/usr/local/                 PLUMED (autotools default)
/usr/local/gromacs/         GROMACS (its own default)
```

The two tiers below exist only because *this* machine keeps many configurations side by side.

---

## Staging — `/smithlab/opt/<pkg>-dev/`

Every dimension that identifies a build is a directory level, so any one of them can be rebuilt and tested
independently. Built by `~/scripts/build-all-workflow.sh`.

```
kenref-dev/<kenref_ver>/<buildtype>/<accel>
kenref-gmx-dev/<kenref_ver>/<gmx_year>/<gmx_ver>/<buildtype>/<accel>
gromacs-4-kenref-dev/<gmx_year>/<gmx_ver>/<buildtype>/<accel>
plumed-dev/<plumed_ver>/<kenref_ver>/<buildtype>/<accel>
gromacs-4-plumed-dev/<gmx_year>/<gmx_ver>/<buildtype>/<accel>
modulefiles/plumed-dev/<plumed_ver>/<kenref_ver>/<buildtype>/<accel>
```

Dimensions: kenref `1.0.0`/`2.0.0`; buildtype `debug`/`release`/`relwithdebinfo`; accel
`AVX_512`/`AVX2_256`/`AVX_256`; GROMACS `2025.4`; PLUMED `master`. **63 builds.**

Two asymmetries worth knowing:

- **`kenref-gmx` is keyed by GROMACS version but `kenref` is not.** The core is engine-independent; the gmx
  interface compiles against GROMACS's *internal* headers and its generated `config.h`, so it is specific
  to one GROMACS build.
- **`gromacs-4-plumed` carries no plumed/kenref dimension.** `plumed patch --mode runtime` does not link
  PLUMED into GROMACS — the kernel is `dlopen`'d at run time via `$PLUMED_KERNEL`. One batched GROMACS
  therefore serves every plumed/kenref combination; adding those dimensions would force byte-identical
  rebuilds.

Staging keeps **stable** version levels rather than `git describe` output, so smoke tests and module loads
don't churn on every commit.

---

## Deployment — `/smithlab/opt/<pkg>/`

What users load. Built by `~/scripts/deploy-all-workflow.sh`. Three rules define it:

1. **Release only.** The build type is not a dimension and never appears in a name.
2. **Named, not nested.** The config is encoded in one directory name.
3. **Append-only.** Nothing here is ever deleted or overwritten. A failed or superseded build stays where
   it is and its replacement is installed *beside* it; the scripts refuse to write into an existing prefix
   rather than clobber a live install.

### The grammar

```
<primary-version>_<accel>[-<relation><version>]…
```

`_` separates config components of the same product; `-` introduces a **relation** to another product,
carrying only what the leading part does not already imply (the GROMACS relation is just its version — the
accel is shared). Later modifiers append to the leading part: `2.0.0_AVX_512_cuda13`.

> ACCEL tokens contain `_` themselves (`AVX_512`), so a name is **positional** — never split it on `_`.

### The trees

```
kenref/<kenref>_<accel>                    KEnRef core + the kenref-gmx interface, one prefix
kenref/<kenref>_<accel>-gmx-<gmxver>       the stock GROMACS that build was made against
plumed/<plumed>_<accel>-kenref<kenref>     PLUMED (its kernel LINKS that kenref)
gromacs-4-plumed/<gmxver>_<accel>          the PLUMED-batched GROMACS
modulefiles/{kenref,plumed}/<same name>    modulefiles, flat, symlinked into the install
```

A worked example, with `2.0.0` and a PLUMED fork 778 commits past `v2.10.1`:

```
/smithlab/opt/kenref/2.0.0_AVX_512
/smithlab/opt/kenref/2.0.0_AVX_512-gmx-2025.4
/smithlab/opt/plumed/2.10.1-778-gb1dc713b2_AVX_512-kenref2.0.0
/smithlab/opt/gromacs-4-plumed/2025.4_AVX_512
```

The kenref pair sorts adjacent on purpose: their **ACCELs must match**, or Eigen's alignment differs
between them and Eigen objects crossing the boundary are misaligned. The install-time SIMD/Eigen ABI gate
catches it, but the naming is what makes the pairing obvious in the first place.

### Versions come from git

`git describe --tags --always --abbrev=9`, with the leading `v` stripped:

| HEAD is… | token | means |
|---|---|---|
| at a tag | `2.0.0` | a release |
| past a tag | `2.0.0-7-gabc123456` | 7 commits past `v2.0.0`, at `abc123456` — a **snapshot** |
| no tag reachable | `abc123456` | an untagged line |

That distinction only works if the tags sit on the branch where work happens — which is why `master` is the
trunk and the release branches carry the tags (see `CLAUDE.md`). The PLUMED fork is permanently far past
its last upstream tag, so a describe string *is* its normal release name, not an anomaly.

**A dirty tree is refused.** A directory named after a commit must contain that commit, or the name is a
claim nobody can reproduce. `KN_ALLOW_DIRTY=1` overrides and appends `-dirty`, so the name still tells the
truth. Staging has no such rule — it is meant to be built from work in progress.

---

## Checking an install

```bash
# staging
module load plumed-dev/master/2.0.0/release/AVX_512
source /smithlab/opt/kenref-gmx-dev/2.0.0/2025/2025.4/release/AVX_512/env.sh

# deployment
module load plumed/2.10.1-778-gb1dc713b2_AVX_512-kenref2.0.0
source /smithlab/opt/kenref/2.0.0_AVX_512/env.sh
```

`./smoke_test_installs.sh [staging|deploy]` drives the committed `res/` sets through live MD from **both**
sides (the KEnRef force provider and `gmx mdrun -plumed`) against whichever tier you name, and resolves
these paths itself. `./default_install_check.sh` answers the separate question of whether the *documented
default* install works out of the box, at `/usr/local`, with no options at all.
