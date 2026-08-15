# `md-00_PBC-BROKEN.tpr` — a protein split across a periodic boundary

Every committed GB3 set (`res/sigma/md_{single,double}`, `res/plateau/md_{single,double}`,
`res/relax/md_double`) ships the **same** tpr — byte-identical, md5 `a41df146f7b1b7a3a3d41fd4dff91e05`
— and in it **the protein is split across the periodic boundary in x**.

Renamed from `md-00.tpr` on 2026-08-15 so the name states the fact. It is kept, not fixed: a real
structure that is genuinely broken by periodicity is a valuable and hard-to-fabricate test fixture, and
it is the input that found the bug described below.

## The measurement

From `gmx dump -s`, protein atoms only (global 0..861), box `6.122 x 6.122 x 4.329` nm:

| axis | span | largest interior gap |
|---|---|---|
| x | **6.1164 nm** (the whole box) | **2.8607 nm**, between x=0.196 and x=3.057 |
| y | 2.3049 nm | 0.0499 nm |
| z | 2.9013 nm | 0.0554 nm |

y and z are a normal ~2.5 nm protein. x spans the entire box with a 2.86 nm void in the middle and
atoms piled at both ends: the molecule wraps. Measured directly on the bonds, the longest bond in the
raw coordinates is **6.07 nm**; after make-whole it is **0.184 nm**.

## Why it matters

KEnRef fits and restrains a *whole* molecule: the Kabsch superposition and the pair distances are
global, with no cutoff, so a split molecule silently produces a different energy rather than an error.
Until make-whole was added (2026-08), KEnRef consumed these coordinates as-is, so **every energy
computed from these tprs was computed on a torn protein**. On `res/sigma/md_single` that was the
difference between 0.00125246 (broken) and 3.65155e-05 (repaired).

**The published ubiquitin results are NOT affected.** `alef`/`baaa`/`taaa` and
`res/10nsstart+fitting/t000.00.tpr` were checked with the same method and are whole (largest interior
gaps 0.04–0.15 nm), and make-whole is a verified no-op on them — bond lengths and energies unchanged,
serial and domain-decomposed runs agreeing exactly at step 0.

Note the asymmetry that hid this for so long: **without** domain decomposition GROMACS hands out the
raw, wrapped state, so the molecule arrives broken; **with** domain decomposition its halo
communication delivers locally consistent coordinates, so it arrives whole. The serial path — the
supported one — was the wrong one.

## Regenerated reference energies (2026-08-15)

Step-0 energies from the GROMACS side with make-whole active, `-ntomp 1`, the same parameters
`smoke_test_installs.sh` uses. These supersede any value recorded before that date.

| set | model | replicas | step-0 energy |
|---|---|---|---|
| `res/plateau/md_single` | PLATEAUS | 1 | 0.0018212 |
| `res/plateau/md_double` | PLATEAUS | 2 | 2.87581 |
| `res/sigma/md_single`   | SIGMA    | 1 | 3.65155e-05 |
| `res/sigma/md_double`   | SIGMA    | 2 | 3.65155e-05 |
| `res/relax/md_double`   | RELAX    | 2 | 34.5337 |

For reference, `res/sigma/md_single` gave **0.00125246** before the repair — that is the value of the
torn structure, and it should not be used as a baseline again.

The **offline** figures in `res/sigma/GENERATED_README.md` (301.79 / 236.61 / 124.88) and
`res/relax/GENERATED_README.md` (4033.69) are NOT regenerated here: `energycalc`/`s2calc` read
coordinate frames directly and never pass through the GROMACS force provider, so make-whole does not
apply to them. They still agree with R, which reads the same frames. Whether those *frames* are
themselves periodically broken is a separate, still-open question about the offline path.

## Consequence for the PLUMED side (open)

`KEnRefBias` takes coordinates from `getPositions()` and never calls `makeWhole()`
(`src/kenref/KEnRefBias.cpp`), so the PLUMED engine still computes on the torn structure. Until it is
given the same repair, the two engines will **disagree** on these five sets — which is exactly what
`smoke_test_installs.sh` asserts they must not do. Expect that check to fail on GB3 until PLUMED's
side is fixed; it is a real divergence, not a flake.

## What the fixture is used for

`google_tests/testMoleculeGraph.cpp` reads this tpr directly and asserts:

1. the topology yields one connected fragment of 862 atoms (the protein);
2. the raw coordinates really are broken (longest bond > 5 nm);
3. after `makeWhole()` every bond is a chemical bond (< 0.25 nm);
4. `makeWhole()` is idempotent — a second application changes nothing;
5. `makeWhole()` on already-whole coordinates is a **bit-for-bit** no-op, which is what allows it to
   run unconditionally without perturbing results on well-formed inputs.

If these tprs are ever regenerated, regenerate them **broken**, or move the tests to a fixture that
still is.
