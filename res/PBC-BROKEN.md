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

## The reference PDB was broken too — repaired 2026-08-23

`GB3_27_10us.pdb` was byte-identical in all five sets (md5 `693ea2af8dc920348450c52cd61843d2`): one
reference structure, copied five times, and **split in the same way as the tpr** — the same 28.6 Å void
in x. Exactly **12 atoms** were on the wrong side, all side-chain atoms of LYS 10 and THR 11 (including
THR 11's methyl `CG2`/`HG2x`, the kind of proton the restraint data uses).

The 33 guide C-alphas were **not** among them: they span 19.5 / 14.0 / 19.9 Å with a largest gap of
3.2 Å, i.e. whole. That is the only reason the Kabsch fit survived — the fit reads guide atoms only. It
was luck, not design: had one guide atom been in that wrapped fragment, the superposition itself would
have been garbage rather than merely the pair distances.

It is now repaired. Each set ships:

| file | md5 | what it is |
|---|---|---|
| `GB3_27_10us.pdb`            | `3367bafedcddd8a8874f28de45041818` | the repaired reference — what the `.toml`/`.dat` files name |
| `GB3_27_10us_PBC-BROKEN.pdb` | `693ea2af8dc920348450c52cd61843d2` | the original, preserved, named for what it is |

The repair moved **exactly those 12 atoms, by exactly one box vector** (+61.224 Å in x). Every other
atom — and every header, `CRYST1` and B-factor column — is byte-identical, so the `diff` between the two
files is precisely twelve lines. Regenerate it, if it ever needs regenerating, with the tpr for
connectivity:

```
gmx convert-tpr -s md-00_PBC-BROKEN.tpr -n KEnRefAtomIndex.ndx -o protein.tpr   # select Protein
gmx trjconv -s protein.tpr -f GB3_27_10us_PBC-BROKEN.pdb -pbc whole -o GB3_27_10us.pdb
```

Use `-pbc whole`, not `-pbc mol`: `mol` also re-centres molecules, which would move every atom.

**Two mechanisms, cross-checked.** `trjconv` alone is not evidence. `KEnRefMoleculeGraph::makeWhole` —
the bond graph the GROMACS engine uses — was run over the same coordinates and reproduces the result
atom for atom, and that comparison is a committed test rather than a one-off:
`MoleculeGraphTest.BondGraphRepairOfReferencePdbMatchesTrjconv` reads the preserved broken PDB, repairs
it, and requires agreement with the shipped file to within PDB write precision, with every other row
bit-identical. `MoleculeGraphTest.ShippedReferencePdbIsWhole` guards the shipped file itself.

**Why this mattered, per engine:**

* **GROMACS — a no-op, provably.** The live coordinates come from the tpr and are repaired from the
  topology. The reference is used only for the guide atoms, which were already whole, and only their
  coordinates reach the numerics (`KEnRefForceProvider` centres them for Kabsch and nothing else). The
  step-0 energies are unchanged: 3.65155e-05 on the broken GB3 tpr, 7.18869e-08 on ubiquitin.
* **PLUMED — this was the enabler.** PLUMED is never handed a topology, so it must build its spanning
  tree from the *reference*. With a broken reference the wrapped fragment is linked by a 29.17 Å edge
  against a 30.61 Å half-box, minimum-image concludes there is nothing to do, and the break survives by
  1.44 Å. A whole reference makes every tree edge short and the repair correct. Making this file whole
  was the prerequisite for the PLUMED-side fix below, not a cosmetic tidy-up.

## The PLUMED side now repairs too (2026-08-23)

`KEnRefBias` used to take coordinates straight from `getPositions()`, so the PLUMED engine computed on
the torn structure and the two engines disagreed on all five sets. It now mirrors the GROMACS design:
repair once per step, over the whole requested set, before anything is computed.

* At setup, `KEnRefBias_setup.cpp` reads reference positions for **all** of `atoms_` (guide + sub) and
  builds a Euclidean **minimum spanning tree** over them (Prim, O(n²), a few hundred atoms, once). The
  longest reference edge is logged: a few Å is healthy, a large value means the reference is still
  broken and the repair will silently do nothing.
* Each step, `makeRequestedAtomsWhole()` walks that tree parents-before-children, placing each atom at
  the image nearest its parent using PLUMED's own `pbcDistance`. It writes only when the image actually
  changes, so already-whole input comes back **bit-for-bit identical** — the same property that lets the
  GROMACS repair run unconditionally.
* Why an MST and not something simpler: "nearest image to my predecessor in the list" and "nearest image
  to a shared anchor" were both implemented on the GROMACS side and both broke already-whole structures
  (see the notes in `KEnRefMoleculeGraph.h`). The requested atoms are sparse, so list order says nothing
  about proximity; and GB3 reaches further from its centroid than half its shortest box vector, past
  which "nearest image to an anchor" picks the wrong image. An MST joins genuine spatial neighbours, so
  every edge is short and every minimum image along it is unambiguous.

* At the first step it also **refuses** if the longest tree edge is not less than half the smallest box
  width, because past that the nearest image of a child to its parent is simply not the right one and
  the walk re-images atoms that were correctly placed. The box is not known at setup — PLUMED hands it
  over per step — which is why the edge is logged there and judged here. On the committed GB3 sets the
  repaired reference gives 2.73 Å against a 21.65 Å half-width (a factor of 7.9 of margin), and the
  preserved broken one gives 31.49 Å and is rejected.

  This guard matters more than it looks. Without it, the broken reference on these fixtures happens to
  produce the *right* answer anyway: its long edge is 31.49 Å, which exceeds half the 61.224 Å box
  vector it spans, so minimum-image re-images that fragment and lands on the correct result by luck.
  A slightly different geometry — PLUMED's own `makeWhole()` over all 862 protein atoms finds a 29.17 Å
  link instead — falls on the other side and silently produces a torn structure. The guard is what turns
  that coin-flip into a refusal.

The frozen frame in the PLUMED fork changes in three places only: `getLocalModelX` calls
`makeRequestedAtomsWhole()` first, and `getGuideAtomsX`/`fillSubAtomsX` read the repaired buffer.

**Verified 2026-08-23**, `plumed driver` against a step-0 `.trr`:

| case | reference | result |
|---|---|---|
| broken GB3 tpr | repaired `GB3_27_10us.pdb` | runs, **3.65155e-05** — the GROMACS value, to the last digit |
| whole ubiquitin | `00ns.pdb` | runs, **7.18869e-08** — the GROMACS value; the repair is a no-op |
| broken GB3 tpr | `GB3_27_10us_PBC-BROKEN.pdb` | **refuses**: 31.49 Å edge vs a 21.65 Å half-width |

Both gates are committed: `rt_validate.sh` (GROMACS) and `rt_validate_plumed.sh` (PLUMED), 3/3 each.

**Feed `plumed driver` a `.trr`, not a `.gro`.** GRO carries 0.001 nm, and on these fixtures the step-0
SIGMA energy is a near-zero residual: that rounding alone moves 3.65155e-05 to 3.4e-03, a factor of 90,
with nothing actually wrong. `rt_validate_plumed.sh` produces the frame by running the tpr for zero
steps so it carries the tpr's own single-precision coordinates.

## The split-check is a tear test, not a size test (2026-08-23)

`KEnRefDriver`'s refusal is the backstop behind both repairs — repair first, refuse second. It used to
ask whether any atom sat further than half the smallest box vector from the set's centroid. That is a
test of **size**, and it fails on structures that are perfectly whole: after the bond-graph repair the
restrained GB3 atoms reach 21.657 Å against a 21.6561 Å half-box — over by 0.004% — so a 500-step run
aborted at about **step 120** on a molecule that was never torn, identically in serial and under 4-rank
domain decomposition.

It now tests for a **tear**. For each lattice direction the atoms are projected onto their fractional
coordinate, and two extents are compared: the extent as delivered, and the smallest extent any choice of
periodic images could give (1 minus the largest *cyclic* gap between neighbours). For a whole set these
are equal exactly — with nothing wrapped, the largest cyclic gap is the empty space outside the molecule.
They differ only when the largest gap is an interior void, i.e. when the set is split with atoms piled at
both ends and nothing in the middle, and the difference is the size of the tear.

Working in fractional coordinates is what makes it right for the triclinic cells these sets use: whatever
the cell shape, a periodic image differs by a whole number in exactly one fractional component.

Extent was never the property that mattered. Once a set is whole, nothing downstream uses minimum-image
reasoning over the whole molecule: the pair distances have no cutoff and no PBC at all, the Kabsch fit is
global, and `restoreNoJump` images each atom against *its own* position in the previous frame — a
per-step displacement, not a molecular radius. Half the box bounds the reach of minimum-image reasoning,
and both repairs deliberately avoid needing it: they walk short edges, never a molecular radius.

Pinned by `google_tests/testPeriodicSplitCheck.cpp`, entirely on synthetic data so the rule is
independent of any fixture — including
`PeriodicSplitCheckTest.AcceptsALopsidedMoleculeReachingPastHalfTheBox`, which reproduces the step-120
abort, and `HandlesATriclinicCell`, which tears along the third lattice vector of the GB3 dodecahedron.

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
