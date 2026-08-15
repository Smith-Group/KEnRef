# Open questions for the PI

Written 2026-08-15. Everything here is **verified in the current source or measured**, with the
evidence given so each item can be checked rather than taken on trust. Nothing here is speculative
design; where something is uncertain, it says so.

Ordered by how much scientific damage the item can do, not by effort.

---

## 1. The restraint contributes NOTHING to the virial — NPT pressure is wrong today

**Status: live defect, affects any constant-pressure run, no domain decomposition required.**

`KEnRefForceProvider::addLocalModelDerivatives` writes forces into
`currentOutput_->forceWithVirial_.force_` and never calls `addVirialContribution`. Traced through
GROMACS 2025.4 `sim_util.cpp`: `calc_virial` (:353) computes the virial from `forceWithShiftForces`
only; forces placed in `forceWithVirial.force_` are summed into the total afterwards (:403) but
contribute to the virial **only** through what a module explicitly adds (:409). KEnRef adds nothing —
`grep -c enerd_ src/gmxinterface/` returns 0, and there is no `addVirialContribution` call anywhere.

**Consequence.** In NPT the barostat sees a pressure that omits the restraint entirely, so the box
relaxes to the wrong volume/density while the restraint is doing work. NVT and NVE are unaffected.

**Questions.** Has any published or in-flight result used NPT with the restraint active? Is NPT a
supported mode at all, or should KEnRef refuse it until this is fixed?

**Fix sketch.** Accumulate `Σ r ⊗ f` over the atoms each rank owns and call `addVirialContribution`
once. It has a domain-decomposition shape (per-rank accumulation over owned atoms only) and would
reuse the machinery already built. Not attempted yet.

---

## 2. The restraint energy never reaches GROMACS

`ForceProviderOutput` carries `enerd_` for exactly this purpose; KEnRef never touches it. The energy
is printed to stdout every 10 steps and nowhere else.

**Consequence.** The restraint energy is absent from `md.edr`, from the conserved-quantity check, and
from every analysis tool that reads energies. A run cannot be assessed by the usual means, and energy
drift caused by the restraint is invisible.

**Question.** Should the restraint appear as its own `md.edr` term? That is the normal expectation for
a bias and would make runs self-documenting. PLUMED's side already reports `kenref.energy`, so the two
engines currently differ in what they expose.

---

## 3. Restarts are silently wrong — no-jump history is not checkpointed

`KEnRefDriver` holds `lastGuide_`, `lastSub_` and `firstStep_` (`KEnRefDriver.h:77-79`) as plain
members. They are the previous frame used for the no-jump correction. Nothing writes them to the
checkpoint and nothing reads them back.

**Consequence.** On restart `firstStep_` is true again, so the first step after a restart performs no
no-jump correction and the history restarts from that frame. A continued run therefore does not
reproduce the uninterrupted one. It fails silently — no warning, plausible numbers.

**Question.** How long are production runs relative to the checkpoint interval, and have any results
come from restarted runs? That determines whether this is urgent or merely untidy.

**Note.** GROMACS provides `MDModulesWriteCheckpointData` for exactly this, but it is delivered only to
built-in modules — which ties this item to §6.

---

## 4. Periodic wholeness — fixed on one engine, refused on the other (context, not a question)

Resolved 2026-08-15; recorded because it reframes what "validated" has meant. Full detail in
`res/PBC-BROKEN.md`.

The protein in all five committed GB3 sets is **split across a periodic boundary**, and KEnRef
consumed those coordinates as given. Every serial run against them was refining a torn protein
(longest bond 6.07 nm raw against 0.184 nm repaired). On `res/sigma/md_single` that is the difference
between 0.00125246 (torn) and 3.65155e-05 (whole).

**The published ubiquitin results are NOT affected** — that system is whole (largest interior gaps
0.04–0.15 nm) and the repair is a verified no-op on it, with serial, 2-rank and 4-rank runs agreeing
exactly. This was checked specifically because the answer mattered.

The uncomfortable part is *why it was not caught*: every check in place — engine-vs-engine agreement,
byte-identity against the previous build, R ground truth — compared two computations that shared the
same broken input. **Agreement between two wrong numbers is indistinguishable from correctness.** The
lasting change is the testing rule now in `CLAUDE.md`: always contrast a broken fixture against a
known-whole one.

**Still open:** the offline path (`energycalc`/`s2calc`) does not go through the force provider and so
gets no repair. Whether the trajectory frames it reads are themselves broken has not been checked. If
they are, the R-ground-truth figures (301.79 / 236.61 / 124.88 / 4033.69) validate the *port* — R reads
the same files — but not the physics. **This is the most likely place for a further surprise.**

---

## 5. Two engines that no longer agree

PLUMED's `KEnRefBias` cannot repair periodic breaks: it has no topology, and PLUMED's own
`makeWhole()` does not help. Measured on our reference: with MOLINFO it walks a spanning tree built
from the *reference* coordinates, which is broken the same way as the frames, so the wrapped fragment
is linked by a 29.17 Å edge against a 30.61 Å half-box — minimum image leaves it exactly where it was,
by 1.44 Å.

So PLUMED now **refuses** rather than biasing a torn structure. Correct, but it means the two engines
no longer produce comparable results on those inputs, and the engine-agreement smoke test cannot run
there until the inputs and reference are made whole.

**Question.** Is two-engine parity still a goal worth its cost? Maintaining it has repeatedly forced
work to be done twice and constrains what can be fixed where (the frozen-frame split exists precisely
to manage this). The alternative is to designate one engine as primary.

---

## 6. Should KEnRef become a GROMACS module? (deferred pending this discussion)

Externally added modules are second-class in GROMACS: `MDModules::add()` stores the module, but every
`subscribe*` dispatcher hard-codes the five built-ins, so `subscribeToSimulationSetupNotifications` is
never called for us. KEnRef works around this by subscribing itself from its own hacked `mdrun.cpp` —
which is how it now gets both the `LocalAtomSetManager` and the topology. That workaround is load-bearing
for domain-decomposition support.

Two separable asks:

- **The small one (recommended, and independent).** Ask GROMACS to iterate `impl_->modules_` in the
  `subscribe*` dispatchers. `IMDModule` declares those methods for all implementers and `add()` accepts
  any `IMDModule`, yet only built-ins receive them — an external module silently gets a subset of the
  interface it implements. There is even a `\todo include field_ in modules_` beside the member. Small,
  self-contained, needs no scientific dependency from GROMACS, and benefits every external module. It
  would also remove most of our `mdrun.cpp` hack.
- **The large one.** Become a first-class built-in like PLUMED. Requirements: an
  `applied_forces/kenref/` module plus a stub; **parameters moved into the `.tpr` via `grompp`**, which
  is a user-facing redesign, not a port; checkpointing (§3) becomes mandatory; a decision on vendoring
  Eigen against the `ACCEL`/`EIGEN_MAX_ALIGN_BYTES` coupling; their test framework, house style and
  licensing; and a named maintainer. Risks: they may simply decline an external scientific dependency;
  KEnRef would then live in-tree in *both* GROMACS and PLUMED, doubling the maintenance surface; and
  the `.mdp`/`.tpr` design is irreversible for users once shipped.

**Question.** Distributed *with* GROMACS, or alongside it? Only the second ask depends on that; the
first is worth filing regardless.

---

## 7. Smaller items worth a decision

- **One core per replica.** The ensemble is one replica per rank by construction. DD now lifts the
  single-rank restriction on the GROMACS side, but the scaling ceiling is really the ensemble size.
  Worth stating what problem sizes are actually targeted.
- **The GB3 fixtures are a poor test set.** All five are the same byte-identical tpr with a broken
  molecule, and the reference PDB is broken too. Regenerating them whole (keeping the broken one as
  the deliberate fixture it now is) would make the smoke tests meaningful.
- **A likely upstream GROMACS bug found in passing.** At `forcerec.cpp:770-780` the `gmx_fatal` for
  "orientation restraints + domain decomposition" appears unreachable: it sits inside
  `if (... && !havePPDomainDecomposition(commrec))` and then tests `if (havePPDomainDecomposition(commrec))`.
  If that reading is right, GROMACS silently runs orientation restraints on broken molecules under DD —
  the same class of bug as §4, in their tree. Not exhaustively verified; worth reporting upstream if
  confirmed, and it would strengthen our standing there.
