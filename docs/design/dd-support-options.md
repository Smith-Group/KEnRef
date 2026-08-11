# Domain-decomposition support in `KEnRefForceProvider` — options

Status: **investigation only**, nothing implemented. Task #23 (plan 3.4b).
Written 2026-08-04 against GROMACS **2025.4** (`~/git/gromacs-4-kenref-2025.4`) and **2026.3**
(`~/git/gromacs-4-kenref`).

## The two defects, precisely

`src/gmxinterface/KEnRefForceProvider.cpp` currently refuses DD outright (`gmx_fatal`, line ~203).
Removing that guard is the acceptance test. Behind it sit two independent problems:

1. **Unguarded `ga2la->findHome()`**, at three sites — `fillSubAtomsX` (:345), `getGuideAtomsX` (:373)
   and `addLocalModelDerivatives` (:289). Each takes a *global* atom index and dereferences the result
   without a null check. Under DD an atom that is not home on this rank returns `nullptr`, so the
   read is undefined behaviour. Worse, the code assumes **every** sub-atom and guide atom is
   resolvable on **this** rank, which is false by construction once atoms are spread over ranks.
2. **`MPI_Gather`/`MPI_Scatter` on `mainRanksComm_`** (:309, :334). That communicator is only valid on
   a simulation's main rank. With >1 PP rank per simulation the non-main ranks reach those calls with
   an invalid communicator.

Both are the *same* root assumption: **one rank per simulation**, so "local == global" and "every rank
is a main rank". The failure is silent then fatal — `nan` energy at step 0, then a crash with nothing
naming KEnRef.

## What GROMACS already provides (verified in both trees)

- **`gmx::LocalAtomSet`** (`src/gromacs/domdec/localatomset.h`) — the purpose-built answer. You register
  a set of global indices once with `LocalAtomSetManager::add(...)`; GROMACS then maintains, across
  every repartition:
  - `localIndex()` — local indices of the set's atoms **this rank owns**
  - `collectiveIndex()` — for each local atom, its slot in the **set's own global ordering**
  - `globalIndex()`, `numAtomsLocal()`, `numAtomsGlobal()`
  `LocalAtomSetManager*` is already delivered to `subscribeToSimulationSetupNotifications`, which
  KEnRef already implements.
- **`MDModulesAtomsRedistributedSignal`** — emitted after each DD repartition (and without DD after
  atoms are put in the box). Note: it is carried on the **`simulationSetupNotifier_`**, *not* a
  separate run notifier — the plan's wording ("delivered through `subscribeToSimulationRunNotifications`")
  is wrong on this point. With `LocalAtomSet` you may not need this signal at all.
- **The in-tree idiom for "collect a subset across ranks"**, from `applied_forces/nnpot`
  (`nnpotforceprovider.cpp`) — the closest analogue to KEnRef, since an NN potential also needs a
  whole atom set before it can compute:
  ```cpp
  gIdx = inpAtoms_->globalIndex()[inpAtoms_->collectiveIndex()[i]];
  // ... zero a full-size array, write ONLY locally-owned atoms into their collective slot ...
  gmx_sum(3 * numInput, positions_.data()->as_vec(), cr_);   // 2025.4
  mpiComm_.sumReduce(3 * numInput, positions_.data()->as_vec());  // 2026.3
  ```
  Because non-owned slots are zero, an **allreduce-sum reconstructs the complete set on every rank**.
  No `Gatherv`, no displacement arrays, no root rank. Forces come back the same way.

## Options

### A. `LocalAtomSet` + allreduce-sum (the nnpot idiom) — **recommended**

Register the sub-atom set and the guide-atom set as two `LocalAtomSet`s at setup. Per step: zero the
full-size sub/guide arrays, fill only `localIndex()` entries at their `collectiveIndex()` slot,
allreduce-sum over the simulation's PP ranks. Every rank then holds the full, correctly-ordered set —
which is exactly what the existing pipeline (no-jump → Kabsch → cross-replica gather) expects.

For forces, symmetric: every rank has all derivatives, and each adds only the atoms in its own
`localIndex()`. The `addLocalModelDerivatives` loop becomes a loop over `localIndex()`, so the null
deref disappears rather than being null-checked.

The cross-replica `MPI_Gather`/`MPI_Scatter` on `mainRanksComm_` must then be guarded to run **only on
each simulation's main rank**, with the result broadcast down the simulation's own communicator.

- **+** Idiomatic and in-tree-proven; drops all three `ga2la` pokes and the repartition bookkeeping.
- **+** Portable across 2025.4/2026.x behind the existing `GMX_VERSION` shim (`gmx_sum` vs `sumReduce`).
- **+** Reduces GROMACS-internal coupling (risk R4) — `LocalAtomSet` is `\inlibraryapi`, `ga2la` is not.
- **−** One allreduce per step per simulation over `3*(nSub+nGuide)` reals. Negligible at KEnRef's
  subset sizes; measure before assuming.
- **−** Allreduce ordering is not bitwise-reproducible across rank counts — see "Numerical" below.

### B. `LocalAtomSet` + `Gatherv`/`Scatterv` to the simulation main rank

Same registration, but assemble on the main rank only. Less bandwidth than A and keeps the "main rank
computes" shape of the current code.

- **+** Slightly cheaper; only the main rank holds the full set.
- **−** Needs displacement/count arrays rebuilt at every repartition; more code, more to get wrong.
- **−** Forces need an explicit scatter back. A already gets this for free.

### C. Collect the whole global state (`dd_collect_vec`-style)

Let GROMACS gather all coordinates onto the main rank, then index globally as today.

- **+** Smallest diff to the existing indexing logic.
- **−** Moves **every** atom every step to serve a subset — wrong complexity for large systems.
- **−** Leans on internal DD APIs; worsens exactly the coupling risk R4 flags.

### D. Keep the guard (status quo)

Document one-rank-per-replica as the only supported mode; scale with `-ntomp` and `-multidir`.

- **+** Zero risk, zero work; correct today.
- **−** Leaves the known one-core-per-replica ceiling in place, which is a standing risk-register item.

## Numerical consequence — must be decided before implementing

An allreduce (A) or a gather (B) sums/orders contributions in a rank-count-dependent way. Any
reduction KEnRef performs afterwards can therefore differ in the last ULP between 1-rank and N-rank
runs. That collides with the current acceptance bar used in 3.4a ("energies identical **and** final
`md.gro` byte-identical" vs the 2025.4 build).

So the DD acceptance criterion has to be restated deliberately: byte-identity is achievable *at equal
rank count*, but across rank counts the honest bar is agreement to a stated tolerance. Decide this
with the user before writing code — it changes what the validation harness asserts.

## Findings 2026-08-11 (verified) — three defects, and two things that make stage 1 easier

**There is a THIRD defect, not two.** Identity is *not* broken on the GROMACS side: the force provider
already reads `ms->numSimulations_` / `ms->simulationIndex_` (`KEnRefForceProvider.cpp:146-147`), and
those are plain ints valid on **every** rank (`mainRanksComm_` is the only main-rank-only member). But
that is exactly why `if (simulationIndex_ == 0)` at `:311` and `:330` is wrong under DD: it is true on
*every rank of simulation 0*, so all of them would redundantly run the model and drive the cross-replica
collectives. It must become "replica 0 **and** this simulation's main rank".

| # | defect | fix |
|---|---|---|
| 1 | `ga2la->findHome()` dereferenced unguarded at `:289`, `:345`, `:373` | `LocalAtomSet` + cross-rank assembly |
| 2 | `MPI_Gather`/`MPI_Scatter` on `mainRanksComm_` from all ranks (`:309`, `:334`) | gate on being the simulation's main rank |
| 3 | `simulationIndex_ == 0` true on every rank of replica 0 | add the main-rank condition |

**Two facts that make stage 1 cheaper than feared:**

- **The indexing has no step-0 dependency.** `fillParamsStep0(homenr, numSimulations, forceProviderInput)`
  **uses neither `homenr` nor `forceProviderInput`** — they appear only in the signature. Everything
  comes from `Settings::` and file reads. So it can move to setup time, which is where
  `LocalAtomSetManager*` is delivered and the only time `add()` may be called.
- **`subscribeToSimulationSetupNotifications` is an empty TODO** (`KEnRefMDModule.cpp`), so there is no
  existing behaviour to preserve there — the `LocalAtomSetManager*` subscription is new code, not a
  rewrite.

**This work is NOT blocked on plumed/plumed2#1443.** That issue concerns PLUMED's identity plumbing;
GROMACS already hands identity over as per-rank scalars. (It is also evidence *for* option D upstream:
the consumer that reads scalars is the one that is correct.)

## CORRECTION 2026-08-11 — `LocalAtomSet` is not reachable by default, and is not required

**The recommendation above was wrong about the mechanism.** `KEnRefMDModule::subscribeToSimulationSetupNotifications`
is **never called by GROMACS**. `MDModules::add()` stores the module in `impl_->modules_`, and that
vector is iterated in exactly one function, `initForceProviders`. Every `subscribe*` dispatcher in
`mdmodules.cpp` hard-codes the five built-ins (`densityFitting_`, `qmmm_`, `colvars_`, `plumed_`,
`nnpot_`). PLUMED receives notifications because it is a built-in; an externally added module never
will. KEnRef's own `mdrun.cpp` cannot supply the manager either — it holds no `LocalAtomSetManager`
and calls no `subscribe*`; the manager is created inside `Mdrunner::mdrunner()`.

**Step 1 (implemented) — subscribe ourselves.** The notifier is a plain broadcast: subscribing only
registers callbacks and nothing checks who registered, and `MDModules::notifiers()` returns the very
object the runner emits on (`runner.cpp:1150`, `setupNotifier.notify(&atomSets)` at `:1828`). So
`mdrun.cpp` now calls `subscribeToSimulationSetupNotifications` itself, with a `const_cast` because
`notifiers()` returns `const&` while subscribing mutates. This unlocks the whole setup payload
(`LocalAtomSetManager*`, `t_commrec&`, `gmx_multisim_t*`, `MDLogger`, `gmx_mtop_t&`,
`MDModulesAtomsRedistributedSignal`) and could later **replace** `setSimulationContext`, shrinking the
mdrun hack rather than growing it.

**`LocalAtomSet` is an optimisation, not a necessity.** The zero-fill + reduce design works with plain
**null-checked `findHome()`**: each rank fills only the slots it owns, the rest stay zero, and the sum
reconstructs the set. Exactly one writer per slot (findHome is home-only) keeps it bit-exact, and a
fresh lookup each step is repartition-safe by construction. The current code already calls `findHome`
every step — it just dereferences the result unguarded.

Hence the two milestones:

- **A** (no `mdrun.cpp` dependency): guarded `findHome` + zero-fill + reduce. Fixes correctness; must be
  bit-identical at one rank per simulation.
- **B** (needs step 1): swap to `localIndex()`/`collectiveIndex()`. Complexity drops from
  O(nSub) lookups *per rank* to O(nSub) *total*, and it moves us off `ga2la` (risk R4). Acceptance test
  is bit-identity against A. **Trap:** `localIndex_`/`collectiveIndex_` are vectors cleared and
  re-`push_back`ed at every repartition, so re-fetch the `ArrayRef` each step and hoist only *within* a
  step — never cache it across steps.

## Suggested sequence

1. Register the two `LocalAtomSet`s; keep the `gmx_fatal` guard in place — no behaviour change yet.
2. Rewrite the three `ga2la` sites in terms of `localIndex()`/`collectiveIndex()`; still guarded. At
   one rank per simulation this must be **bit-identical** to today — that is the free regression test.
3. Add the allreduce and the main-rank guard on the cross-replica collectives.
4. Remove `gmx_fatal`. Validate with the `dd_validate.sh` harness (>1 rank *per simulation*, all three
   models, 2025.4 and 2026.x).

Step 2 is the valuable one: it is verifiable without any DD run at all.
