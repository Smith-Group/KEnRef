# Generated GB3 RELAX test set (`md_double/`)

Self-consistent RELAX (relaxation-rate restraint) exp-data for the GB3 system, runnable through the
offline tool (`energycalc -m RELAX`) and live `gmx mdrun` (KEnRef force provider / PLUMED `KEnRefBias`).
Mirrors `res/sigma/md_double` (same GB3 system + shared `md-00.tpr`/`KEnRefAtomIndex.ndx`/`GB3_27_10us.pdb`)
and the `ke` relaxation pipeline of `demo/relax_deriv_check.R`. Targets are SYNTHETIC (computed from the
structure), so the set is self-validating.

Regenerate: `Rscript gen_relax_gb3.R` (writes `md_double/{1-1..3-3}_{atom_relax,a_coef,groupings,lambda_coef}.csv`
+ `relax_energy.txt`).

- Ensemble = 2 copies of `GB3_27_10us.pdb` (proton-only), i.e. num_models = 2 (run with **2 replicas**).
- Rates = RelaxModel built-in defaults `kens=5e8, kmethyl=1e12, karo=1e4, Dx=Dy=Dz=2.5e8` (no proton_mhz).
- **R ground-truth energy = 4033.6946910272513** (`md_double/relax_energy.txt`).

## Validation (2026-07-21)
- `energycalc -m RELAX` on `md_double` (2 xtc frames) → **4033.69** = R.
- live `gmx mdrun` KEnRef-native (2 replicas, step 0) → 4033.84; PLUMED `-plumed` (rep0) → 4033.837 — agree
  with each other; the ~0.15 offset from R is `md-00.tpr` single-precision coords vs the exact PDB.

## Note: no `md_single`
RELAX needs **>= 2** ensemble members (it restrains relaxation arising from ensemble interconversion). A
1-member set is degenerate — it emits no energy through the tool or live engine — so only `md_double` exists
here (unlike `res/sigma`, which has a meaningful `md_single`).
