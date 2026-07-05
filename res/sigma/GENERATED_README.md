# Generated SIGMA test input sets (`gen_<N>member/`)

Clean, self-consistent SIGMA spectral-density input folders for development testing & debugging,
regenerated from the `ke` demo (`~/git/ke/demo/sigma_deriv_check.R`) via `gen_multi.R` (copied here).

Each `gen_<N>member/expdata/` folder is built for an **N-member ensemble** and, crucially, must be
driven with **exactly N input models** (num_models == N): `energycalc --input m1.xtc m2.xtc ... mN.xtc`,
or a PLUMED run with N replicas (`mpirun -n N -multidir ...`). Supplying the wrong number of models is
now rejected cleanly (`SigmaModel` validates once: "number of input models does not match the number of
ensemble members ..."), instead of the previous heap corruption.

Coordinates: the ensemble is `ke`'s `extdata/gb3/2lum_subset.pdb.gz` (proton-only, 423 atoms), which is
also checked in as `res/google_tests/2lum_subset_proton.pdb`. Energy is computed over ensemble members
`c(1,3,5,...)[1:N]`; synthetic cross-relaxation rates were built from members `c(2,4,6,...)[1:N]`.

## R ground-truth energies (proton_mhz=700, kens=1/2e-9, kc=1/4e-9, default k/n)

| folder            | members (num_models) | grouping width (1-1) | R sigma_energy |
|-------------------|----------------------|----------------------|----------------|
| `gen_1member/`    | 1 (single sim)       | 1                    | 301.7856126    |
| `gen_2member/`    | 2 (2-member ensemble)| 2                    | 236.610388     |
| `gen_3member/`    | 3                    | 3                    | 124.8784614    |

`gen_3member` reproduces the unit-test ground truth (`testCoordArrayToSigmaEnergy` = 124.8785).

## "single" vs "double" (naming clarification)
`single`/`double` = **number of ensemble members** (1 vs 2 exchanging conformational states), NOT
"double-exponential". A 2-member set adds the `kens` interconversion component (`a_coef` gains a 2nd
column, grouping width 2); a 1-member set has only the rigid `"0"` component (width 1).

Regenerate: `Rscript gen_multi.R` (writes these folders).
