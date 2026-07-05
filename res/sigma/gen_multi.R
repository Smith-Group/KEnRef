# Generate clean, self-consistent SIGMA exp-data folders for 1, 2, 3 ensemble members.
# Writes each to res/sigma/gen_<m>member/{expdata, ensemble.pdb} + prints the R ground-truth energy.
suppressMessages(pkgload::load_all("/home/amr/git/ke", quiet=TRUE))
OUTROOT <- "/home/amr/CLionProjects/KEnRef/res/sigma"

pdb_path <- system.file("extdata", "gb3", "2lum_subset.pdb.gz", package = "ke")
pdb2lum  <- read_ensemble(pdb_path, proton_only = TRUE)

gen_one <- function(m) {
  members_energy    <- seq(1, by = 2, length.out = m)   # 1 / 1,3 / 1,3,5
  members_synthetic <- seq(2, by = 2, length.out = m)   # 2 / 2,4 / 2,4,6
  coord_array     <- pdb2lum[, , members_energy,    drop = FALSE]
  coord_synthetic <- pdb2lum[, , members_synthetic, drop = FALSE]

  base_rates <- c(kens = 1/2e-9); kc <- 1/4e-9; proton_mhz <- 700
  base_rate_mat <- rate_mat_simple(base_rates, dimnames(coord_synthetic)[[3]])
  ke_data <- make_ke_data(coord_synthetic, base_rate_mat, base_rates, kc, proton_mhz=proton_mhz, mix_times=numeric())

  dist_list    <- apply(coord_synthetic, 3, function(x) as.matrix(dist(t(x))), simplify=FALSE)
  min_dist_mat <- Reduce(pmin, dist_list)   # element-wise min across members; preserves dimnames for any m
  eq <- ke_data$equiv_list
  emd <- matrix(NA_real_, length(eq), length(eq), dimnames=list(names(eq), names(eq)))
  eq_single <- unlist(eq[sapply(eq, length) == 1])
  emd[names(eq_single), names(eq_single)] <- min_dist_mat[eq_single, eq_single]
  eq_multi <- eq[sapply(eq, length) > 1]
  emp <- expand.grid(names(eq), names(eq_multi), stringsAsFactors=FALSE)
  for (i in seq_len(nrow(emp))) {
    md <- min(min_dist_mat[eq[[emp[i,1]]], eq[[emp[i,2]]]])
    emd[emp[i,1],emp[i,2]] <- emd[emp[i,2],emp[i,1]] <- md
  }
  idx <- which(emd < 5 & upper.tri(emd), arr.ind=TRUE); idx <- idx[order(idx[,1], idx[,2]),]
  sigma_pairs <- matrix(rownames(emd)[idx], ncol=2)
  sdl <- make_spec_den_data(ke_data, sigma_pairs, perm_internal=TRUE)
  sig_syn <- coord_array_to_sigma(aperm(coord_synthetic, c(2,1,3)), ke_data$rates, sdl, proton_mhz)
  for (i in seq_along(sdl)) {
    sdl[[i]]$atom_pairs <- data.frame(sdl[[i]]$atom_pairs, sigma=NA_real_)
    sdl[[i]]$atom_pairs[seq_along(sig_syn[[i]]), "sigma"] <- sig_syn[[i]]
  }
  EXP <- file.path(OUTROOT, sprintf("gen_%dmember", m), "expdata")
  dir.create(EXP, showWarnings=FALSE, recursive=TRUE)
  for (i in seq_along(sdl)) {
    tn <- names(sdl)[i]
    write.csv (sdl[[i]][["atom_pairs"]],  file.path(EXP, paste0(tn, "_atom_pairs.csv")),  row.names=FALSE, na="")
    write.table(sdl[[i]][["groupings"]],  file.path(EXP, paste0(tn, "_groupings.csv")),   sep=",", row.names=FALSE, col.names=FALSE)
    write.csv (sdl[[i]][["a_coef"]],      file.path(EXP, paste0(tn, "_a_coef.csv")),      row.names=FALSE)
    write.csv (sdl[[i]][["lambda_coef"]], file.path(EXP, paste0(tn, "_lambda_coef.csv")))
  }
  res <- coord_array_to_sigma_energy(aperm(coord_array, c(2,1,3)), ke_data$rates, sdl, proton_mhz, gradient=FALSE)
  se <- if (is.list(res)) res[[1]] else res
  cat(sprintf("MEMBERS=%d  R_SIGMA_ENERGY=%.10g  N_PREFIXES=%d  N_ATOMS=%d\n", m, se, length(sdl), dim(coord_array)[2]))
}
for (m in 1:3) try(gen_one(m))
