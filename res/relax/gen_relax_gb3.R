# Generate a self-consistent GB3 RELAX exp-data set (res/relax/md_{single,double}) that is
# runnable through the offline tool (energycalc -m RELAX) and live gmx/plumed, mirroring how
# res/sigma/md_{single,double} were built (GB3 ensemble = N copies of GB3_27_10us.pdb) and the
# relaxation pipeline of ke demo/relax_deriv_check.R. Targets are synthetic (from the structure).
# Rates (kens=5e8,kmethyl=1e12,karo=1e4,Dx=Dy=Dz=2.5e8) == RelaxModel's built-in defaults.
suppressWarnings(suppressMessages(pkgload::load_all("/home/amr/PycharmProjects/ke", quiet=TRUE)))

REPO <- "/home/amr/CLionProjects/KEnRef"
PDB  <- file.path(REPO, "res/sigma/md_double/GB3_27_10us.pdb")

gen_relax <- function(n_members, out_dir) {
  dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
  ensemble        <- read_ensemble(rep(PDB, n_members), proton_only=TRUE)
  coord_array     <- ensemble[,,1:n_members, drop=FALSE]
  coord_synthetic <- coord_array

  base_rates <- c(kens = 1/2e-9); kmethyl <- 1/1e-12; karo <- 1/100e-6; d_rot <- 1/4e-9; proton_mhz <- 700
  base_rate_mat <- rate_mat_simple(base_rates, dimnames(coord_synthetic)[[3]])
  ke_data <- make_ke_data(coord_synthetic, base_rate_mat, base_rates, kc=d_rot, proton_mhz=proton_mhz, mix_times=numeric())

  dist_array   <- simplify2array(apply(coord_synthetic, 3, function(x) as.matrix(dist(t(x))), simplify=FALSE))
  min_dist_mat <- apply(dist_array, 1:2, min)
  equiv_min_dist_mat <- matrix(NA_real_, length(ke_data$equiv_list), length(ke_data$equiv_list),
                               dimnames=list(names(ke_data$equiv_list), names(ke_data$equiv_list)))
  equiv_single <- unlist(ke_data$equiv_list[sapply(ke_data$equiv_list, length)==1])
  equiv_min_dist_mat[names(equiv_single), names(equiv_single)] <- min_dist_mat[equiv_single, equiv_single]
  equiv_multi       <- ke_data$equiv_list[sapply(ke_data$equiv_list, length)>1]
  equiv_multi_pairs <- expand.grid(names(ke_data$equiv_list), names(equiv_multi), stringsAsFactors=FALSE)
  for (i in seq_len(nrow(equiv_multi_pairs))) {
    md <- min(min_dist_mat[ke_data$equiv_list[[equiv_multi_pairs[i,1]]], ke_data$equiv_list[[equiv_multi_pairs[i,2]]]])
    equiv_min_dist_mat[equiv_multi_pairs[i,1], equiv_multi_pairs[i,2]] <-
      equiv_min_dist_mat[equiv_multi_pairs[i,2], equiv_multi_pairs[i,1]] <- md
  }
  idx <- which(equiv_min_dist_mat < 5 & upper.tri(equiv_min_dist_mat), arr.ind=TRUE)
  idx <- idx[order(idx[,1], idx[,2]), , drop=FALSE]
  sigma_pairs <- matrix(rownames(equiv_min_dist_mat)[idx], ncol=2)

  coord_array_aperm     <- aperm(coord_array, c(2,1,3))
  coord_synthetic_aperm <- aperm(coord_synthetic, c(2,1,3))
  rates <- c(base_rates, kmethyl=kmethyl, karo=karo, Dx=d_rot, Dy=d_rot, Dz=d_rot)
  perm_rates_by_type <- list("1-1"=c(NA_real_,NA_real_), "1-2"=c(NA_real_,karo=karo), "1-3"=c(NA_real_,kmethyl=kmethyl),
                             "2-2"=c(karo=karo,karo=karo), "2-3"=c(karo=karo,kmethyl=kmethyl), "3-3"=c(kmethyl=kmethyl,kmethyl=kmethyl))
  spec_den_to_relax <- function(spec_den_data, type_name) {
    sigma_value <- spec_den_data$atom_pairs[, "sigma"]; n_relax_rates <- sum(!is.na(sigma_value))
    make_spec_den_relax_data(
      atom_relax_data=list(atom_pairs=spec_den_data$atom_pairs, unit=FALSE,
        relax_data_list=list(sigma=list(value=unname(sigma_value[!is.na(sigma_value)]),
          spec_den_term_array=make_sigma_spec_den_term_array(n_pairs=n_relax_rates, proton_mhz=proton_mhz)))),
      base_rate_mat=base_rate_mat, base_rates=base_rates, perm_rates=perm_rates_by_type[[type_name]])
  }
  spec_den_data_list <- make_spec_den_data(ke_data, sigma_pairs, perm_internal=TRUE)
  sigma_synthetic <- coord_array_to_sigma(coord_synthetic_aperm, c(base_rates, kmethyl=kmethyl, karo=karo, kc=d_rot), spec_den_data_list, proton_mhz)
  for (i in seq_along(spec_den_data_list)) {
    spec_den_data_list[[i]]$atom_pairs <- data.frame(spec_den_data_list[[i]]$atom_pairs, sigma=NA_real_)
    spec_den_data_list[[i]]$atom_pairs[seq_along(sigma_synthetic[[i]]), "sigma"] <- sigma_synthetic[[i]]
  }
  spec_den_relax_data_list <- Map(spec_den_to_relax, spec_den_data_list, names(spec_den_data_list))
  names(spec_den_relax_data_list) <- names(spec_den_data_list)

  for (type_name in names(spec_den_relax_data_list)) {
    sd <- spec_den_relax_data_list[[type_name]]
    ap <- sd$atom_pairs[, 1:2, drop=FALSE]
    relax_df <- spec_den_term_array_list_to_atom_relax_df(sd$relax_data_list)
    M <- nrow(ap); N <- nrow(relax_df)
    pad <- relax_df[rep(NA_integer_, M-N), , drop=FALSE]
    atom_relax <- cbind(ap, rbind(relax_df, pad))
    write.csv (atom_relax,        file.path(out_dir, paste0(type_name, "_atom_relax.csv")),  row.names=FALSE, na="")
    write.table(sd$groupings,     file.path(out_dir, paste0(type_name, "_groupings.csv")),   sep=",", row.names=FALSE, col.names=FALSE)
    write.csv (sd$a_int_coef,     file.path(out_dir, paste0(type_name, "_a_coef.csv")),      row.names=FALSE)
    write.csv (sd$lambda_int_coef,file.path(out_dir, paste0(type_name, "_lambda_coef.csv")))
  }
  relax_energy <- coord_array_to_relax_energy(coord_array_aperm, rates, spec_den_relax_data_list, gradient=FALSE)
  writeLines(format(as.numeric(relax_energy), digits=17), file.path(out_dir, "relax_energy.txt"))
  cat(sprintf("  %-40s %d members  energy = %s  (%d prefixes)\n",
      out_dir, n_members, format(as.numeric(relax_energy), digits=17), length(spec_den_relax_data_list)))
}

cat("Generating GB3 RELAX set (md_double; RELAX needs >=2 ensemble members - 1-member is degenerate):\n")
gen_relax(2, file.path(REPO, "res/relax/md_double"))
cat("done\n")
