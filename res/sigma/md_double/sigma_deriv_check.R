source("~/Documents/Wesleyan/R/ke/R/ke.R")
source("~/Documents/Wesleyan/R/ke/R/ke_prep.R")

ensemble <- read_ensemble(c("GB3_27_10us.pdb","GB3_27_10us.pdb"), proton_only = TRUE)

# coordinates against which derivatives will be calculated
coord_array <- ensemble[,,1:2,drop=FALSE]

# coordinates for creating synthetic cross-relaxation rates
coord_synthetic <- ensemble[,,1:2,drop=FALSE]

# ensemble member interconversion eigenvalue
base_rates <- c(kens = 1/2e-9)

# molecular tumbling rate
kc <- 1/4e-9

# proton field strength in MHz
proton_mhz <- 700

# transition rate matrix for ensemble members
base_rate_mat <- rate_mat_simple(base_rates, dimnames(coord_synthetic)[[3]])

# data structure to hold data related to kinetic ensemble calculation
ke_data <- make_ke_data(coord_synthetic, base_rate_mat, base_rates, kc, proton_mhz=proton_mhz, mix_times=numeric())

# calculate an array with interatomic distances for three ensemble members
dist_array <- simplify2array(apply(coord_synthetic, 3, function(x) as.matrix(dist(t(x))), simplify=FALSE))

# calculate matrix with minimum distances across three ensemble members
min_dist_mat <- apply(dist_array, 1:2, min)

# create empty matrix with minimum distances between equivalent atom groups
equiv_min_dist_mat <- matrix(NA_real_, nrow=length(ke_data$equiv_list), ncol=length(ke_data$equiv_list), dimnames=list(names(ke_data$equiv_list), names(ke_data$equiv_list)))

# create a named vector mapping for equivalent atom groups with just a single atom
equiv_single <- unlist(ke_data$equiv_list[sapply(ke_data$equiv_list, length) == 1])

# copy distances from equivalent groups with just a single atom
equiv_min_dist_mat[names(equiv_single),names(equiv_single)] <- min_dist_mat[equiv_single,equiv_single]

# find equivalent groups with multiple atoms
equiv_multi <- ke_data$equiv_list[sapply(ke_data$equiv_list, length) > 1]

# create pairs of all possible combinations of all groups with multi groups
equiv_multi_pairs <- expand.grid(names(ke_data$equiv_list), names(equiv_multi), stringsAsFactors=FALSE)

# compute minimum distances for multi group pairs
for (i in seq_len(nrow(equiv_multi_pairs))) {
	min_dist <- min(min_dist_mat[ke_data$equiv_list[[equiv_multi_pairs[i,1]]],ke_data$equiv_list[[equiv_multi_pairs[i,2]]]])
	equiv_min_dist_mat[equiv_multi_pairs[i,1],equiv_multi_pairs[i,2]] <- equiv_min_dist_mat[equiv_multi_pairs[i,2],equiv_multi_pairs[i,1]] <- min_dist
}

# select pairs of equivalent atom groups whose minimum distances are less than 5 angstroms
idx <- which(equiv_min_dist_mat < 5 & upper.tri(equiv_min_dist_mat), arr.ind=TRUE)
idx <- idx[order(idx[,1], idx[,2]),]
sigma_pairs <- matrix(rownames(equiv_min_dist_mat)[idx], ncol=2)

# create data structures for calculating spectral density functions
spec_den_data_list <- make_spec_den_data(ke_data, sigma_pairs, perm_internal=TRUE)

# calculate synthetic cross relaxation rates
sigma_synthetic <- coord_array_to_sigma(aperm(coord_synthetic, c(2,1,3)), ke_data$rates, spec_den_data_list, proton_mhz)

# insert synthetic cross relaxation rates into spec_den_data_list
for (i in seq_along(spec_den_data_list)) {
	spec_den_data_list[[i]]$atom_pairs <- data.frame(spec_den_data_list[[i]]$atom_pairs, sigma=NA_real_)
	spec_den_data_list[[i]]$atom_pairs[seq_along(sigma_synthetic[[i]]), "sigma"] <- sigma_synthetic[[i]]
}

# write out input data if input_dir is set to something other than NULL
input_dir <- "."
if (!is.null(input_dir)) {
	for (i in seq_along(spec_den_data_list)) {
	
		type_name <- names(spec_den_data_list)[i]
		write.csv(spec_den_data_list[[i]][["atom_pairs"]], file.path(input_dir, paste0(type_name, "_atom_pairs.csv")), row.names=FALSE, na="")
		write.table(spec_den_data_list[[i]][["groupings"]], file.path(input_dir, paste0(type_name, "_groupings.csv")), sep=",", row.names=FALSE, col.names=FALSE)
		write.csv(spec_den_data_list[[i]][["a_coef"]], file.path(input_dir, paste0(type_name, "_a_coef.csv")), row.names=FALSE)
		write.csv(spec_den_data_list[[i]][["lambda_coef"]], file.path(input_dir, paste0(type_name, "_lambda_coef.csv")))
	}
}
stop()

sigma_energy <- coord_array_to_sigma_energy(aperm(coord_array, c(2,1,3)), ke_data$rates, spec_den_data_list, proton_mhz, gradient=TRUE)

if (!"sigma_deriv_check" %in% ls()) {
	sigma_deriv_check <- deriv_check(coord_array_to_sigma_energy, aperm(coord_array, c(2,1,3)), dv=1e-6, vdims=1:3, gdims=1:3, rates=ke_data$rates, spec_den_data_list=spec_den_data_list, proton_mhz=700)
}

gradient <- attr(sigma_energy, "gradient")
gradient_fd <- gradient - sigma_deriv_check

par(mfrow=c(1, 1))

plot(gradient_fd, gradient, xlab="Finite Difference Gradient", ylab="Analytical Gradient", )
abline(0, 1, col="red")

spec_den_data_ind_list <- make_spec_den_data(ke_data, sigma_pairs, perm_internal=FALSE)

spec_den_data_int_list <- make_spec_den_data(ke_data, matrix(character(), ncol=2), perm_internal=TRUE)
names(spec_den_data_int_list) <- c("2 Internal", "3 Internal")

spec_den_data_ind_list <- c(spec_den_data_ind_list, spec_den_data_int_list)

# calculate synthetic cross relaxation rates
sigma_synthetic_ind <- coord_array_to_sigma(aperm(coord_synthetic, c(2,1,3)), ke_data$rates, spec_den_data_ind_list, proton_mhz)

# insert synthetic cross relaxation rates into spec_den_data_list
for (i in seq_along(spec_den_data_ind_list)) {
	spec_den_data_ind_list[[i]]$atom_pairs <- data.frame(spec_den_data_ind_list[[i]]$atom_pairs, sigma=NA_real_)
	spec_den_data_ind_list[[i]]$atom_pairs[seq_along(sigma_synthetic_ind[[i]]), "sigma"] <- sigma_synthetic_ind[[i]]
}

sigma_energy_ind <- coord_array_to_sigma_energy(aperm(coord_array, c(2,1,3)), ke_data$rates, spec_den_data_ind_list, proton_mhz, gradient=TRUE)

if (! "gradient_list" %in% ls()) {

	gradient_list <- lapply(seq_along(spec_den_data_ind_list), function(i) {

		sigma_energy <- coord_array_to_sigma_energy(aperm(coord_array, c(2,1,3)), ke_data$rates, spec_den_data_ind_list[i], proton_mhz, gradient=TRUE)
	
		sigma_deriv_check <- deriv_check(coord_array_to_sigma_energy, aperm(coord_array, c(2,1,3)), dv=1e-6, vdims=1:3, gdims=1:3, rates=ke_data$rates, spec_den_data_list=spec_den_data_ind_list[i], proton_mhz=700)

		gradient <- attr(sigma_energy, "gradient")
		gradient_fd <- gradient - sigma_deriv_check
	
		list(analytical=gradient, finite_difference=gradient_fd)
	})
	names(gradient_list) <- names(spec_den_data_ind_list)
}

par(mfrow=c(2, 4))

for (type in names(gradient_list)) {

	gradient_fd <- gradient_list[[type]]$finite_difference
	gradient <- gradient_list[[type]]$analytical

	slope <- lsfit(gradient_fd, gradient, intercept=FALSE)$coef
	
	plot(gradient_fd, gradient, asp=1, xlab="Finite Difference Gradient", ylab="Analytical Gradient", main=type)
	
	abline(0, slope, col="red")
	abline(0, 1, col="blue")
	
	legend("topleft", legend=c("y = x", sprintf("y = %0.6fx", slope)), lwd=1, col=c("blue", "red"), bty="n")
}
