# Outputs Seal Encounter Rate & Search Volume
# Look into what good default values for f0, h, & q should be with seals (these values are needed to define search volume)
# https://sizespectrum.org/mizer/reference/get_h_default.html
# https://sizespectrum.org/mizer/reference/get_f0_default.html
# https://sizespectrum.org/mizer/reference/setSearchVolume.html?q=sear#null

getSealEncounter <- function(params , search_vol = params@other_params$sealParams$search_vol , n = params@initial_n, n_pp = params@initial_n_pp , ...) {
  sp <- params@other_params$sealParams
  if (is.null(sp)) {
    sp <- list(...)
  }
  n_all <- matrix(NA , nrow = nrow(n)+1 , ncol = ncol(n))
  n_all[1:nrow(n),] <- n
  n_all[nrow(n_all),] <- n_pp[params@w_full %in% params@w]
  
  prey <- array(0 , dim = c(1 , length(sp$w)))
  prey[1:length(params@w)] <- c(sp$interaction_seal , sp$resource_interaction_seal) %*% n_all
  
  prey <- sweep(prey, 2, sp$w * sp$dw, "*")
  avail_energy <- Re(base::t(mvfft(base::t(sp$ft_pred_kernel_e) * mvfft(base::t(prey)), inverse = TRUE)))/length(sp$w)
  avail_energy[avail_energy < 1e-18] <- 0
  encounter <- search_vol * avail_energy
  return(encounter)
}

getSealGamma <- function(params , w , q , h , f0 , ...) {
  sp <- list(w = w , q = q , h = h , f0 = f0)
  search_vol <- sapply(sp$w , function(x) x^sp$q , simplify = T)
  # and setting a power-law prey spectrum
  n <- params@initial_n ; n[] <- 0
  n_pp <- params@initial_n_pp
  n_pp <- params@resource_params$kappa * params@w_full^(-params@resource_params$lambda)
  avail_energy <- getSealEncounter(params , search_vol , n = n , n_pp = n_pp , w = w, ...)[, length(sp$w)] /
    sp$w[length(sp$w)] ^ (2 + sp$q - params@resource_params$lambda)
  # Now set gamma so that this available energy leads to f0
  gamma <- (sp$h/ avail_energy) * (sp$f0 / (1 - sp$f0))
  return(gamma)
}

getSealSearchVol <- function(params , w = params@other_params$sealParams$w , q = params@other_params$sealParams$q , gamma = params@other_params$sealParams$gamma) {
  search_vol <- (w^q)*gamma
  return(search_vol)
}