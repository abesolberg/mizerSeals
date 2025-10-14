# Outputs Seal Encounter Rate & Search Volume
# Look into what good default values for f0, h, & q should be with seals (these values are needed to define search volume)
# https://sizespectrum.org/mizer/reference/get_h_default.html
# https://sizespectrum.org/mizer/reference/get_f0_default.html
# https://sizespectrum.org/mizer/reference/setSearchVolume.html?q=sear#null
getSealGamma <- function(params) {
  params@other_params$sealParams$gamma <- 1
  # params <- setSearchVolume(params)
  sweep(outer(params@species_params[["q"]], params@w,
              function(x, y) y^x), 1, params@species_params$gamma,
        "*")

  # and setting a power-law prey spectrum
  params@initial_n[] <- 0
  if (defaults_edition() < 2) {
    # See issue #238
    params@species_params$interaction_resource <- 1
  }
  params@initial_n_pp[] <- params@resource_params$kappa *
    params@w_full^(-params@resource_params$lambda)
  avail_energy <- getEncounter(params)[, length(params@w)] /
    params@w[length(params@w)] ^
    (2 + params@species_params[["q"]] - params@resource_params$lambda)
  # Now set gamma so that this available energy leads to f0
  gamma_default <- (species_params[["h"]] / avail_energy) *
    (species_params$f0 / (1 - species_params$f0))
}

getSealSearchVol <- function(params , n = params@initial_n, n_pp = params@initial_n_pp , ...) {

  sp <- params@other_params$sealParams
  n_all <- matrix(NA , nrow = nrow(n)+1 , ncol = ncol(n))
  n_all[1:nrow(n),] <- n
  n_all[nrow(n_all),] <- n_pp[params@w_full %in% params@w]

  prey <- array(0 , dim = c(1 , length(sp$w)))
  prey[1:length(params@w)] <- c(sp$interaction_seal , sp$resource_interaction_seal) %*% n_all

  prey <- sweep(prey, 2, sp$w * sp$dw, "*")
  avail_energy <- Re(base::t(mvfft(base::t(sp$ft_pred_kernel_e) * mvfft(base::t(prey)), inverse = TRUE)))/length(sp$w)
  #avail_energy <- avail_energy[, idx_sp, drop = FALSE]
  avail_energy[avail_energy < 1e-18] <- 0
  gamma <- (sp$h/avail_energy) * (sp$f0/(1 - sp$f0))
  search_vol <- (sp$w^sp$q ) * gamma
  return(search_vol)
}

getSealEncounter <- function(params , search_volume , ...) {
  sp <- params@other_params$sealParams
  avail_energy <- sp$h/(search_volume/(sp$w^sp$q)/(sp$f0/(1 - sp$f0)))
  encounter <- search_volume * avail_energy
  return(encounter)
}
