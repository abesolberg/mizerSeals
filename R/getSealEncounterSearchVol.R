# Outputs Seal Encounter Rate & Search Volume
# Look into what good default values for f0, h, & q should be with seals (these values are needed to define search volume)
# https://sizespectrum.org/mizer/reference/get_h_default.html
# https://sizespectrum.org/mizer/reference/get_f0_default.html
# https://sizespectrum.org/mizer/reference/setSearchVolume.html?q=sear#null
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
