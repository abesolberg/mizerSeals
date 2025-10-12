# Outputs Seal Encounter Rate & Search Volume
# Look into what good default values for f0, h, & q should be with seals (these values are needed to define search volume)
# https://sizespectrum.org/mizer/reference/get_h_default.html
# https://sizespectrum.org/mizer/reference/get_f0_default.html
# https://sizespectrum.org/mizer/reference/setSearchVolume.html?q=sear#null
getSealEncounterSearchVol <- function(params , w_seal , dw_seal , ft_pred_kernel_e , n , n_pp , seal_interaction , seal_resource_interaction , f0 = .6 , h = 30 , q = .8) {
  
  n_all <- matrix(NA , nrow = nrow(n)+1 , ncol = ncol(n))
  n_all[1:nrow(n),] <- n
  n_all[nrow(n_all),] <- n_pp[params@w_full %in% params@w]
  
  prey <- array(0 , dim = c(1 , length(w_seal)))
  prey[1:length(params@w)] <- c(seal_interaction , seal_resource_interaction) %*% n_all
  
  prey <- sweep(prey, 2, w_seal * dw_seal, "*")
  avail_energy <- Re(base::t(mvfft(base::t(ft_pred_kernel_e) * mvfft(base::t(prey)), inverse = TRUE)))/length(w_seal)
  #avail_energy <- avail_energy[, idx_sp, drop = FALSE]
  avail_energy[avail_energy < 1e-18] <- 0
  gamma <- (h/avail_energy) * (f0/(1 -f0))
  search_vol <- (w_seal^q ) * gamma
  encounter <- search_vol * avail_energy
  
  return(list(encounter = encounter , search_vol = search_vol))
}
