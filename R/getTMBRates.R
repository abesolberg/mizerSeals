# Get TMB Rates
getTMBPredRate <- function (params, n, n_pp, n_other, t, feeding_level, ...) {
  no_sp <- dim(params@interaction)[1]
  no_w <- length(params@w)
  no_w_full <- length(params@w_full)
  if (!is.null(comment(params@pred_kernel))) {
    n_total_in_size_bins <- sweep(n, 2, params@dw, "*", check.margin = FALSE)
    pred_rate <- sweep(params@pred_kernel, c(1, 2), (1 - 
                                                       feeding_level) * params@search_vol * n_total_in_size_bins, 
                       "*", check.margin = FALSE)
    pred_rate <- colSums(aperm(pred_rate, c(2, 1, 3)), dims = 1)
    return(pred_rate)
  }
  ft_pred_kernel_p <- params@ft_pred_kernel_p_real+1i*params@ft_pred_kernel_p_imag
  idx_sp <- (no_w_full - no_w + 1):no_w_full
  Q <- matrix(0, nrow = no_sp, ncol = no_w_full)
  Q[, idx_sp] <- sweep((1 - feeding_level) * params@search_vol * 
                         n, 2, params@dw, "*")
  pred_rate <- Re(base::t(mvfft(base::t(ft_pred_kernel_p) * 
                                  mvfft(base::t(Q)), inverse = TRUE)))/no_w_full
  pred_rate[pred_rate < 1e-18] <- 0
  return(pred_rate * params@ft_mask)
}

getTMBEncounterRate <- function(params , n, n_pp, n_other, t, ... ){
  idx_sp <- (length(params@w_full) - length(params@w) + 1):length(params@w_full)
  if (!is.null(comment(params@pred_kernel))) {
    n_eff_prey <- sweep(params@interaction %*% n, 2, params@w * 
                          params@dw, "*", check.margin = FALSE)
    phi_prey_species <- rowSums(sweep(params@pred_kernel[, 
                                                         , idx_sp, drop = FALSE], c(1, 3), n_eff_prey, "*", 
                                      check.margin = FALSE), dims = 2)
    phi_prey_background <- params@species_params$interaction_resource * 
      rowSums(sweep(params@pred_kernel, 3, params@dw_full * 
                      params@w_full * n_pp, "*", check.margin = FALSE), 
              dims = 2)
    encounter <- params@search_vol * (phi_prey_species + 
                                        phi_prey_background)
  }
  else {
    ft_pred_kernel_e <- params@ft_pred_kernel_e_real+1i*params@ft_pred_kernel_e_imag
    
    prey <- outer(params@species_params$interaction_resource, 
                  n_pp)
    prey[, idx_sp] <- prey[, idx_sp] + params@interaction %*% 
      n
    prey <- sweep(prey, 2, params@w_full * params@dw_full, 
                  "*")
    avail_energy <- Re(base::t(mvfft(base::t(ft_pred_kernel_e) * 
                                       mvfft(base::t(prey)), inverse = TRUE)))/length(params@w_full)
    avail_energy <- avail_energy[, idx_sp, drop = FALSE]
    avail_energy[avail_energy < 1e-18] <- 0
    encounter <- params@search_vol * avail_energy
  }
  for (i in seq_along(params@other_encounter)) {
    encounter <- encounter + do.call(params@other_encounter[[i]], 
                                     list(params = params, n = n, n_pp = n_pp, n_other = n_other, 
                                          component = names(params@other_encounter)[[i]], 
                                          ...))
  }
  return(encounter + params@ext_encounter)
}
