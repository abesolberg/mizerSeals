# Outputs Seal Predation Rate & Prey Mortality Rate (by size)
# Prey Mortality Rate = interaction * pred rate
getSealPredRate <- function(params ,
                            n = params@initial_n,
                            n_pp = params@initial_n_pp,
                            n_other = params@initial_n_other,
                            t , feeding_level , search_vol , ... ) {

  sp <- params@seal_params

  ft_pred_kernel_p <- sp$ft_pred_kernel_p_real+1i*sp$ft_pred_kernel_p_imag
  
  no_w <- length(params@w)
  no_w_full <- length(sp$w)
  Q <- ((1 - feeding_level) * search_vol * c(n_other$seals)) * sp$dw
  pred_rate <- Re(base::t(mvfft(base::t(ft_pred_kernel_p) *
                                  mvfft(base::t(Q)), inverse = TRUE)))/no_w_full
  pred_rate[pred_rate < 1e-18] <- 0
  return(pred_rate)
}

getSealMortRate <- function(params , pred_rate , ... ) {
  sp <- params@seal_params
  pred_mort <- sp$interaction_seal %*% pred_rate[, 1:length(params@w), drop = FALSE]
  return(pred_mort)
}

getSealResourceMortRate <- function(params , pred_rate , ... ) {
  sp <- params@seal_params
  pred_mort <- sp$resource_interaction_seal %*% pred_rate[,1:length(params@w), drop = FALSE]
  return(pred_mort)
}
