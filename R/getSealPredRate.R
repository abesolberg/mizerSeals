# Outputs Seal Predation Rate & Prey Mortality Rate (by size)
# Prey Mortality Rate = interaction * pred rate
getSealPredRate <- function(n_seal , w_seal , dw_seal , ft_pred_kernel_p , seal_interaction , params , feeding_level , search_vol) {

  no_w <- length(params@w)
  no_w_full <- length(w_seal)
  Q <- ((1 - feeding_level) * search_vol * c(n_seal)) * dw_seal
  pred_rate <- Re(base::t(mvfft(base::t(ft_pred_kernel_p) *
                                  mvfft(base::t(Q)), inverse = TRUE)))/no_w_full
  pred_rate[pred_rate < 1e-18] <- 0
  pred_mort <- seal_interaction %*% pred_rate[, 1:length(params@w), drop = FALSE]

  return(list(pred_rate = pred_rate , pred_mort = pred_mort))

}
