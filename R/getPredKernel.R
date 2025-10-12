# Get values needed for pred kernel
# get_phi <- function (ppmr, beta, sigma) {
#   Beta <- log(beta)
#   phi <- exp(-(log(ppmr) - Beta)^2/(2 * sigma^2))
#   return(phi)
# }
# 
# # Get Pred Kernel
# 
# ppmr <- params@w_full/params@w_full[1]
# phis <- get_phi(species_params, ppmr)
# phis[, 1] <- 0
# pred_kernel <- array(0, dim = c(no_sp, no_w, no_w_full), 
#                      dimnames = list(sp = species_params$species, w_pred = signif(params@w, 
#                                                                                   3), w_prey = signif(params@w_full, 3)))
# for (i in 1:no_sp) {
#   min_w_idx <- no_w_full - no_w + 1
#   for (k in seq_len(no_w)) {
#     pred_kernel[i, k, (min_w_idx - 1 + k):1] <- phis[i, 
#                                                      1:(min_w_idx - 1 + k)]
#   }
# }
# return(pred_kernel)