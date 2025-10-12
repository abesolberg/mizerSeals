# Compile Seal Predation Kernel
sealPredKernel <- function(params , w , beta = 1000 , sigma = 1){

  no_w_full <- length(w)
  ft_pred_kernel_e <- array(NA, dim = c(1, no_w_full), dimnames = list(sp = c('harp seals') , k = 1:no_w_full))
  ft_pred_kernel_p <- ft_pred_kernel_e

  ppmr <- w/w[1]
  phi <- get_phi(ppmr , beta , sigma)
  phi[1] <- 0

  pred_kernel <- array(
    0,
    dim = c(1, length(w), length(w)),
    dimnames = list(
      sp = 'Seal',
      w_pred = signif(w, 3),
      w_prey = signif(w, 3)
    )
  )

  for (k in 1:dim(pred_kernel)[2]) {
    pred_kernel[1, k, k:1] <- phi[1:k]
  }

  ft_pred_kernel_e[1,] <- fft(phi)
  ri <- min(max(which(phi > 0)), no_w_full - 1)
  phi_p <- rep(0, no_w_full)
  phi_p[(no_w_full - ri + 1):no_w_full] <- phi[(ri + 1):2]
  ft_pred_kernel_p[1, ] <- fft(phi_p)

  return(
    list(
      ft_pred_kernel_p = ft_pred_kernel_p ,
      ft_pred_kernel_e = ft_pred_kernel_e ,
      pred_kernel = pred_kernel[1,,,drop=T]
    ))

}
