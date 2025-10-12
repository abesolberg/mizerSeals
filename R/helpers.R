# Plot Pred Kernel
plotSealPredKernel <- function(params , idx = NULL) {
  sp <- params@other_params$sealParams
  if(is.null(idx)) idx <- length(sp$w)
  plot(sp$pred_kernel[idx,])
}

getSealMortSim <- function(sim){
  if(!inherits(sim , 'MizerSim')) stop('Not a simulation.')
  idx <- seq(1 , dim(sim@n)[1])
  ret <- sim@n
  ret[] <- 0
  for(i in idx) ret[i,,] <- getSealMort(params = sim@params , n = sim@n[i,,] , n_pp = sim@n_pp[i,],setNames(sim@n_other[i,] , 'seals') , t = i)
  return(ret)
}


getSealDietSim <- function(sim , type = c('diet' , 'total_consumption' , 'diet_by_sp' , 'percent_by_sp')) {
  type <- match.arg(type , several.ok = F)
  if(!inherits(sim , 'MizerSim')) stop('Not a simulation.')
  idx <- seq(1 , dim(sim@n)[1])
  ret <- sapply(idx , function(i) getSealDiet(params = sim@params , n = sim@n[i,,] , n_pp = sim@n_pp[i,] , n_seal = sim@n_other[i,][[1]] , idx_sp = 1:length(params@w) , t = i)[[type]] , simplify = F)
  return(ret)
}

list2array <- function(list){
  dims <- c(dim(list[[1]]), length(list))
  out <- array(NA , dim = dims)
  for(i in 1:dims[3]) out[,,i] <- list[[i]]
  return(out)
}

