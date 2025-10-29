# Plot Pred Kernel
plotSealPredKernel <- function(params , idx = NULL) {
  sp <- params@seal_params
  if(is.null(idx)) idx <- length(sp$w)
  plot(sp$pred_kernel[idx,])
}

getSealMortSim <- function(sim){
  if(!inherits(sim , 'MizerSim')) stop('Not a simulation.')
  idx <- seq(1 , dim(sim@n)[1])
  ret <- sim@n
  ret[] <- 0
  for(i in idx) ret[i,,] <- getSealMort(params = sim@params , n = sim@n[i,,] , n_pp = sim@n_pp[i,], n_other = setNames(sim@n_other[i,] , 'seals') , t = i)
  return(ret)
}


getSealDietSim <- function(sim , type = c('diet' , 'total_consumption' , 'diet_by_sp' , 'percent_by_sp')) {
  type <- match.arg(type , several.ok = F)
  if(!inherits(sim , 'MizerSim')) stop('Not a simulation.')
  idx <- seq(1 , dim(sim@n)[1])
  ret <- sapply(idx , function(i) getSealDiet(params = sim@params , n = sim@n[i,,] , n_pp = sim@n_pp[i,] , n_other = setNames(sim@n_other[i,] , 'seals') , idx_sp = 1:length(params@w) , t = i)[[type]] , simplify = F)
  return(ret)
}

list2array <- function(list){
  dims <- c(dim(list[[1]]), length(list))
  out <- array(NA , dim = dims)
  for(i in 1:dims[3]) out[,,i] <- list[[i]]
  return(out)
}

getNewVals <- function(params , newvals = list()) {
  for(i in names(newvals)) {
    val <- slot(params , i)
    if (is.list(newvals[[i]])) {
      .newvals <- newvals[[i]]
      for (k in names(.newvals)) {
        val[[k]][] <- .newvals[[k]]
      }
    } else val[] <- newvals[[i]]
    slot(params , i , check = F) <- val
  }
  return(params)
}

get_val_at_t <- function(val , t, ...) {
  if (is.vector(val)) {
    arr <- matrix(val , ncol = 1)
  } else arr <- val
  if (is.null(rownames(arr))) 
    rownames(arr) <- 1:nrow(arr)
  if (!floor(t) %in% as.numeric(rownames(arr))) 
    t <- t%%as.numeric(rownames(arr)[nrow(arr)]) + as.numeric(rownames(arr)[1])
  index_vec <- as.numeric(rownames(arr)) - t
  index <- which(index_vec == max(index_vec[index_vec <= 0]))
  var_at_t <- arr[index, , drop = F]
  return(var_at_t[1,])
}


