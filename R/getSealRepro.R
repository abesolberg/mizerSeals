## Get Seal Reproduction -- Exported
getSealRepro <- function(params, n = params@initial_n,
                        n_pp = params@initial_n_pp ,
                        n_other = params@initial_n_other,
                        t , dt , ...) {
  sp <- params@seal_params
  if(is.null(sp)) stop('Must add seal parameters to other params.')
  if(!sp$dynamicSeals) return(getSealN(params , t = t))

  # Set Prey & Bm0
  prey <- getSealDiet(params , n = n , n_seal = n_other$seals)$total_consumption
  H <- getSealH(params , t)*dt # Change harvest level to fit with time step
  alpha <- getCoeff(params , par = 'alpha' , where = 'seal_params' , t) # Random Effect for Seals
  list <- append(list(N = n_other$seals , prey = prey , H = H , .alpha = alpha) , sp)
  ret <- do.call(setBioMod , list)
  return(ret)
}
