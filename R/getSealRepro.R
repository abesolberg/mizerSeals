## Get Seal Reproduction -- Exported
getSealRepro <- function(params, n = params@initial_n,
                        n_pp = params@initial_n_pp ,
                        n_other = params@initial_n_other,
                        t , dt , ...) {
  sp <- params@other_params$sealParams
  if(is.null(sp)) stop('Must add seal parameters to other params.')
  if(!sp$dynamicSeals) return(getSealN(params , t = t))

  # Set Prey & Bm0
  prey <- getSealDiet(params , n = n , n_seal = n_other$seals)$total_consumption
  Bm0 <- sum(n_other$seals)
  H <- getSealH(params , t)*dt # Change harvest level to fit with time step

  list <- append(c(prey = prey , Bm0 = Bm0 , H = H) , sp)
  ret <- do.call(setBioMod , list)
  return(ret)
}
