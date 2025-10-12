# outputs matrix of prey (amount) consumed by seals at size
getSealDiet <- function(params , n = params@initial_n , n_pp = params@initial_n_pp , n_seal = params@initial_n_other[['seals']], idx_sp = 1:length(params@w) , t = NULL , ...) {

  seal_params <- params@other_params$sealParams
  #n_seal <- sealParams$initialSealN[t,]

  encounterSearchVol <- getSealEncounterSearchVol(
    params ,
    w_seal = seal_params$w ,
    dw = seal_params$dw ,
    ft_pred_kernel_e = seal_params$ft_pred_kernel_e ,
    n = n,
    n_pp = n_pp , 
    seal_params$interaction_seal ,
    seal_resource_interaction = seal_params$resource_interaction_seal ,
    f0 = seal_params$f0  ,
    h = seal_params$h  ,
    q = seal_params$q
  )

  feedingLevel <- getSealFeedingLevel(
    w_seal = seal_params$w ,
    seal_encounter = encounterSearchVol$encounter ,
    h = seal_params$h ,
    n = seal_params$n
  )

  seal_interaction <- c(seal_params$interaction_seal , seal_params$resource_interaction_seal)
  search_vol <- encounterSearchVol$search_vol

  no_sp <- nrow(params@species_params) + 1
  no_w_full <- length(seal_params$w)
  
  n_all <- matrix(NA , nrow = nrow(n)+1 , ncol = ncol(n))
  n_all[1:nrow(n),] <- n
  n_all[nrow(n_all),] <- n_pp[params@w_full %in% params@w]
  
  prey <- matrix(0, nrow = no_sp , ncol = no_w_full)
  prey[1:no_sp, idx_sp] <- sweep(n_all, 2, params@w * params@dw, "*")

  ft <- array(0 , dim = c(1,  no_w_full, no_sp))
  for(i in 1:dim(ft)[3]) ft[,,i] <- seal_params$ft_pred_kernel_e*mvfft(t(prey))[,i]

  ft <- matrix(aperm(ft , c(2, 1, 3)) , nrow = no_w_full)

  ae <- array(Re(mvfft(ft, inverse = TRUE)/no_w_full),
              dim = c(no_w_full, 1, no_sp))
  ae <- aperm(ae, c(2, 1, 3))
  ae[ae < 1e-18] <- 0
  ae <- ae[, , 1:no_sp , drop = F]

  diet <- t(ae[1,,,drop = T])

  for(i in 1:nrow(diet)) {
    diet[i,] = diet[i,]*encounterSearchVol$search_vol*seal_interaction[i]
  }
  diet <- diet * matrix(data = (1 - feedingLevel[1,]) , nrow = no_sp , ncol = no_w_full , byrow = T)
  diet[,n_seal <= 0] <- 0

  return(
    list(
      diet = diet ,
      total_consumption = sum(sweep(diet , 2 , n_seal*seal_params$dw , "*")) ,
      diet_by_sp = rowSums(sweep(diet , 2 , n_seal*seal_params$dw , "*")) ,
      percent_by_sp = rowSums(sweep(diet , 2 , n_seal*seal_params$dw , "*"))/sum(sweep(diet , 2 , n_seal*seal_params$dw , "*"))
    )
  )
}
