## Get Seal Mortality -- Exported
getSealMort <- function(params, n = params@initial_n,
                        n_pp = params@initial_n_pp,
                        n_other = params@initial_n_other,
                        t , ...) {

  if(is.null(params@other_params$sealParams)) stop('Must add seal parameters to other params.')
  sp <- params@other_params$sealParams

  encounterSearchVol <- getSealEncounterSearchVol(
    params ,
    w_seal = sp$w ,
    dw_seal = sp$dw ,
    ft_pred_kernel_e = sp$ft_pred_kernel_e ,
    n = n,
    n_pp = n_pp , 
    seal_interaction = sp$interaction_seal ,
    seal_resource_interaction = sp$resource_interaction_seal ,
    f0 = sp$f0  ,
    h = sp$h  ,
    q = sp$q
  )

  feedingLevel <- getSealFeedingLevel(w_seal = sp$w , seal_encounter = encounterSearchVol$encounter , h = sp$h , n = sp$n)
  sealRates <- getSealPredRate(
    params = params ,
    n_seal = n_other$seals ,
    w_seal = sp$w ,
    dw_seal = sp$dw ,
    ft_pred_kernel_p = sp$ft_pred_kernel_p  ,
    seal_interaction = sp$interaction_seal ,
    feeding_level = feedingLevel ,
    search_vol = encounterSearchVol$search_vol
  )

  mu_s <- sealRates$pred_mort
  return(mu_s)
}

