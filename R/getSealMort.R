## Get Seal Mortality -- Exported

getSealMort <- function(params, n = params@initial_n,
                        n_pp = params@initial_n_pp,
                        n_other = params@initial_n_other,
                        t , ...) {
  if(is.null(params@sealParams)) stop('Must add seal parameters to other params.')
  sp <- params@sealParams

  encounter <- getSealEncounter(params , search_vol = sp$search_vol)
  feedingLevel <- getSealFeedingLevel(params , encounter)
  pred_rate <- getSealPredRate(params , n , n_pp , n_other , t , feeding_level = feedingLevel , search_vol = sp$search_vol )

  mu_seal <- getSealMortRate(params , pred_rate)
  return(mu_seal)
}

getSealResourceMort <- function(params, n = params@initial_n,
                                n_pp = params@initial_n_pp,
                                n_other = params@initial_n_other,
                                t ,...) {
  if(is.null(params@sealParams)) stop('Must add seal parameters to other params.')
  sp <- params@sealParams

  encounter <- getSealEncounter(params , search_vol = sp$search_vol)
  feedingLevel <- getSealFeedingLevel(params , encounter)
  pred_rate <- getSealPredRate(params , n , n_pp , n_other , t , feeding_level = feedingLevel , search_vol = sp$search_vol )

  mu_seal <- getSealResourceMortRate(params , pred_rate)

  f <- function (params, n, n_pp, n_other, t, pred_rate, ...) {
    as.vector(params@species_params$interaction_resource %*% pred_rate)
  }
  mort <- f(params, n = n, n_pp = n_pp, n_other = n_other, t = t, pred_rate = mizer::getPredRate(params, n = n, n_pp = n_pp, n_other = n_other, t = t))
  names(mort) <- names(params@initial_n_pp)
  idx <- tail(1:length(mort) , length(params@w))
  mort[idx] <- mort[idx] + c(mu_seal)
  return(mort)
}

