# List of Rates

RTMBRates <- list(
  Rates = "mizerRatesRTMB" ,
  Encounter = "getEncounterRTMB",
  FeedingLevel = "getFeedinLevelRTMB",
  EReproAndGrowth = "getEReproAndGrowthRTMB",
  ERepro = "getEReproRTMB",
  EGrowth = "getEGrowthRTMB",
  PredRate = "getPredRateRTMB",
  PredMort = "getPredMortRTMB",
  FMort = "getFMortRTMB",
  SealEncounter = "getSealEncounterRTMB",
  SealFeedingLevel = "getSealFeedingLevelRTMB",
  SealPredRate = "getSealPredRateRTMB",
  SealMort = "getSealMortRTMB",
  Mort = "getMortRTMB",
  RDI = "getRDIRTMB",
  RDD = "getRDDRTMB",
  SealResourceMortRate = "getSealResourceMortRateRTMB",
  ResourceMortRate = "getResourceMortRateRTMB",
  ResourceMort = "getResourceMortRTMB",
  SealDiet = "getSealDietRTMB",
  SealRepro = "getSealReproRTMB",
  ResourceDynamics = "getResourceDynamicsRTMB" ,
  TrawlCatch = "getTrawlCatchRTMB" ,
  Yield = "getYieldRTMB"
)

# Mizer Rates Function

mizerRatesRTMB <- function(n , n_pp , n_other , dt , effort , effort_trawl , rdd , z , alpha_seal , H , dat) {
  list2env(dat , environment())
  r <- list()
  r$encounter <- rates_fns$Encounter(n = n , n_pp = n_pp , n_other = n_other)
  r$feeding_level <- rates_fns$FeedingLevel(r$encounter)
  r$e <- rates_fns$EReproAndGrowth(encounter = r$encounter, feeding_level = r$feeding_level)
  r$e_repro <- rates_fns$ERepro(e = r$e)
  r$e_growth <- rates_fns$EGrowth(e_repro = r$e_repro, e = r$e)
  r$pred_rate <- rates_fns$PredRate(n = n, n_pp = n_pp, n_other = n_other, feeding_level = r$feeding_level)
  r$pred_mort <- rates_fns$PredMort(pred_rate = r$pred_rate)
  r$f_mort <- rates_fns$FMort()
  r$encounter_seal <- rates_fns$SealEncounter(n = n , n_pp = n_pp)
  r$feeding_level_seal <- rates_fns$SealFeedingLevel(r$encounter_seal)
  r$pred_rate_seal <- rates_fns$SealPredRate(n = n , n_pp = n_pp , n_other = n_other , feeding_level = r$feeding_level_seal)
  r$seal_mort <- rates_fns$SealMort(pred_rate = r$pred_rate_seal)
  r$mort <- rates_fns$Mort(pred_mort = r$pred_mort , f_mort = r$f_mort , seal_mort =  r$seal_mort , z = z)
  r$rdi <- rates_fns$RDI(e_repro = r$e_repro , n = n)
  r$rdd <- rates_fns$RDD(rdd)
  r$resource_mort_rate_seal <- rates_fns$SealResourceMortRate(pred_rate = r$pred_rate_seal)
  r$resource_mort_rate <- rates_fns$ResourceMortRate(pred_rate = r$pred_rate)
  r$resource_mort <- rates_fns$ResourceMort(resource_mort_rate = r$resource_mort_rate , resource_mort_rate_seal = r$resource_mort_rate_seal)
  r$seal_diet <- rates_fns$SealDiet(n = n , n_pp = n_pp , n_other = n_other , feeding_level = r$feeding_level_seal)
  r$seal_repro <- rates_fns$SealRepro(n_other = n_other , prey = r$seal_diet  , alpha_seal = alpha_seal)
  r$n_pp_new <- rates_fns$ResourceDynamics(n_pp = n_pp , resource_mort = r$resource_mort)
  r$t_catch <- rates_fns$TrawlCatch(n = n , effort_trawl = effort_trawl)
  r$yield <- rates_fns$Yield(f_mort = r$f_mort , n = n)
  return(r)
}

rates_vars <- c(
  "alpha","catchability","catchability_trawl","cc_pp","dt","dw","dw_full","dw_seal","e_seal","erepro","ext_encounter",
  "ft_mask","ft_pred_kernel_e_imag","ft_pred_kernel_e_imag_seal","ft_pred_kernel_e_real","ft_pred_kernel_e_real_seal",
  "ft_pred_kernel_p_imag","ft_pred_kernel_p_imag_seal","ft_pred_kernel_p_real","ft_pred_kernel_p_real_seal","h_seal","idx_seal",
  "intake_max","interaction","interaction_resource","interaction_seal","Jmax_seal","M_seal","metab","Msrr_seal","mu_b","n_seal",
  "psi","R_max","resource_interaction_seal","resource_seal","rr_pp","sc","sealHarvest_seal","search_vol","search_vol_seal",
  "selectivity","selectivity_trawl","w","w_full","w_min_idx","w_seal"
)

f <- function() {
  args <- mget(ls(envir = parent.frame()) , envir = parent.frame())
  if(any(sapply(args , is.null))) {
    missing <- mget(names(args)[sapply(args , is.null)] , env = parent.frame(2) , mode = 'any')
    args <- append(args[!sapply(args , is.null)] , missing)
    list2env(args , envir = parent.frame())
  }
}

# Encounter
getEncounterRTMB <- function(n , n_pp , n_other , w = NULL , w_full = NULL , dw_full = NULL ,
                             ft_pred_kernel_e_real = NULL , ft_pred_kernel_e_imag = NULL ,
                             interaction = NULL , interaction_resource = NULL ,  search_vol = NULL) {
  f()
  idx_sp <- (length(w_full) - length(w) + 1):length(w_full)
  ft_pred_kernel_e <- as.complex(ft_pred_kernel_e_real) + (0+1i) *
    as.complex(ft_pred_kernel_e_imag)
  prey <- outer(interaction_resource, n_pp)
  prey[, idx_sp] <- prey[, idx_sp] + interaction %*% n
  prey <- sweep(prey, 2, w_full * dw_full, "*")
  m <- AD(matrix(0 , nrow = nrow(n) , ncol = length(w_full)))
  for(i in 1:nrow(n)) {
    avail_energy <- Re((RTMB::fft((
      as.complex(ft_pred_kernel_e_real[i, ]) + (0 + 1i) * as.complex(ft_pred_kernel_e_imag[i, ])
    ) *
      RTMB::fft(base::t(prey[i, ])) / (length(w_full)),
    inverse = TRUE
    )))
    m[i, ] <- avail_energy
  }
  avail_energy <- as.matrix(m[, idx_sp, drop = FALSE])
  #avail_energy[avail_energy < 1e-18] <- 0
  encounter <- search_vol * avail_energy
  return(encounter)
}

# Feeding Level
getFeedinLevelRTMB <- function(encounter, intake_max = NULL) {
  f()
  return(encounter/(encounter + intake_max))
}

## E Repro & Growth ----

# EreproAndGrowth
getEReproAndGrowthRTMB <- function(encounter , feeding_level , metab = NULL , alpha = NULL) {
  f()
  sweep((1 - feeding_level) * encounter, 1, alpha,
        "*", check.margin = FALSE) - metab
}

# Erepro
getEReproRTMB <- function(e , psi = NULL) {
  f()
  #e[e<0] <- 0
  psi * e
}

# EGrowth
getEGrowthRTMB <- function(e , e_repro) {
  f()
  #e[e < 0] <- 0
  e - e_repro
}

## Predation ----

getPredRateRTMB <- function (n , n_pp , n_other ,feeding_level ,
                             w = NULL, w_full = NULL , dw = NULL, interaction = NULL,
                             ft_pred_kernel_p_real = NULL , ft_pred_kernel_p_imag = NULL,
                             search_vol = NULL , ft_mask = NULL) {
  f()
  no_sp <- dim(interaction)[1]
  no_w <- length(w)
  no_w_full <- length(w_full)

  # ft_pred_kernel_p <- as.complex(ft_pred_kernel_p_real) + (0+1i) *
  #   as.complex(ft_pred_kernel_p_imag)
  idx_sp <- (no_w_full - no_w + 1):no_w_full
  Q <- matrix(0, nrow = no_sp, ncol = no_w_full)
  Q[, idx_sp] <- sweep((1 - feeding_level) * search_vol *
                         n, 2, dw, "*")

  m <- AD(matrix(0 , nrow = nrow(n) , ncol = length(w_full)))
  for(i in 1:nrow(n)) {
    pred_rate <- Re((RTMB::fft((
      as.complex(ft_pred_kernel_p_real[i, ]) + (0 + 1i) * as.complex(ft_pred_kernel_p_imag[i, ])
    ) *
      RTMB::fft(base::t(Q[i, ])) / (length(w_full)),
    inverse = TRUE
    )))
    m[i, ] <- pred_rate
  }
  pred_rate <- as.matrix(m)

  # pred_rate <- Re(base::t(mvfft(base::t(ft_pred_kernel_p) *
  #                                 mvfft(base::t(Q)), inverse = TRUE)))/no_w_full
  # pred_rate[pred_rate < 1e-18] <- 0
  return(pred_rate * ft_mask)
}

getPredMortRTMB <- function(pred_rate , w = NULL , w_full = NULL, interaction = NULL) {
  f()
  idx_sp <- (length(w_full) - length(w) + 1):length(w_full)
  return((base::t(interaction) %*% pred_rate[, idx_sp, drop = FALSE]))
}

## Fishing Mortality ----

getFMortRTMB <- function(effort = NULL , catchability = NULL, selectivity = NULL) {
  f()
  out <- selectivity
  out[] <- effort * c(catchability) * c(selectivity)
  return(colSums(out))
}

getTrawlCatchRTMB <- function(n , effort_trawl , w = NULL , dw = NULL , catchability_trawl = NULL, selectivity_trawl = NULL) {
  f()
  out <- selectivity_trawl
  out[] <- effort_trawl * c(catchability_trawl) * c(selectivity_trawl)
  biomass <- sweep(n , 2 , w*dw , '*')
  out <- apply(sweep(out , c(2, 3), biomass, "*"), c(1, 2), sum)
  return(out)
}

getYieldRTMB <- function(f_mort , n , w = NULL , dw = NULL , dt = NULL) {
  f()
  biomass <- sweep(n , 2 , w*dw , '*')
  out <- apply(f_mort * biomass, 1 , sum)*dt
}


## Seal Predation Rates ----

getSealEncounterRTMB <- function(n , n_pp , search_vol_seal = NULL ,
                                 w_full = NULL , w = NULL, w_seal = NULL , dw_seal = NULL,
                                 interaction_seal = NULL , resource_interaction_seal = NULL ,
                                 ft_pred_kernel_e_real_seal = NULL , ft_pred_kernel_e_imag_seal = NULL) {
  f()
  n_all <- matrix(0, nrow = nrow(n) + 1, ncol = ncol(n))
  n_all[1:nrow(n), ] <- n
  n_all[nrow(n_all), ] <- n_pp[w_full %in% w]
  prey <- array(0, dim = c(1, length(w_seal)))
  prey[1:length(w)] <- c(interaction_seal, resource_interaction_seal) %*% n_all
  #ft_pred_kernel_e <- as.complex(ft_pred_kernel_e_real_seal) + (0+1i) * as.complex(ft_pred_kernel_e_imag_seal)
  prey <- sweep(prey, 2, w_seal * dw_seal, "*")
  avail_energy <- Re((
    RTMB::fft(
      (as.complex(ft_pred_kernel_e_real_seal) + (0 + 1i) * as.complex(ft_pred_kernel_e_imag_seal)) *
        RTMB::fft((prey)) / length(w_seal),
      inverse = TRUE
    )
  ))
  encounter <- c(search_vol_seal * avail_energy)
  return(encounter)
}

getSealFeedingLevelRTMB <- function(encounter , w_seal = NULL , h_seal = NULL, n_seal = NULL) {
  f()
  intake_max <- h_seal * w_seal^n_seal
  return(encounter/(encounter + intake_max))
}

getSealPredRateRTMB <- function(n , n_pp , n_other , feeding_level ,
                                w = NULL , w_seal = NULL, dw_seal = NULL,
                                ft_pred_kernel_p_real_seal = NULL , ft_pred_kernel_p_imag_seal = NULL,
                                search_vol_seal = NULL) {
  f()
  #ft_pred_kernel_p <- as.complex(ft_pred_kernel_p_real_seal) + (0+1i) * as.complex(ft_pred_kernel_p_imag_seal)
  no_w <- length(w)
  no_w_full <- length(w_seal)
  Q <- ((1 - feeding_level) * search_vol_seal * c(n_other)) *
    dw_seal
  pred_rate <- Re(base::t(RTMB::fft((as.complex(ft_pred_kernel_p_real_seal) + (0+1i) * as.complex(ft_pred_kernel_p_imag_seal)) *
                                      RTMB::fft((Q)), inverse = TRUE)))/no_w_full
  #pred_rate[pred_rate < 1e-18] <- 0
  return(pred_rate)
}

getSealMortRTMB <- function(pred_rate , interaction_seal = NULL , w = NULL) {
  f()
  pred_mort <- interaction_seal %*% pred_rate[, 1:length(w), drop = FALSE]
  return(pred_mort)
}

## Total Mortality ----

getMortRTMB <- function(pred_mort , f_mort , seal_mort , z = NULL , mu_b = NULL) {
  f()
  mort <- pred_mort + (mu_b * z) + f_mort + seal_mort
  return(mort)
}

## RDI & RDD ----

getRDIRTMB <- function(e_repro , n , w = NULL , dw = NULL, erepro = NULL, w_min_idx = NULL) {
  f()
  e_repro_pop <- drop((e_repro * n) %*% dw)
  rdi <- 0.5 * (e_repro_pop * erepro)/w[w_min_idx]
  return(rdi)
}

getRDDRTMB <- function(rdd) {
  return(rdd)
}

## Resource Mort ----

getSealResourceMortRateRTMB <- function(pred_rate , w = NULL, resource_interaction_seal = NULL) {
  f()
  pred_mort <- resource_interaction_seal %*% pred_rate[,1:length(w), drop = FALSE]
}

getResourceMortRateRTMB <- function(pred_rate , interaction_resource = NULL) {
  f()
  as.vector(interaction_resource %*% pred_rate)
}

getResourceMortRTMB <- function(resource_mort_rate , resource_mort_rate_seals , w = NULL) {
  f()
  idx <- tail(1:length(resource_mort_rate), length(w))
  resource_mort_rate[idx] <- resource_mort_rate[idx] + c(resource_mort_rate_seals)
  return(resource_mort_rate)
}

## Seal Reproduction ----

getSealDietRTMB <- function(n , n_pp , n_other , feeding_level ,
                            interaction_seal = NULL , resource_interaction_seal = NULL , w_full = NULL ,
                            w = NULL , dw = NULL , w_seal = NULL , dw_seal = NULL, search_vol_seal  = NULL,
                            ft_pred_kernel_e_real_seal = NULL, ft_pred_kernel_e_imag_seal = NULL , idx_seal = NULL ,
                            type = 'diet_by_sp') {
  f()
  idx_sp <- 1:length(w)
  seal_interaction <- c(interaction_seal, resource_interaction_seal)
  no_sp <- nrow(n) + 1
  no_w_full <- length(w_seal)
  n_all <- matrix(NA, nrow = nrow(n) + 1, ncol = ncol(n))
  n_all[1:nrow(n), ] <- n
  n_all[nrow(n_all), ] <- n_pp[w_full %in% w]
  prey <- matrix(0, nrow = no_sp, ncol = no_w_full)
  prey[1:no_sp, idx_sp] <- sweep(n_all, 2, w * dw,
                                 "*")
  ft_pred_kernel_e <- as.complex(ft_pred_kernel_e_real_seal) + (0+1i) * as.complex(ft_pred_kernel_e_imag_seal)
  ft <- array(0, dim = c(1, no_w_full, no_sp))
  diet <- matrix(0 , no_sp , no_w_full)
  for(i in 1:no_sp) {
    ae <- Re(RTMB::fft((as.complex(ft_pred_kernel_e_real_seal) + (0+1i) * as.complex(ft_pred_kernel_e_imag_seal)) *
                         RTMB::fft(t(prey)[,i]) , inverse = T)/no_w_full)
    diet[i,] <- ae
  }

  for (i in 1:nrow(diet)) {
    diet[i, ] = diet[i, ] * search_vol_seal * seal_interaction[i]
  }
  diet <- diet * matrix(data = (1 - feeding_level), nrow = no_sp,
                        ncol = no_w_full, byrow = T)
  diet[,-idx_seal] <- 0
  out <- list(
    diet = diet,
    total_consumption = sum(sweep(diet, 2, n_other * dw_seal, "*")),
    diet_by_sp = rowSums(sweep(diet, 2, n_other * dw_seal, "*")),
    percent_by_sp = rowSums(sweep(diet, 2, n_other * dw_seal, "*")) / sum(sweep(diet, 2, n_other *
                                                                                  dw_seal, "*"))
  )
  return(out[[type]])
}


getSealReproRTMB <- function(n_other , prey , alpha_seal , w_seal = NULL , dw_seal = NULL ,
                             Msrr_seal = NULL , Jmax_seal = NULL , e_seal = NULL , resource_seal = NULL ,
                             M_seal = NULL , H = NULL , idx_seal = NULL) {
  f()
  N <- n_other
  Bm0 <- N*w_seal*dw_seal
  tot <- sum(c(prey,resource_seal), na.rm = T)
  num <- tot^2 ; den <- Jmax_seal+tot^2

  ret <- (Bm0 + Bm0*(-Msrr_seal+e_seal*Jmax_seal*(num/den)) - M_seal*(Bm0)^2 - H) + alpha_seal
  ret <- ret/(w_seal*dw_seal)
  ret[-idx_seal] <- 0
  return(ret)
}

getResourceDynamicsRTMB <- function(n_pp, resource_mort , dt = NULL , rr_pp = NULL , cc_pp = NULL) {
  f()
  mur <- rr_pp + resource_mort
  n_steady <- rr_pp * cc_pp/mur
  n_pp_new <- n_steady + (n_pp - n_steady) * exp(-mur * dt)
  #sel <- !is.finite(n_pp_new)
  #n_pp_new[sel] <- n_pp[sel]
  n_pp_new
}


# ratesToAD <- function(rates) {
#   rates <- deparse(rates_fns$Rates)
#   rates <- paste(rates , collapse = '\n')
#   rates <- gsub(', \\n' , ', ' , rates)
#   rates <- gsub(' {2,}' , ' ' , rates)
#   rates <- gsub('rates_fns\\$(?!E|Fe|M)' , 'RTMB::AD(rates_fns\\$' , rates , perl = T)
#   rates <- strsplit(rates , '\\n')[[1]]
#   rates[sapply(rates , grepl , pattern = 'AD')] <- gsub('\\.{3}\\)$' , '\\.\\.\\.\\)\\)' , rates[sapply(rates , grepl , pattern = 'AD')])
#   text <- paste0('Rates <-' , paste(rates , collapse = '\n'))
#   eval(parse(text = text))
#   return(Rates)
# }
