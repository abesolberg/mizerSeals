## BioEnergetic Model for Seals

setBioMod <- function(N , Msrr , Jmax , e , resource , M , H , prey , .alpha , w , dw , idx ,...) {
  Bm0 <- N*w*dw
  tot <- sum(c(prey,resource), na.rm = T)
  num <- tot^2 ; den <- Jmax+tot^2
  ret <- (Bm0 + Bm0*(-Msrr+e*Jmax*(num/den)) - M*(Bm0)^2 - H[1,]) + .alpha[1,]
  ret <- ret/(w*dw)
  ret[-idx] <- 0
  return(ret)
}

# Parameters to add to setSeals

# https://www.int-res.com/articles/meps_oa/m379p279.pdf
# https://www.int-res.com/articles/meps_oa/m511p265.pdf

# Basal metabolic rate(293*BMi ^-0.75; Kleiber 1975
# https://www-jstor-org.qe2a-proxy.mun.ca/stable/pdf/2462335.pdf?refreqid=fastly-default%3A6046f102ad5851ec920f0b6897a15c52&ab_segments=&initiator=&acceptTC=1

#At = 54.9 for endotherms
#Aj = 89.4 for Endotherms

# Mssr = 54.9*(params@other_params$seal_params$w[idx])^-0.25
# Jmax = 89.4*(params@other_params$seal_params$w[idx])^-0.25
# e = 0.89
# M = .2
# a = Not needed in this model
  
# Seal Params:
# H = Harvest (given)
# Msrr = Mass Specific Respiration Rate Using Average Sized Seal (given) ... see eq 5 & 6 https://www.int-res.com/articles/meps_oa/m511p265.pdf
# Jmax = Maximum mass-specific ingestion rate for the stock (given)
# e = 0.83 Fraction of ingested energy/biomass available at the metabolizable level (this could be done for multiple species) (https://publications.gc.ca/collections/collection_2013/mpo-dfo/Fs70-5-2012-156-eng.pdf)
# a = functional response coefficient (“attack rate”) - Not needed because predation is being calculated
# M = mortality coefficent
# Resource = 0 (given, could be estimated but we now have resource predation)

# Msrr = 54.9*(params@other_params$seal_params$w[idx])^-0.25
# Jmax = 89.4*(params@other_params$seal_params$w[idx])^-0.25
# e = 0.89
# M = .000001
# 
# setBioMod <- function(Bm0 , Msrr , Jmax , e , resource , M , H , prey , ...) {
#   tot <- sum(c(prey,resource), na.rm = T)
#   num <- tot^2 ; den <- Jmax+tot^2
#   ret <- Bm0 + Bm0*(-Msrr+e*Jmax*(num/den)) - M*(Bm0)^2 - H
#   return(ret)
# }
# 
# out <- c(seals$BiomassKT[1])
# for (i in 2:nrow(seals)) {
#   out <- c(
#     out ,
#     setBioMod(
#       seals$BiomassKT[i - 1] * 1000 ,
#       Msrr = Msrr ,
#       Jmax = Jmax ,
#       a = 1 ,
#       resource = 0,
#       e = e ,
#       M = M ,
#       H = seals$CatchBMKT[i - 1] * 1000 ,
#       prey = seals$FoodEstimateKT[i - 1] * 1000
#     )
#   )
# }
# 
# sum((log(seals$BiomassKT[-1]) - log(out[-1]/1000))^2)