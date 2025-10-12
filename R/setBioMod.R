## BioEnergetic Model for Seals

setBioMod <- function(Bm0 , Msrr , Jmax , e , a , resource , M , H , prey , ...) {
  tot <- sum(c(prey,resource), na.rm = T)
  num <- a*tot^2 ; den <- Jmax+a*tot^2
  ret <- Bm0 + Bm0*(-Msrr+e*Jmax*(num/den)) - M(Bm0)^2 - H
  return(ret)
}

# Parameters to add to setSeals

# H
# Msrr
# Jmax
# e
# a
# resource
#
