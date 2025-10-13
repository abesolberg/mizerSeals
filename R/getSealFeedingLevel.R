# Outputs seal feeding level
# https://sizespectrum.org/mizer/reference/setMaxIntakeRate.html?q=intake#null

getSealFeedingLevel <- function(params , encounter) {
  sp <- params@other_params$sealParams
  intake_max <- sp$h*sp$w^sp$n
  return(encounter/(encounter + intake_max))
}
