# Outputs seal feeding level
# https://sizespectrum.org/mizer/reference/setMaxIntakeRate.html?q=intake#null
getSealFeedingLevel <- function(w_seal , seal_encounter , h = 30 , n = .7) {
  intake_max <- h*w_seal^n
  return(seal_encounter/(seal_encounter + intake_max))
}