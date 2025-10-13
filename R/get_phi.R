# Get values needed for pred kernel -- Not Exported (Rename this to not conflict with mizer)
get_phi <- function (ppmr, beta, sigma) {
  Beta <- log(beta)
  phi <- exp(-(log(ppmr) - Beta)^2/(2 * sigma^2))
  return(phi)
}
