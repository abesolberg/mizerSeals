# Set Seals in Mizer -- Exported
setSeals <- function(params, sealParams) {
  if(!validObject(params)) params <- validParams(params)
  params@other_params$sealParams <- sealParams
  # Hook into mizer
  params@initial_n_other[["seals"]] <- sealParams$initialSealN[1,,drop=F]
  params@other_dynamics[["seals"]] <- "getSealRepro"
  params@other_mort[["seals"]] <- "getSealMort"
  return(params)
}
