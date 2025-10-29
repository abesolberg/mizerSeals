# Set Seals in Mizer -- Exported
setSeals <- function(params, seal_params) {
  if(!validObject(params)) params <- validParams(params)
  # Hook into mizer
  params <- `slot<-`(params , 'seal_params' , check = F , NULL)
  params@initial_n_other[["seals"]] <- seal_params$initialSealN[1,,drop=F]
  params@other_dynamics[["seals"]] <- "getSealRepro"
  params@other_mort[["seals"]] <- "getSealMort"
  params <- setRateFunction(params, "ResourceMort", "getSealResourceMort")
  params <- setRateFunction(params, "PredRate" , 'getTMBPredRate')
  params <- setRateFunction(params , "Encounter" , 'getTMBEncounterRate')
  params <- do.call(new , append(list(Class = 'MizerSealParams') , S4toList(params)))
  params@seal_params <- seal_params # Note this change throughout
  params@ft_pred_kernel_e_real <- Re(params@ft_pred_kernel_e)
  params@ft_pred_kernel_e_imag <- Im(params@ft_pred_kernel_e)
  params@ft_pred_kernel_p_real <- Re(params@ft_pred_kernel_p)
  params@ft_pred_kernel_p_imag <- Im(params@ft_pred_kernel_p)
  return(params)
}
