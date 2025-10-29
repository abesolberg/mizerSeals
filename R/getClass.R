## MizeRTMB Predation

setClass(
  'MizerSealParams' ,
  slots = c(
    getSlots('MizerParams') ,
    seal_params = "list" ,
    trawl_params = "list" ,
    ft_pred_kernel_e_real = 'array',
    ft_pred_kernel_e_imag = 'array',
    ft_pred_kernel_p_real = 'array' ,
    ft_pred_kernel_p_imag = 'array'
  ) , contains = 'MizerParams'
)

S4toList <- function(obj) {
  sn <- slotNames(obj)
  structure(lapply(sn, slot, object = obj), names = sn)
}

setGeneric("validParams")
setMethod("validParams" , "MizerSealParams" , definition = function(params) return(params))
