## MizeRTMB Predation 

setClass(
  'MizerSealParams' ,
  slots = c(
    getSlots('MizerParams') ,
    sealParams = "list" ,
    trawlParams = "list" ,
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

