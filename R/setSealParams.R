# Set Seal Life History
setSealParams <- function(
    params ,
    w_max_seal ,
    w_min_seal ,
    interaction_seal,
    resource_interaction_seal,
    a = .35, q = .8, n = .75,
    beta = 50, sigma = 3, f0 = 0.6, h = 30 ,
    gamma = NULL , 
    time_steps = 1 ,
    initialSealN = NULL ,
    dynamicSeals = F ,
    ...
) {

  dx <- log10(params@dw/params@w +1)[1]
  w <- c(params@w[-length(params@w)] , 10^(seq(from = log10(max(params@w)), to = log10(w_max_seal) , by = dx)))
  dw <- (10^dx - 1) * w

  if (is.null(initialSealN)) initialSealN <- array(0 , dim = c(time_steps , length(w)))
  if (!all.equal(dim(initialSealN), c(time_steps , length(w)))) stop('Seal Initial N does not have proper dimensions.')
  
  pred <- sealPredKernel(params , w , beta , sigma)
  
  if (is.null(gamma)) {
    gamma <- getSealGamma(
      params ,
      w ,
      q ,
      h ,
      f0 ,
      interaction_seal = interaction_seal ,
      resource_interaction_seal = resource_interaction_seal ,
      ft_pred_kernel_e = pred$ft_pred_kernel_e ,
      dw = dw
    )
  }
  
  search_vol <- getSealSearchVol(
    params ,
    w = w ,
    q = q ,
    gamma = gamma 
  )
 
  return(
    list(
      ft_pred_kernel_e = pred$ft_pred_kernel_e ,
      ft_pred_kernel_p = pred$ft_pred_kernel_p ,
      pred_kernel = pred$pred_kernel ,
      intake_max = h*w^n ,
      search_vol = search_vol , 
      w = w ,
      dw = dw ,
      initialSealN = array(0 , dim = c(time_steps , length(w))) ,
      interaction_seal = interaction_seal,
      resource_interaction_seal = resource_interaction_seal ,
      a = a,
      q = q ,
      n= n ,
      beta = beta ,
      sigma = sigma,
      f0 = f0 ,
      h = h,
      gamma = gamma , 
      dynamicSeals = dynamicSeals ,
      ...
    )
  )
}

