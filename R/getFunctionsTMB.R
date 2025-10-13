library(TMB)
compile ("src/predation_dev.cpp",libinit = FALSE ) ## notice flag
dyn.load(dynlib ("src/predation_dev"))


# Feeding Level
tmb_getFeedingLevel <- function(intake_max , encounter){
  .Call("call_getFeedingLevel",intake_max , encounter)
}

# Compute Prey
tmb_computePrey <- function(params , n = params@initial_n , n_pp = params@initial_n_pp , interaction_resource = params@species_params$interaction_resource, 
                            interaction = params@interaction , w_full = params@w_full , w = params@w , dw_full = params@dw_full){
  .Call("call_computePrey" ,interaction_resource, n_pp, n, interaction, w_full, w, dw_full)
}

# Q Matrix
tmb_computeQMatrix <- function(params , feeding_level , search_vol = params@search_vol , n = params@initial_n , w = params@w , w_full = params@w_full , dw = params@dw ) {
  .Call("call_computeQMatrix" , n, w_full , w , feeding_level , search_vol , dw)
}

tmb_getPredRate <- function(params , feeding_level , n = params@initial_n , search_vol = params@search_vol , dw = params@dw , ft_pred_kernel_p_real = NULL ,
                            ft_pred_kernel_p_imag = NULL , ft_mask = params@ft_mask) {
  .Call("call_getPredRate", n , w_full , w , feeding_level , search_vol , dw , ft_pred_kernel_p_real , ft_pred_kernel_p_imag , ft_mask)
}


library(mizer)
params <- NS_params

fl <- tmb_getFeedingLevel(NS_params@intake_max , getEncounter(params))
pr <- tmb_computePrey(params)
Q <- tmb_computeQMatrix(params , feeding_level = fl)

Q

all(fl == getFeedingLevel(params))
pr
