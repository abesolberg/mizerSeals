library(TMB)
library(mizer)

sapply(list.files('R/', , full.names=T)[!grepl('TMB' , list.files('R/'))] , source)

params <- mizer::NS_params

sealParams <- setSealParams(
  mizer::NS_params , w_max_seal = 150000 , w_min_seal = 11000 ,
  interaction_seal = rep(1 , length(params@species_params$species)) ,
  beta = 5000 , sigma = 3 , resource_interaction_seal = 100 ,
  time_step = 1 , dynamicSeals = F , h = 250
)

sealParams$initialSealN[,107] <- rnorm(1 , 20000000 , 0)/sealParams$dw[107]
params <- setSeals(mizer::NS_params , sealParams)

ft_pred_kernel_real <- Re(params@ft_pred_kernel_e)
ft_pred_kernel_imag <- Im(params@ft_pred_kernel_e)
ft_pred_kernel_p_imag <- Im(params@ft_pred_kernel_p)
ft_pred_kernel_p_real <- Re(params@ft_pred_kernel_p)


dat <- list(
  # Initial Values
  n = params@initial_n ,
  n_pp = params@initial_n_pp ,
  
  # Bins
  w_full = params@w_full ,
  dw_full = params@dw_full , 
  w = params@w ,
  dw = params@dw ,
  ft_mask = params@ft_mask ,
  
  # Pred Kernel
  ft_pred_kernel_p_real = ft_pred_kernel_p_real ,
  ft_pred_kernel_p_imag = ft_pred_kernel_p_imag,
  ft_pred_kernel_real = ft_pred_kernel_real ,
  ft_pred_kernel_imag = ft_pred_kernel_imag,
  
  # Predation Data Values
  search_vol = params@search_vol ,
  intake_max = params@intake_max ,
  
  # Interactions
  species_interaction_resource = params@species_params$interaction_resource ,
  species_interaction = params@interaction ,
  
  # FFT Dims
  nRows = as.integer(length(params@w_full)) , 
  nCols = as.integer(length(params@w_full))
)

compile('src/model_test_full.cpp')

compile ("dev/predation_dev.cpp",libinit = FALSE , PKG_CXXFLAGS = "-DEIGEN_FFTW_DEFAULT", PKG_LIBS = "-lfftw3 -lfftw3f") ## notice flag
dyn.load(dynlib ("dev/predation_dev"))

obj <- MakeADFun(dat, parameters = list(A = 1) ,  DLL="predation_dev")
dyn.unload('predation_dev')



# Feeding Level -- works for both species and seals
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

# # Pred Rate
# matrix<Type> n,
# vector<Type> w_full,
# vector<Type> w,
# matrix<Type> feeding_level,
# matrix<Type> search_vol,
# vector<Type> dw,
# matrix<Type> ft_pred_kernel_p_real,
# matrix<Type> ft_pred_kernel_p_imag,
# matrix<Type> ft_mask

tmb_getPredRate <- function(params , feeding_level , n = params@initial_n , search_vol = params@search_vol , w = params@w , dw = params@dw , ft_pred_kernel_p_real = NULL ,
                            ft_pred_kernel_p_imag = NULL , ft_mask = params@ft_mask) {
  .Call("call_getPredRate", n , w_full , w , feeding_level , search_vol , dw , ft_pred_kernel_p_real , ft_pred_kernel_p_imag , ft_mask)
}

# Encounter
tmb_getEncounter <- function(params , n = params@initial_n , n_pp = params@initial_n_pp , w_full = params@w_full , w = params@w, 
                             dw_full = params@dw_full , ft_pred_kernel_real , ft_pred_kernel_imag , search_vol = params@search_vol ,
                             species_params_interaction_resource = params@species_params$interaction_resource , 
                             species_params_interaction = params@interaction) {
  .Call("call_getEncounter" , species_params_interaction_resource , n_pp , n , species_params_interaction ,w_full , w ,
                    dw_full ,ft_pred_kernel_real , ft_pred_kernel_imag , search_vol)
}


## Helper Funs ----
tmb_convertNppUp <- function(npp , w , w_new) {
  .Call("call_convertDimsUpNpp" , npp , w , w_new)
}
tmb_convertNppDown <- function(npp , w , npp_convert) {
  .Call("call_convertDimsDownNpp" , npp , w , npp_convert)
}

tmb_convertNUp <- function(n , w_new) {
  .Call("call_convertDimsUpN" , n , w_new)
}
tmb_convertNDown <- function(n , n_convert) {
  .Call("call_convertDimsDownN" , n , n_convert)
}

## Compare Outputs ----

library(mizer)

sapply(list.files('R/', , full.names=T)[!grepl('TMB' , list.files('R/'))] , source)

params <- mizer::NS_params

sealParams <- setSealParams(
  mizer::NS_params , w_max_seal = 150000 , w_min_seal = 11000 ,
  interaction_seal = rep(1 , length(params@species_params$species)) ,
  beta = 5000 , sigma = 3 , resource_interaction_seal = 100 ,
  time_step = 1 , dynamicSeals = F , h = 250
)

sealParams$initialSealN[,107] <- rnorm(1 , 20000000 , 0)/sealParams$dw[107]
params <- setSeals(mizer::NS_params , sealParams)

make_complex <- function(x, inverse=F){
  if(inverse==F){
    tmp <- ncol(x)/2
    return(matrix(complex(real=x[,1:tmp],imaginary = x[,-(1:tmp)]),nrow=nrow(x),ncol=tmp))
  }
  return(cbind(Re(x),Im(x)))
}

ft_pred_kernel_real <- Re(params@ft_pred_kernel_e)
ft_pred_kernel_imag <- Im(params@ft_pred_kernel_e)
ft_pred_kernel_p_imag <- Im(params@ft_pred_kernel_p)
ft_pred_kernel_p_real <- Re(params@ft_pred_kernel_p)
ft_pred_kernel_e <- make_complex(params@ft_pred_kernel_e,TRUE)



tmb_getEncounter(params , ft_pred_kernel_real = ft_pred_kernel_real , ft_pred_kernel_imag = ft_pred_kernel_imag)

tmb_getFeedingLevel(NS_params@intake_max , tmb_getEncounter(params , ft_pred_kernel_real = ft_pred_kernel_real , ft_pred_kernel_imag = ft_pred_kernel_imag))

tmb_computePrey(params)

tmb_computeQMatrix(params , feeding_level = tmb_getFeedingLevel(NS_params@intake_max , getEncounter(params)))

encounter = tmb_getEncounter(params , ft_pred_kernel_real = ft_pred_kernel_real , ft_pred_kernel_imag = ft_pred_kernel_imag)

plot(encounter[1,])
plot(getEncounter(params)[1,])
dim(getEncounter(params))

mat <- as.matrix(apply(params@ft_mask , 1:2 , as.numeric))

tmb_getPredRate <- function(params , feeding_level , n = params@initial_n , search_vol = params@search_vol , w = params@w , dw = params@dw , ft_pred_kernel_p_real  ,
                            ft_pred_kernel_p_imag , ft_mask = params@ft_mask) {
  .Call("call_getPredRate", n , w_full , w , feeding_level , search_vol , dw , ft_pred_kernel_p_real , ft_pred_kernel_p_imag , ft_mask)
}

tmb_getPredRate(params = params , 
                feeding_level = tmb_getFeedingLevel(params@intake_max , tmb_getEncounter(params , ft_pred_kernel_p_real = ft_pred_kernel_p_real , ft_pred_kernel_p_imag = ft_pred_kernel_p_imag))
                                                                      ))

apply(params@ft_mask , 1:2 , as.numeric)

tmb_getPredRate <- function(params , feeding_level , n = params@initial_n , search_vol = params@search_vol , dw = params@dw , ft_pred_kernel_p_real = NULL ,
                            ft_pred_kernel_p_imag = NULL , ft_mask = params@ft_mask) {
  .Call("call_getPredRate", n , w_full , w , feeding_level , search_vol , dw , ft_pred_kernel_p_real , ft_pred_kernel_p_imag , apply(params@ft_mask , 1:2 , as.numeric))
}

tmb_getPredRate(params , feeding_level = tmb_getFeedingLevel(params@intake_max , getEncounter(params)) ,
                ft_pred_kernel_p_real = ft_pred_kernel_p_real , ft_pred_kernel_p_imag = ft_pred_kernel_p_imag)

## Testing Functions for Seals ----
fl_seals <- tmb_getFeedingLevel(matrix(params@other_params$sealParams$intake_max , nrow = 1), 
                    matrix(getSealEncounter(params) , nrow = 1))

# This works for seals
tmb_computePrey(
  params ,
  n =  tmb_convertNUp(params@initial_n , params@other_params$sealParams$w) , 
  n_pp = tmb_convertNppUp(params@initial_n_pp , params@w , params@other_params$sealParams$w) , 
  interaction_resource = params@other_params$sealParams$resource_interaction_seal ,
  interaction = matrix(params@other_params$sealParams$interaction_seal , nrow = 1) ,
  w_full = params@other_params$sealParams$w ,
  dw_full = params@other_params$sealParams$dw ,
  w = params@other_params$sealParams$w
)

prey <- tmb_computePrey(params)

plot(prey[6,])

