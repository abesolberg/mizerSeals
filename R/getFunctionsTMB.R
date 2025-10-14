library(TMB)
compile ("src/predation_dev.cpp",libinit = FALSE ) ## notice flag
dyn.load(dynlib ("src/predation_dev"))


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

tmb_getPredRate <- function(params , feeding_level , n = params@initial_n , search_vol = params@search_vol , dw = params@dw , ft_pred_kernel_p_real = NULL ,
                            ft_pred_kernel_p_imag = NULL , ft_mask = params@ft_mask) {
  .Call("call_getPredRate", n , w_full , w , feeding_level , search_vol , dw , ft_pred_kernel_p_real , ft_pred_kernel_p_imag , ft_mask)
}


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

fl <- tmb_getFeedingLevel(NS_params@intake_max , getEncounter(params))
pr <- tmb_computePrey(params)
Q <- tmb_computeQMatrix(params , feeding_level = fl)

fl_seals <- tmb_getFeedingLevel(matrix(params@other_params$sealParams$intake_max , nrow = 1), 
                    matrix(getSealEncounter(params) , nrow = 1))

all(fl == getFeedingLevel(params))
all(fl_seals == getSealFeedingLevel(params , getSealEncounter(params)))

n <- as.matrix(cbind(params@initial_n , matrix(0 , nrow = nrow(params@initial_n) , ncol = 7)))
n_pp <- c(tail(params@initial_n_pp , 100) , rep(0 , 7)) 

x <- getSealPredRate(params)

tmb_computePrey()

tmb_computePrey(
  params ,
  n =  as.matrix(cbind(params@initial_n , matrix(0 , nrow = nrow(params@initial_n) , ncol = 7))) , 
  n_pp = c(tail(params@initial_n_pp , 100) , rep(0 , 7)) , 
  interaction_resource = params@other_params$sealParams$resource_interaction_seal ,
  interaction = matrix(params@other_params$sealParams$interaction_seal , nrow = 1) ,
  w_full = params@other_params$sealParams$w ,
  dw_full = params@other_params$sealParams$dw ,
  w = params@other_params$sealParams$w
)

sp <- params@other_params$sealParams
list2env(sp , globalenv())
n <- params@initial_n ; n_pp <- params@initial_n_pp


computePrey(resource_interaction_seal,
            n_pp,
            n,
            interaction_seal,
            w_full,
            w,
            dw_full) {
  ret <- matrix(NA , nrow = length(resource_interaction_seal) , ncol = length(n_pp))
  matrix<Type> ret(species_params_interaction_resource.size(), n_pp.size());
  for(i in 1:nrow(ret)){
    for(j in 1:ncol(ret)){
      ret[i,j] = resource_interaction_seal[i] * n_pp[j]
    }
  }
  
  ret.block(0, w_full.size() - w.size(), ret.rows(), w.size()) += species_params_interaction * n;
  vector<Type> tmp = w_full * dw_full;
  for (int i = 0; i < ret.rows(); ++i) {
    ret.row(i) = vector<Type>(ret.row(i)) * tmp;
  }
  return ret;
}
