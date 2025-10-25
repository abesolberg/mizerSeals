library(mizer)
library(mizerSeals)
library(ggplot2)

# install.packages('mizerSeals_0.1.0.tar.gz')

# https://publications.gc.ca/collections/collection_2020/mpo-dfo/fs70-5/Fs70-5-2020-014-eng.pdf
# https://www.frontiersin.org/journals/marine-science/articles/10.3389/fmars.2021.579946/full

# sapply(list.files('R/', , full.names=T) , source)

## To Do:

# Create function to make all params
  # Function to sort effort and annual rates
  # Function should basically set everything as a paramater and then map out whatever isn't needed
    # Map when catch = 0, sex = M, ...
  # Need to create helper function to put everything back into the param object. 
  # Params:
    # Rmax, MuB multiplier, Recruitment anomaly
    # Interaction matrix, resource interaction (so you can potentially fit to diet)
    # Catchability, Effort
    # Seal Diet interaction & resource interaction
    # Seal bioenergetic rates

getNewVals <- function(params , newvals = list()) {
  for(i in names(newvals)) {
    val <- slot(params , i)
    if (is.list(newvals[[i]])) {
      .newvals <- newvals[[i]]
      for (k in names(.newvals)) {
        val[[k]][] <- .newvals[[k]]
      }
    } else val[] <- newvals[[i]]
    slot(params , i) <- val
  }
  return(params)
  #sp_vals <- names(params@species_params)
  #other_vals <- names(params@other_params)
  #seal_params <- names(params@seal_params)
  #trawl_params <- names(trawl_params)
}

params <- mizer::NS_params
params@species_params$q

newvals <- list(
  species_params = list(p = 1 , q = runif(12 , 0 , 1)) ,
  catchability = runif(12*4 , 0 , 1)
)

test <- getNewVals(params , newvals)
test@species_params$p
test@species_params$q
test@catchability


sealParams <- setSealParams(
  params , w_max_seal = 150000 , w_min_seal = 11000 ,
  interaction_seal = rep(1 , length(params@species_params$species)) ,
  beta = 5000 , sigma = 3 , resource_interaction_seal = 100 ,
  time_step = 5 , dynamicSeals = F , h = 250
)
sealParams$initialSealN[,107] <- rnorm(5 , 20000000 , 25000)/sealParams$dw[107]
#sealParams$initialSealN[,107] <- 1:100
params <- setSeals(mizer::NS_params , sealParams)
sim <- project_simple(params , dt = .1 , steps = 1)

plot(sim)
setBioMod

steps = 100 ; dt = .1

n_pp <- params@initial_n_pp
n_other <- params@initial_n_other
n <- params@initial_n
effort <- params@initial_effort

fun <- function(params, n = params@initial_n, n_pp = params@initial_n_pp, 
                n_other = params@initial_n_other, effort = params@initial_effort, 
                t = 1, dt = 0.1, steps = 100 , resource_dynamics_fn = get(params@resource_dynamics), 
                other_dynamics_fns = lapply(params@other_dynamics, get), 
                rates_fns = lapply(params@rates_funcs, get), ...) { #rmax , mu_b , effort , seal_re , ...) { # Need to add params, t, dt, steps ...
  #params@species_params$R_max <- R_max # New Parameter
  #params@mu_b <- params@mu_b * mu_b # New Scaling Factor
  
  #params <- getNewVals(params , newvals = list())
  
  no_sp <- nrow(params@species_params)
  no_w <- length(params@w)
  idx <- 2:no_w
  years <- steps*dt # steps should be calculated by length of effort
  
  f_at_dt <- n_at_dt <- array(0 , dim = c(no_sp , no_w , steps))
  seals_dt <- n_other$seals
  
  w_min_idx_array_ref <- (params@w_min_idx - 1) * no_sp + (1:no_sp)
  a <- matrix(0, nrow = no_sp, ncol = no_w)
  b <- matrix(0, nrow = no_sp, ncol = no_w)
  S <- matrix(0, nrow = no_sp, ncol = no_w)
  for (i_time in 1:steps) {
    # Set Effort, Seals RE, etc... 
    r <- rates_fns$Rates(params, n = n, n_pp = n_pp, n_other = n_other, 
                         t = t, effort = effort, rates_fns = rates_fns , ...)
    n_other_new <- list()
    for (component in names(params@other_dynamics)) {
      n_other_new[[component]] <- other_dynamics_fns[[component]](params, 
                                                                  n = n, n_pp = n_pp, n_other = n_other, rates = r, 
                                                                  t = t, dt = dt, component = component , ...)
    }
    
    n_pp <- resource_dynamics_fn(params, n = n, n_pp = n_pp, 
                                 n_other = n_other, rates = r, t = t, dt = dt, resource_rate = params@rr_pp, 
                                 resource_capacity = params@cc_pp , ...)
    a[, idx] <- sweep(-r$e_growth[, idx - 1, drop = FALSE] * 
                        dt, 2, params@dw[idx], "/")
    b[] <- 1 + sweep(r$e_growth * dt, 2, params@dw, "/") + 
      r$mort * dt
    S[, idx] <- n[, idx, drop = FALSE]
    n[w_min_idx_array_ref] <- (n[w_min_idx_array_ref] + r$rdd * 
                                 dt/params@dw[params@w_min_idx])/b[w_min_idx_array_ref]
    n <- mizer:::inner_project_loop(no_sp = no_sp, no_w = no_w, n = n, 
                            A = a, B = b, S = S, w_min_idx = params@w_min_idx)
    n_other <- n_other_new
    t <- t + dt
    
    n_at_dt[,,i_time] <- n # Can figure out how to do this otherwise
    f_at_dt[,,i_time] <- n*sim$rates$f_mort*dt # Can figure out how to do this otherwise
    seals_dt <- c(seals_dt , n_other$seals)
  }
  n_at_y <- n_at_dt[,,seq(1 , i_time , by = years)]
  y <- rep(seq_along(1:years), each = years*dt)
  f_at_y <- sapply(1:years , function(x) (apply(f_at_dt[,,which(y == x)] , 1 , sum))) # This puts out the wrong format, also not aggregated
  # Trawl at T
  
  # Seal Consumption Reports (maybe needed, but can just do it after...)
  
  ## Fitting Functions
  # Fit to catch
  
  # Fit to trawl
  
  # Fit RE F (normal or lognormal, allow for switch)
  
  # Fit RE Recruitment (Normal or lognormal distribution)
  
  # Fit RE Mu_b (Normal distribution, random walk)
  
  # Fit to seal biomass
  
  # Fit Seal RE (Normal or Lognormal distribution)
  
  
  return(list(
    n = n_at_y, 
    f = f_at_y ,
    n_other = seals_dt 
  ))
   
}

dim(f_at_dt)
n_at_y

x <- fun(params)

x$n_other == params@initial_n_other$seals[,107]

dim(x$n)
dim(x$f)

a <- array(1 , dim = c(5 , 10 , 100))

x <- c(3, 5, 2)
y <- c(3, 2, 1, 1, 2, 3, 4, 5, 4, 5)
z <- c(2, 4, 8, 1, 5)

g <- rep(seq_along(1:10), each = 10/1)
sapply(1:10 , function(x)apply(a[,,which(g == x)] , 1 , sum))

tapply(a[,,], g, rowSums)



x = 1
plot(sim)

diet <- getSealDiet(params)
diet$diet_by_sp
diet$percent_by_sp ; sum(diet$percent_by_sp)
mort <- getSealMort(params)
resource_mort <- getSealResourceMort(params)

plot(resource_mort)
plot(mort[2,])

plotSealMort(params , return_data = F)

max(mort)

sealParams2 <- setSealParams(
  mizer::NS_params , w_max_seal = 150000 , w_min_seal = 11000 ,
  interaction_seal = rep(1 , length(params@species_params$species)) ,
  beta = 1000 , sigma = 3 , resource_interaction_seal = 0 ,
  time_step = 1 , dynamicSeals = F
)
sealParams2$initialSealN[,107] <- rnorm(1 , 20000000 , 0)/sealParams$dw[107]
#sealParams$initialSealN[,107] <- 1:100
params2 <- setSeals(mizer::NS_params , sealParams2)
sim2 <- project(params2 , t_max = 5)

diet2 <- getSealDiet(params2)
mort2 <- getSealMort(params2)
resource_mort2 <- getSealResourceMort(params2)

all.equal(resource_mort2 , resource_mort)

ggplot() +
  geom_line(aes(x = params@w_full , y = resource_mort) , color = 'black') +
  geom_line(aes(x = params@w_full , y = resource_mort2) , color = 'red') +
  scale_x_log10()

getSealDiet(params)$total_consumption
getSealDiet(params2)$total_consumption

diet2$diet[13,107]

plot(resource_mort2 , resource_mort)



## Test
sapply(list.files('R/',full.names=T) , source)

sealParams <- setSealParams(
  params , w_max_seal = 150000 , w_min_seal = 11000 ,
  interaction_seal = rep(1 , length(params@species_params$species)) ,
  beta = 1000 , sigma = 3 , resource_interaction_seal = 0 ,
  time_step = 1 , dynamicSeals = F
)
sealParams$initialSealN[,107] <- rnorm(1 , 20000000 , 0)/sealParams$dw[107]
#sealParams$initialSealN[,107] <- 1:100
params <- setSeals(params , sealParams)
sim <- project(params , t_max = 5)

diet_res0 <- getSealDiet(params)
mort_res0 <- getSealMort(params)

diet_res0$total_consumption == diet$total_consumption
all(diet$diet == diet_res0$diet[1:12,])

plot(sim , return_data = F , only_mature = F)

sapply(1:6 , function(x)sim@n_other[x,][[1]][,107])

plotSealPredKernel(params)
getSealDiet(params)$total_consumption/1e12
unlist(getSealDietSim(sim , 'total_consumption'))/1e12
unlist(getSealDietSim(sim2 , 'total_consumption'))/1e12



plot(sim)

outMort <- getSealMortSim(sim)
outDiet <- getSealDietSim(sim)
list2array(outDiet)

getPredMort



