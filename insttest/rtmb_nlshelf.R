## NL Shelf 

## Load Libraries

library(mizer)
library(mizerSeals)
library(RTMB)
compiler::enableJIT(0);

## Hack Valid Params ----

# assignInNamespace('validParams' , mizerSeals::validParams , ns = 'mizer')

## Read Data ----

path <- 'data/inputs/'

init_params <- readRDS(paste0(path , 'initial_steady_params.RDS')) # Steady State Params
catch <- readRDS(paste0(path , 'commercial_catch_matrix.RDS')) # Catch
trawl_indices <- readRDS(paste0(path , 'survey_trawl_catch.RDS')) # Trawl
trawl_params <- readRDS(paste0(path , 'trawl_parameters.RDS')) # Trawl Gear Params
trawl_effort <- readRDS(paste0(path , 'survey_trawl_effort.RDS')) # Trawl Effort (1/0 matrix)
seals <- readRDS(paste0(path , 'seal_observations.RDS')) # Seal Observations 

## Set Up Seal Params ----

sealParams <- setSealParams(
  init_params , w_max_seal = 80000 , w_min_seal = 11000 ,
  interaction_seal = c(1 , 1 , 1 , 0 , 0 , 0 , 1 ) ,
  beta = 5000 , sigma = 3 , resource_interaction_seal = 1 ,
  time_step = 35 , h = 250 , 
  dynamicSeals = T , # Below Params only needed if dynamicSeals = T
  sealHarvest = matrix(seals$CatchKg , nrow = 35 , ncol = 1), 
  alpha = matrix(0 , nrow = 35 , ncol = 1) ,
  Msrr = 54.9*80000^-0.25 , 
  Jmax = 89.4*80000^-0.25 , 
  e = 0.89 ,  M = 1e-11
)

idx <- length(sealParams$w) # Put all seals in their average size bin -- there are other ways to consider doing this
sealParams$initialSealN[1,idx] <- seals$N[seals$Year == 1985]/sealParams$dw[idx] # Set Initial Seal Numbers

sealParams$idx <- idx # Overwrite this because you're just using the average seal size: this is important for the repro function

params <- setSeals(init_params , sealParams)

## Set Up Trawl Params (Optional) ----
params <- setTrawlParams(params , trawl_params , trawl_effort)

## Test sim
# sim <- project(params , t_max = 35)

## Write Calibration Functions

# Function to plug new values into params
getNewVals <- function(params , newvals = list()) {
  for(i in names(newvals)) {
    val <- slot(params , i)
    if (is.list(newvals[[i]])) {
      .newvals <- newvals[[i]]
      for (k in names(.newvals)) {
        val[[k]][] <- .newvals[[k]]
      }
    } else val[] <- newvals[[i]]
    slot(params , i , check = F) <- val
  }
  return(params)
}


## Fit Steady State ----

catch_avg_85 <- apply(catch[1:5,] , 2 , mean , simplify = T)
trawl_avg_85 <- apply(trawl_indices[1:5,,] , 2:3 , mean , simplify = T)
effort_avg <- apply(trawl_effort[1:5,] , 2 , mean)
effort_avg <- matrix(effort_avg , nrow = 1)
colnames(effort_avg) <- colnames(trawl_effort)

sealParams <- setSealParams(
  init_params , w_max_seal = 80000 , w_min_seal = 11000 ,
  interaction_seal = c(1 , 1 , 1 , 0 , 0 , 0 , 1 ) ,
  beta = 5000 , sigma = 3 , resource_interaction_seal = 1 ,
  time_step = 35 , h = 250 , 
  dynamicSeals = T , # Below Params only needed if dynamicSeals = T
  sealHarvest = matrix(mean(seals$CatchKg[1:5]) , nrow = 1 , ncol = 1), 
  alpha = matrix(0 , nrow = 35 , ncol = 1) ,
  Msrr = 54.9*80000^-0.25 , 
  Jmax = 89.4*80000^-0.25 , 
  e = 0.89 ,  M = 1e-12 ,
  idx = sealParams$idx
)

idx <- sealParams$idx # Put all seals in their average size bin -- there are other ways to consider doing this
sealParams$initialSealN[1,idx] <- mean(seals$N[1:5]/sealParams$dw[idx]) # Set Initial Seal Numbers
params <- setSeals(init_params , sealParams)
params <- setTrawlParams(params , trawl_params , effort = effort_avg)

getTVRDD2 <- function (params, t, ...) {
  rdd <- getCoeff(params, par = "rdd", where = "other_params", t)[1, ]
  return(rdd)
}

params@other_params$z <- matrix(1 , nrow = 1 , ncol = length(params@species_params$species))
params@other_params$a <- matrix(1 , nrow = 1 , ncol = length(params@species_params$species))
params@other_params$rdd <- matrix(getTVRDD(getRDI(params) , species_params = params@species_params , params = params , 1) , nrow = 1)

params <- setRateFunction(params, "Mort" , "getTVMort")
params <- setRateFunction(params, "RDD" , "getTVRDD2")

source('R/RTMBFunsDef.R')
source('R/RTMBFuns.R')
source('R/dat_to_include.R')

rates_fns <- lapply(RTMBRates , get)

dat <- mizerSeals:::S4toList(params)
dat <- dat[sapply(dat , is.numeric)]
dat$ft_mask <- params@ft_mask

for(i in c('interaction_resource' , 'alpha' , 'erepro' , 'R_max')) {
  dat[[i]] <- params@species_params[[i]]
}

names(sealParams)[!grepl('_seal$' , names(sealParams))] <- paste0(names(sealParams)[!grepl('_seal$' , names(sealParams))] , '_seal')
names(params@trawlParams) <- paste0(names(params@trawlParams) , '_trawl')

dat <- append(dat , params@trawlParams)
dat <- append(dat , sealParams)
dat <- append(dat[names(dat) %in% rates_vars] , dat[grepl('_trawl' , names(dat))])

fun <- function(...) {
  
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  
  pars <- list(...)
  pars <- lapply(pars , exp) # Exponentialize all: when you have a random error for seals this will need to change
  list2env(pars , envir = environment())
  
  n <- n*rescale # Rescale initial N for steady state
  
  no_sp <- nrow(n)
  no_w <- length(w)
  idx <- 2:no_w
  years <- steps*dt # steps should be calculated by length of effort
  
  f_at_dt <- n_at_dt <- array(0 , dim = c(no_sp , no_w , steps))
  s_at_dt <- array(0 , dim = c(steps , length(w_seal)))
  rdi_at_dt <- matrix(0 , ncol = no_sp , nrow = steps)
  rdd_at_dt <- matrix(0 , ncol = no_sp , nrow = steps)
  p_at_dt <- matrix(0 , ncol = no_sp+1 , nrow = steps)
  t_at_dt <- array(0 , dim = c(dim(selectivity_trawl)[1:2] , steps))
  y_at_dt <- matrix(0 , nrow = steps , ncol = no_sp)
  
  w_min_idx_array_ref <- (w_min_idx - 1) * no_sp + (1:no_sp)
  a <- matrix(0, nrow = no_sp, ncol = no_w)
  b <- matrix(0, nrow = no_sp, ncol = no_w)
  S <- matrix(0, nrow = no_sp, ncol = no_w)
  
  dat <- mget(rates_vars , inherits = T)
  
  for (i_time in 1:steps) {
    # Set Effort, Seals RE, etc... 
    
    effort_trawl_dt <- get_val_at_t(effort_trawl , t)
    effort_dt <- get_val_at_t(effort , t)
    alpha_seal_dt <- get_val_at_t(alpha_seal , t)
    z_dt <- get_val_at_t(z , t)
    rdd_dt <- get_val_at_t(rdd , t)
    H_dt <- get_val_at_t(H , t)
    
    r <- rates_fns$Rates(n = n, n_pp = n_pp, n_other = n_other, effort = effort_dt , effort_trawl = effort_trawl_dt , 
                         rdd = rdd_dt , z = z_dt , alpha_seal = alpha_seal_dt , H = H_dt , 
                         dt = dt , dat = dat)
    
    n_pp <- r$n_pp_new
    n_other <- r$seal_repro
    
    # Calculate New biomass
    
    a[, idx] <- sweep(-r$e_growth[, idx - 1, drop = FALSE] *
                        dt, 2, dw[idx], "/")
    b[] <- 1 + sweep(r$e_growth * dt, 2, dw, "/") +
      r$mort * dt
    S[, idx] <- n[, idx, drop = FALSE]
    n[w_min_idx_array_ref] <- (n[w_min_idx_array_ref] + r$rdd *
                                 dt/dw[w_min_idx])/b[w_min_idx_array_ref]
    
    for (i in 1:nrow(n)) {
      for (j in (w_min_idx[i]+1):length(w)) {
        n[i,j] = (S[i,j] - a[i,j]*n[i,j-1]) / b[i,j];
      }
    }
    
    # n <- mizer:::inner_project_loop(no_sp = no_sp, no_w = no_w, n = n,
    #                                 A = a, B = b, S = S, w_min_idx = w_min_idx)
    t <- t + dt
    
    n_at_dt[,,i_time] <- n 
    s_at_dt[i_time ,] <- n_other
    rdi_at_dt[i_time,] <- r$rdi
    rdd_at_dt[i_time,] <- r$rdd
    p_at_dt[i_time,] <- r$seal_diet
    t_at_dt[,,i_time] <- r$t_catch
    y_at_dt[i_time,] <- r$yield
  }  
  
  n_at_y <- n_at_dt[,,seq(1 , i_time , by = years)]
  s_at_y <- s_at_dt[seq(1 , i_time , by = years) , ]
  y <- rep(seq_along(1:years), each = years)
  y_at_y <- t(sapply(1:years , function(x) (apply(y_at_dt[which(y == x),] , 2 , sum))))
  t_at_y <- out$t_at_dt[,,seq(1 , i_time , by = years)]
  biomass <- apply(sweep(n_at_y, 2, w*dw, "*"), c(3, 1), sum)
  
  # Fit Model
  ret <- 0
  
  # Seal Consumption Reports (maybe needed, but can just do it after...)
  
  ## Fitting Functions
  # Fit to catch
  log_diff_f <- (log(c(y_at_y[,obs$catch>0])) - log(rep(obs$catch[obs$catch > 0] , each = years)))
  ret <- ret - sum(dnorm(log_diff_f , 0 , 1 , log = T))
  
  # Fit to trawl (This could also just be to biomass)
  pred_t <- c(t_at_y) ; obs_t <- rep(c(obs$trawl) , years)
  log_diff_t <- log(pred_t[obs_t > 0]) - log(obs_t[obs_t > 0])
  ret <- ret - sum(dnorm(log_diff_t , 0 , 5 , log = T))
  
  # Fit RE F (normal or lognormal, allow for switch)
  # Fit RE Recruitment (Normal or lognormal distribution)
  
  R_max_new = -rdi_at_dt/(1-(rdi_at_dt/rdd_at_dt))
  
  # for(i in (1:ncol(R_max_new))[-5]) {
  #   ret <- ret - sum(dnorm(diff(R_max_new[,i]) , 0 , 1 , log = T))
  # }
  
  #params@species_params$R_max == -getRDI(params)/(1-(getRDI(params)/getRDD(params)))
  # rdd = rdi/(1 + rdi/R_max)
  
  
  # Fit RE Mu_b (Normal distribution, random walk)
  
  # Fit to seal biomass
  log_diff_s <- log(rowSums(s_at_y*params@sealParams$dw)) - log(rep(obs$seals , years))
  ret <- ret - sum(dnorm(log_diff_s , 0 , 1 , log = T))
  
  ## Report Outputs
  REPORT(n_at_y)
  REPORT(t_at_y)
  REPORT(s_at_y)
  REPORT(y_at_y)
  REPORT(rdd_at_dt)
  REPORT(rdi_at_dt)
  REPORT(t_at_dt)
  REPORT(y_at_dt)
  REPORT(n_at_dt)
  
  ## Report Params
  #REPORT(R_max)
  REPORT(effort)
  REPORT(rdd)
  REPORT(z)
  REPORT(alpha_seal)
  REPORT(R_max_new)
  REPORT(biomass)
  
  REPORT(log_diff_s)
  REPORT(log_diff_f)
  REPORT(log_diff_t)
  
  REPORT(r)
  
  # Fit Seal RE (Normal or Lognormal distribution)
  return(ret)
  
}

# Set Data Needed to initialize the function 
data <- list(
  n = params@initial_n , 
  n_pp = params@initial_n_pp, 
  n_other = params@initial_n_other$seals , 
  effort = matrix(params@initial_effort , nrow = 1) ,
  alpha_seal = params@sealParams$alpha ,
  H = params@sealParams$sealHarvest , 
  z = params@other_params$z , 
  rdd = params@other_params$rdd ,
  t = 0, dt = 0.1, steps = 100 ,  
  rates_fns = rates_fns , 
  rates_vars = rates_vars , 
  obs = list(
    catch = catch_avg_85 , 
    trawl = trawl_avg_85 , 
    trawl_effort = effort_avg , 
    seals = mean(seals$N[1:5]))
)
data <- append(data , dat)
data$effort_trawl <- matrix(data$effort_trawl , nrow = 1)

params@catchability[6,6] <- 0

pars <- list(
  rdd = params@other_params$rdd , 
  z = params@other_params$z , 
  catchability_trawl = params@trawlParams$catchability, 
  M_seal = params@sealParams$M ,
  resource_interaction_seal = params@sealParams$resource_interaction_seal ,
  catchability = params@catchability ,
  rescale = rep(1 , length(params@species_params$species)) , 
  effort = matrix(params@initial_effort , nrow = 1)
)

log_pars <- lapply(pars , log)
data <- data[!names(data) %in% names(log_pars)]
environment(fun) <- list2env(data)

# out <- do.call(fun , log_pars)

map <- with(
  log_pars ,
  list(
    catchability_trawl = matrix(
      factor(1:length(catchability_trawl) , levels = which(!is.infinite(catchability_trawl))) ,
      nrow = nrow(catchability_trawl) ,
      ncol = ncol(catchability_trawl)
    ) ,
    effort = as.character(factor((1:length(effort)) , levels = (which(!is.infinite(effort))))),
    catchability = matrix(
      factor(1:length(catchability) , levels = which(!is.infinite(catchability))) ,
      nrow = nrow(catchability) ,
      ncol = ncol(catchability)
    ) ,
    #catchability = as.character(factor((1:length(catchability)) , levels = which(params@species_params$species != 'SNOWCRAB.F'))),
    rdd = species_params(params)$main_spec
  )
)


for(i in 1:length(map)) map[[i]] <- factor(map[[i]])

obj <- MakeADFun(function(p)do.call(fun,p), log_pars , map = map)
system.time(opt <- nlminb(obj$par, obj$fn, obj$gr , control = list(iter.max = 10000 , eval.max = 2000 , trace = 0)))

opt$message
cbind(names(opt$par) , obj$gr()[1,])

opt$par

trawl <- obj$report()$t_at_y
yield <- obj$report()$y_at_y
seals <- obj$report()$s_at_y

trawl <- `dimnames<-`(trawl , list(gear = dimnames(obs$trawl)$gear , sp = dimnames(obs$trawl)$sp))

plot(yield[1,] , obs$catch)

trawl_comp <- merge(array2DF(trawl) , array2DF(obs$trawl) , by = c('gear' , 'sp'))
trawl_comp <- trawl_comp[trawl_comp$Value.y > 0,]

plot(trawl_comp$Value.x , trawl_comp$Value.y)
