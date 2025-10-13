library(mizer)
library(mizerSeals)
library(ggplot2)

# https://publications.gc.ca/collections/collection_2020/mpo-dfo/fs70-5/Fs70-5-2020-014-eng.pdf
# https://www.frontiersin.org/journals/marine-science/articles/10.3389/fmars.2021.579946/full

# sapply(list.files('R/',full.names=T) , source)

params <- mizer::NS_params

sealParams <- setSealParams(
  mizer::NS_params , w_max_seal = 150000 , w_min_seal = 11000 ,
  interaction_seal = rep(1 , length(params@species_params$species)) ,
  beta = 5000 , sigma = 3 , resource_interaction_seal = 100 ,
  time_step = 1 , dynamicSeals = F
)
sealParams$initialSealN[,107] <- rnorm(1 , 20000000 , 0)/sealParams$dw[107]
#sealParams$initialSealN[,107] <- 1:100
params <- setSeals(mizer::NS_params , sealParams)
sim <- project(params , t_max = 100)

plot(sim)

diet <- getSealDiet(params)
diet$diet_by_sp
diet$percent_by_sp
mort <- getSealMort(params)
resource_mort <- getSealResourceMort(params)

plot(resource_mort)

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



