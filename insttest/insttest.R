library(mizer)
library(mizerSeals)
library(ggplot2)

# https://publications.gc.ca/collections/collection_2020/mpo-dfo/fs70-5/Fs70-5-2020-014-eng.pdf
# https://www.frontiersin.org/journals/marine-science/articles/10.3389/fmars.2021.579946/full

# sapply(list.files('R/',full.names=T) , source)

params <- mizer::NS_params

sealParams <- setSealParams(
  params , w_max_seal = 150000 , w_min_seal = 11000 ,
  interaction_seal = rep(1 , length(params@species_params$species)) ,
  beta = 1000 , sigma = 3 ,
  time_step = 1 , dynamicSeals = F
)
sealParams$initialSealN[,107] <- rnorm(1 , 20000000 , 250000)/sealParams$dw[107]
#sealParams$initialSealN[,107] <- 1:100
params <- setSeals(params , sealParams)
sim <- project(params , t_max = 5)
sim2 <- project(params, t_max= 5 , t_save = 0.1)

dim(sim2@n_other)

dim(sim@n_other)

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



