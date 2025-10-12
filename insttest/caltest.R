library(mizer)
library(mizerSeals)
library(ggplot2)

# sapply(list.files('R/',full.names=T) , source)

params <- mizer::NS_params

sim1 <- project(params , t_max = 10 , dt = 0.1 , t_save = 1)
sim2 <- project(params, t_max = 10 , dt = 0.1 , t_save = 0.1)
sim3 <- project(params , t_max = 10 , dt = 1)

y1 <- getYield(sim1)
y2 <- getYield(sim2)
y3 <- getYield(sim3)

ggplot() +
  geom_line(aes(x = 0:10 , y = y1[,'Cod'])) +
  geom_line(aes(x = as.numeric(rownames(y2)) , y = y2[,'Cod'])) +
  geom_line(aes(x = 0:10 , y = y3[,'Cod']) , linetype = 'dashed') +
  scale_x_continuous(breaks = 0:10) +
  geom_abline()

ggplot() +
  geom_point(aes(
    x = rep(y1[, 'Cod'] , each = 1) ,
    y = y2[seq(1, 101 , by = 10), 'Cod'] ,
    color = factor(floor(as.numeric(rownames(
      y2
    )[seq(1, 101 , by = 10)])))
  )) +
  geom_point(aes(x = y1[-1, 'Cod'] , y = y1[-1, 'Cod']) , color = 'red' , size = 3) +
  geom_abline() + theme(legend.position = 'bottom' , legend.title = element_blank())

ggplot() +
  geom_point(aes(
    x = rep(y1[2:11, 'Cod'] , each = 10) ,
    y = y2[2:101, 'Cod'] ,
    color = factor(floor(as.numeric(rownames(
      y2
    )[2:101])))
  )) +
  geom_point(aes(x = y1[-1, 'Cod'] , y = y1[-1, 'Cod']) , color = 'red' , size = 3) +
  geom_abline() + theme(legend.position = 'bottom' , legend.title = element_blank())

b1 <- getBiomass(sim1)
b2 <- getBiomass(sim2)
b3 <- getBiomass(sim3)

ggplot() +
  geom_line(aes(x = 0:10 , y = b1[,'Cod'])) +
  geom_line(aes(x = as.numeric(rownames(b2)) , y = b2[,'Cod'])) +
  geom_line(aes(x = 0:10 , y = b3[,'Cod']) , linetype = 'dashed') +
  scale_x_continuous(breaks = 0:10) +
  geom_abline()

ggplot() +
  geom_point(aes(
    x = rep(b1[, 'Cod'] , each = 1) ,
    y = b2[seq(1, 101 , by = 10), 'Cod'] ,
    color = factor(floor(as.numeric(rownames(
      b2
    )[seq(1, 101 , by = 10)])))
  )) +
  geom_point(aes(x = b1[-1, 'Cod'] , y = b1[-1, 'Cod']) , color = 'red' , size = 3) +
  geom_abline() + theme(legend.position = 'bottom' , legend.title = element_blank())

ggplot() +
  geom_point(aes(
    x = rep(b1[2:11, 'Cod'] , each = 10) ,
    y = b2[2:101, 'Cod'] ,
    color = factor(floor(as.numeric(rownames(
      b2
    )[2:101])))
  )) +
  geom_point(aes(x = b1[-1, 'Cod'] , y = b1[-1, 'Cod']) , color = 'red' , size = 3) +
  geom_abline() + theme(legend.position = 'bottom' , legend.title = element_blank())

getDiet
