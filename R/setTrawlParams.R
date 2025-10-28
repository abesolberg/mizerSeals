# Setup Trawl Params
setTrawlParams <- function(params , gear_params , effort) {
  tmp <- params
  mizer::gear_params(tmp) <- gear_params
  gp <- list(
    selectivity = tmp@selectivity ,
    catchability = tmp@catchability ,
    effort = effort[,dimnames(tmp@selectivity)$gear]
  )
  params <- `slot<-`(params , 'trawlParams' , check = F , gp)
  return(params)
}
