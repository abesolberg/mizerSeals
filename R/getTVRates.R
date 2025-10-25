getCoeff <- function (params, t , ...) {
  arr <- params@sealParams$sealHarvest
  if(is.null(rownames(arr))) rownames(arr) <- 1:nrow(arr)
  if (!floor(t) %in% as.numeric(rownames(arr)))
    t <- t%%as.numeric(rownames(arr)[nrow(arr)]) +
      as.numeric(rownames(arr)[1])
  index_vec <- as.numeric(rownames(arr)) - t
  index <- which(index_vec == max(index_vec[index_vec <= 0]))
  #cat(n , '--', t , '--' , index , '\n')
  h_at_t <- arr[index,,drop=F]
  return(h_at_t)
}


## Time Vaying Mortality

## Recruitment Anomalies