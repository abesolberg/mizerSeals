# Get TV Rate
getCoeff <- function (params, par , where = 'other_params' , t , ...) {
  arr <- slot(params , where)[[par]]
  if(is.null(rownames(arr))) rownames(arr) <- 1:nrow(arr)
  if (!floor(t) %in% as.numeric(rownames(arr)))
    t <- t%%as.numeric(rownames(arr)[nrow(arr)]) +
      as.numeric(rownames(arr)[1])
  index_vec <- as.numeric(rownames(arr)) - t
  index <- which(index_vec == max(index_vec[index_vec <= 0]))
  #cat(n , '--', t , '--' , index , '\n')
  var_at_t <- arr[index,,drop=F]
  return(var_at_t)
}


## Time Vaying Mortality
getTVMort <- function (params, n, n_pp, n_other, t, f_mort, pred_mort, ...) {
  
  mort <- pred_mort + (params@mu_b*getCoeff(params , par = 'z' , where = 'other_params', t)[1,]) + f_mort
  for (i in seq_along(params@other_mort)) {
    mort <- mort + do.call(params@other_mort[[i]], list(params = params, 
                                                        n = n, n_pp = n_pp, n_other = n_other, t = t, component = names(params@other_mort)[[i]], 
                                                        ...))
  }
  return(mort)
}


## Recruitment Anomalies
getTVRDD <- function(rdi , species_params , params , t , ... ) {
  if (!("R_max" %in% names(species_params))) {
    stop("The R_max column is missing in species_params.")
  }
  if (!("sex" %in% names(species_params))) {
    stop("The sex column is missing in species_params.")
  }
  if (!("main_spec" %in% names(species_params))) {
    stop("The main_spec column is missing in species_params.")
  }
  if (!("sex_ratio" %in% names(species_params))) {
    stop("The sex_ratio column is missing in species_params.")
  }
  rdd <- (rdi / (1 + rdi / species_params$R_max))*getCoeff(params , par = 'a' , where = 'other_params' , t)[1,]
  rdd <- rdd * species_params$sex_ratio*2
  
  if(any(species_params$sex=="female")){
    rdd_sum<-tapply(rdd, list(species_params$main_spec,species_params$sex) , sum)
    rdd_sum[is.na(rdd_sum)]<-0
    
    ma<-match(species_params$main_spec,row.names(rdd_sum))
    rdd_split<-rdd_sum[ma,"female"] / 2
    
    ma<- match(names(rdd_split), species_params$main_spec)
    ord<-rdd_split[ma]
    rdd[!ord==0]<-ord[!ord==0]
    rdd[is.na(rdd)]<-0 
  }
  
  return(rdd)
}

# Fully Independent Recruitment
getTVRDD2 <- function (params, t, ...) {
  rdd <- getCoeff(params, par = "rdd", where = "other_params", t)[1, ]
  return(rdd)
}
