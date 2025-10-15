# plot Seal Mort

plotSealMort <- function(object , return_data = T , only_mature = T , include = c('F' , 'Ext.' , 'Pred' , 'Seals')) {
  include <- match.arg(include , several.ok = T)
  if (is(object, "MizerSim")) {
    sim <- object
    
    f <- getFMort(sim)
    p <- getPredMort(sim)
    s <- getSealMortSim(sim)
    mu  <- aperm(array(rep(sim@params@mu_b , dim(f)[1]) , dim = c(dim(sim@params@mu_b) , dim(f)[1])), c(3,1,2))
    dimnames(mu) <- dimnames(p) <- dimnames(s) <- dimnames(f)
    dat <- do.call(rbind , list(array2DF(f) , array2DF(p) , array2DF(s) , array2DF(mu)))
    dat$mort <- rep(c('F' , 'Pred' , 'Seals' , 'Ext.') , each = length(f))
    dat$w <- as.numeric(dat$w)
    dat <- dat[with(dat , order(time , sp , mort , w , Value)), ]
    dat <- `rownames<-`(dat , NULL)
    
    n <- array2DF(sim@n)
    n$w <- as.numeric(n$w)
    n <- n[with(n , order(time , sp , w , Value)), ]
    n$dw <- rep(sim@params@dw , nrow(n)/length(sim@params@dw))
    n$N <- n$Value*n$dw
    dat <- merge(dat , n[, c('time' , 'sp' , 'w' , 'N')])
    dat <- merge(dat , sim@params@species_params[,c('species' , 'w_mat')] , by.x = 'sp' , by.y = 'species')
    dat$tmp <- dat$N*dat$Value
    df <- dat[,c('sp','time','w','mort','Value')]
    
    if(only_mature) dat <- dat[ dat$w  >= dat$w_mat,]
    
    dat_p <- aggregate(. ~ time + sp + mort , FUN = sum , data = dat[,c('time' , 'sp' , 'Value' , 'mort' , 'N' , 'tmp')] , na.rm = T)
    dat_p$M <- dat_p$tmp/dat_p$N
    dat_p <- dat_p[dat_p$mort %in% include ,]
   
    print(ggplot2::ggplot(dat_p , ggplot2::aes(x = time , y = M , group = paste(sp , mort) , color = mort)) +
      ggplot2::geom_line() + ggplot2::facet_wrap(~sp , scales = 'free'))
  if(return_data) return(df)
  } else if (is(object, "MizerParams")) {
    params <- object

    f <- getFMort(params)
    p <- getPredMort(params)
    s <- getSealMort(params)
    mu  <- params@mu_b
    dimnames(mu) <- dimnames(p) <- dimnames(s) <- dimnames(f)
    dat <- do.call(rbind , list(array2DF(f) , array2DF(p) , array2DF(s) , array2DF(mu)))
    dat$mort <- rep(c('F' , 'Pred' , 'Seals' , 'Ext.') , each = length(f))
    dat$w <- as.numeric(dat$w)
    dat <- dat[with(dat , order(sp , mort , w , Value)), ]
    dat <- `rownames<-`(dat , NULL)
    
    n <- array2DF(params@initial_n)
    n$w <- as.numeric(n$w)
    n <- n[with(n , order(sp , w , Value)), ]
    n$dw <- rep(sim@params@dw , nrow(n)/length(sim@params@dw))
    n$N <- n$Value*n$dw
    dat <- merge(dat , n[, c( 'sp' , 'w' , 'N')])
    dat <- merge(dat , sim@params@species_params[,c('species' , 'w_mat')] , by.x = 'sp' , by.y = 'species')
    dat$tmp <- dat$N*dat$Value
    df <- dat[,c('sp','w','mort','Value')]
    
    if(only_mature) dat <- dat[ dat$w  >= dat$w_mat,]
    dat_p <- dat[dat$mort %in% include ,]
    print(ggplot2::ggplot(dat , ggplot2::aes(x = w , y = Value , group = paste(sp , mort) , color = mort)) +
            ggplot2::geom_line() + ggplot2::facet_wrap(~sp , scales = 'free') + ggplot2::scale_x_log10())
    if(return_data) return(df)
  } else stop("'object' should be a MizerParams or a MizerSim object")
}


