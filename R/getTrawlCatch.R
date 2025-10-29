# Fit to trawl survey
mizerTMortGear <- function (params , effort) {
  out <- params@trawl_params$selectivity
  out[] <- effort * c(params@trawl_params$catchability) * c(params@trawl_params$selectivity)
  return(out)
}

getTrawlMort <- function (object, effort, time_range) {
  if (is(object, "MizerSim")) {
    sim <- object
    if (missing(time_range)) {
      time_range <- dimnames(sim@effort)$time
    }
    if (nrow(sim@params@trawl_params$effort) != nrow(sim@effort)) {
      missing <- nrow(sim@effort) - nrow(sim@params@trawl_params$effort)
      if(!is.matrix(sim@params@trawl_params$trawl_effort))  sim@params@trawl_params$trawl_effort <- matrix(sim@params@trawl_params$trawl_effort , nrow = 1)
      sim@params@trawl_params$trawl_effort <- cbind(
        matrix(
          sim@params@trawl_params$trawl_effort ,
          nrow = missing ,
          ncol =  ncol(sim@params@trawl_params$trawl_effort) ,
          byrow = T)
      )
    }
    time_elements <- get_time_elements(sim, time_range)
    t_mort_gear <- getTrawlMort(sim@params, sim@params@trawl_params$trawl_effort)
    return(t_mort_gear[time_elements, , , , drop = FALSE])
  }
  else {
    params <- validParams(object)
    if (missing(effort)) {
      effort <- rep(1 , dim(params@trawl_params$selectivity)[1])
    }
    if (is(effort, "numeric")) {
      no_gear <- dim(params@trawl_params$catchability)[1]
      if (length(effort) == 1) {
        effort <- rep(effort, no_gear)
      }
      if (length(effort) != no_gear) {
        stop("Effort must be a single value or a vector as long as the number of gears\n")
      }
      cat(effort , '\n')
      f <- mizerTMortGear(params, effort = effort)
      dimnames(f) <- dimnames(params@trawl_params$selectivity)
      return(f)
    }
    else {
      no_gear <- dim(params@trawl_params$catchability)[1]
      if (dim(effort)[2] != no_gear)
        stop("Effort array must have a single value or a vector as long as the number of gears for each time step\n")
      out <- array(NA, dim = c(dim(params@trawl_params$selectivity),
                               dim(effort)[1]), dimnames = c(dimnames(params@trawl_params$selectivity),
                                                             list(time = dimnames(effort)[[1]])))
      out[] <- apply(effort, 1, function(x) mizerTMortGear(params,
                                                           x))
      out <- aperm(out, c(4, 1, 2, 3))
      return(out)
    }
  }
}

getTrawlCatch <- function(object) {
    if (is(object, "MizerSim")) {
      sim <- object
      biomass <- sweep(sim@n, 3, sim@params@w * sim@params@dw,
                       "*")
      f <- getTrawlMort(sim)
      return(apply(sweep(f, c(1, 3, 4), biomass, "*"),
                   c(1, 2, 3), sum))
    }
    if (is(object, "MizerParams")) {
      params <- object
      biomass <- sweep(params@initial_n, 2, params@w * params@dw,
                       "*")
      f <- getTrawlMort(params)
      return(apply(sweep(f, c(2, 3), biomass, "*"), c(1,
                                                           2), sum))
    }
    stop("'object' should be a MizerParams or a MizerSim object")
}

getTrawlCatchSize <- function(object) {
  if (is(object, "MizerSim")) {
    sim <- object
    biomass <- sweep(sim@n, 3, sim@params@w * sim@params@dw,
                     "*")
    f <- getTrawlMort(sim)
    return(sweep(f, c(1, 3, 4), biomass, "*"))
  }
  if (is(object, "MizerParams")) {
    params <- object
    biomass <- sweep(params@initial_n, 2, params@w * params@dw,
                     "*")
    f <- getTrawlMort(params)
    return(sweep(f, c(2, 3), biomass, "*"))
  }
  stop("'object' should be a MizerParams or a MizerSim object")
}

