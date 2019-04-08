#' @title New criterion for the likelihood of an RAD
#'
#' @description Testing a new criterion
#'
#' @param focal_RAD focal RAD
#' @param pool_RADs pool to compare to, must have same S and N as `focal_RAD`
#' @return list of log p's and their sum
#'
#' @export

new_ct_criterion <- function(focal_RAD, pool_RADs){

  pool_RADs[nrow(pool_RADs) + 1,] <- focal_RAD

  pool_RADs <- dplyr::distinct(pool_RADs)

  ps = vector(length = length(focal_RAD), mode= 'numeric')

  for(i in 1:length(ps)) {
    ps[i] <- pull_discrete_p(value = focal_RAD[i], pool = as.vector(pool_RADs[,i]))
  }

  overall_p = sum(ps)

  return(list(overall_p, ps))
}

pull_discrete_p <- function(value, pool, log = TRUE){
  p = sum(pool == value) / length(pool)
  if(log) p = log(p)
  return(p)
}
