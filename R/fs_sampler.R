#' @title Draw a sample from the feasible set given S and N
#'
#' @description Sample the feasible set.
#'
#' @param s how many species
#' @param n how many individuals
#' @return a single sample from the feasible set for the SAD
#'
#' @export

draw_fs <- function(i, nspec, nind){

  vect <- sample.int(nspec, size = nind - nspec, replace = T)

  return(c(1:nspec, vect))
}


#' @title Draw multiple samples from the feasible set given S and N
#'
#' @description Sample the feasible set.
#'
#' @param s how many species
#' @param n how many individuals
#' @return matrix of samples
#'
#' @export


sample_fs <- function(nsamples, nspec, nind){

  if(nspec > nind) return(NULL)

  if(nind == nspec) {
    vects <- c(1:nspec)
    sads <- as.matrix(table(vects))
    return(list(vects, sads))
  }

  #create cluster
  cl <- parallel::makeCluster(parallel::detectCores()-1)
  # get library support needed to run the code
  parallel::clusterEvalQ(cl,library(MASS))
  # Set a different seed on each member of the cluster (just in case)
  parallel::clusterSetRNGStream(cl)
  #... then parallel replicate...
  vects <- parallel::parSapply(cl, 1:nsamples, FUN = draw_fs, nind = nind, nspec = nspec)
  sads <- parallel::parApply(cl,X =  vects, MARGIN = 2, FUN = table)
  sads <- parallel::parApply(cl, X = sads, MARGIN = 2, FUN = sort)
   #stop the cluster
  parallel::stopCluster(cl)

  return(list(vects, sads))
}

