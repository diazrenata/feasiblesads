#' @title Build the k table
#'
#' @description Build the table of the maximum abundance for the minimum-abundance species in a community of s species and n total individuals (or the maximum increase in abundance compared to the past species, when building an SAD using a uniform sampler...)
#'
#' @param max_s the max number of species
#' @param max_n the max number of individuals
#' @return table of ks for all combinations of s and n up to max_s and max_n.
#'
#' @export

fill_ks <- function(max_s, max_n){
  # ks are the possible values for a given n given number of species and remaining abundance
  kmax <- matrix(ncol = max_n + 1, nrow = max_s)

  colnames(kmax) <- c(0:max_n)

  for(i in 1:nrow(kmax)) {
    for(j in 1:ncol(kmax)) {
      kmax[i, j] <- floor((j-1)/i)
    }
  }

  return(kmax)

}

#' @title Build the partition counts table.
#'
#' @description Recursively build P(s, n) table where P(s, n) is the number of partitions of n into s, including 0s.
#'
#' @param max_s max species
#' @param max_n max individuals
#' @param storeyn TRUE/FALSE whether to store the P table. If FALSE, returns the table.
#' @param storepath path to store the P table.
#' @return matrix of P(s, n) for all combinations of and s and n up to max_s and max_n.
#'
#' @export

fill_ps <- function(max_s, max_n, storeyn = TRUE,
                    storepath = "uniform_fs_sampling") {

ps <- matrix(ncol = max_n + 1, nrow = max_s)

# columns are individuals, rows are species
colnames(ps) <- c(0:max_n)

ps[1, ] = 1 # only one way to put any number of individuals into one species
ps[, "0"] = 1 # only one way to put 0 individuals into any number of species
ps[, "1"] = 1 # only one way to put 1 individual into any number of species

kmax = fill_ks(max_s = max_s, max_n = max_n)

p_lookup <- function(k, s, n) {
  this_p = ps[s-1, n - (s*k) + 1]
  return(gmp::as.bigz(this_p))
}

p_over_k <- function(s, n) {

  sum_p = 0

   this_kmax = kmax[s, n + 1]

  for(k in 0:this_kmax) {
  sum_p = gmp::add.bigz(gmp::as.bigz(sum_p), p_lookup(k, s, n))
  }
   return(sum_p)
}

for(i in 2:nrow(ps)) {
  for(j in 3:ncol(ps)) {
    ps[i, j] <- as.character(p_over_k(i, j - 1))
  }
}

if(storeyn) {
  save(ps, file = file.path(storepath, "p_table.Rds"))
}
return(ps)

}
