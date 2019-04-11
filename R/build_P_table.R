# Recursively build P(s, n) table where P(s, n) is the number of partitions of N into S, including 0s.

library(bit64)

max_s = 20
max_n = 100

fill_ps <- function(max_s, max_n) {

ps <- matrix(ncol = max_n + 1, nrow = max_s)

# columns are individuals, rows are species

colnames(ps) <- c(0:max_n)

ps[1, ] = 1 # only one way to put any number of individuals into one species
ps[, "0"] = 1 # only one way to put 0 individuals into any number of species
ps[, "1"] = 1 # only one way to put 1 individual into any number of species

# ks are the possible values for a given n given number of species and remaining abundance
kmax <- matrix(ncol = max_n + 1, nrow = max_s)

colnames(kmax) <- c(0:max_n)

for(i in 1:nrow(kmax)) {
  for(j in 1:ncol(kmax)) {
    kmax[i, j] <- floor((j-1)/i)
  }
}

p_lookup <- function(k, s, n) {
  this_p = ps[s-1, n - (s*k) + 1]
  return(as.integer(this_p))
}

p_over_k <- function(s, n) {

  sum_p = 0

   this_kmax = kmax[s, n + 1]

  for(k in 0:this_kmax) {
  sum_p = as.integer64(sum_p + p_lookup(k, s, n))
  }
   return(sum_p)
}

for(i in 3:ncol(ps)) {
  for(j in 2:nrow(ps)) {
    ps[j, i] <- p_over_k(j, i - 1)
  }
}

return(ps)

}

system.time(fill_ps(20, 1000))
