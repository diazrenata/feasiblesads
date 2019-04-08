library(feasiblesads)

test <- sample_fs(nsamples = 500, nspec = 5, nind = 30)

system.time(sample_fs(nsamples = 100000, nspec = 16, nind = 4362))

sads <- test[[2]]

sads <- t(sads)

sads <- as.data.frame(sads)
sads <- dplyr::distinct(sads)

toy_empirical <- sample_fs(nsamples = 1, nspec= 5, nind = 30)[[2]]

pool_RADs <- sads

pool_RADs[nrow(pool_RADs) + 1,] <- toy_empirical

pool_RADs <- dplyr::distinct(pool_RADs)

par(mfrow = c(3,ceiling(length(toy_empirical)/3)))
for(i in 1:length(toy_empirical)) {
  hist(pool_RADs[,i])
  #== add observed likelihood
  abline(v=toy_empirical[i],lty=2, col = 'red')
}


# what is the most likely rad given this criterion?
all_RAD_ps <- vector(length = nrow(pool_RADs), mode = 'numeric')
all_RAD_all_ps <- matrix(nrow = nrow(pool_RADs), ncol = ncol(pool_RADs))

for(i in 1:nrow(pool_RADs)) {
  all_RAD_ps[i] <- new_ct_criterion(focal_RAD = as.integer(pool_RADs[i,]), pool_RADs = pool_RADs)[[1]]
all_RAD_all_ps[i, ] <-new_ct_criterion(focal_RAD = as.integer(pool_RADs[i,]), pool_RADs = pool_RADs)[[2]]

  }


empirical_p <- new_ct_criterion(focal_RAD = toy_empirical, pool_RADs = pool_RADs)[[1]]

which(all_RAD_ps== max(all_RAD_ps))
pool_RADs[which(all_RAD_ps== max(all_RAD_ps)),]

all_RAD_all_ps_quantiles <- matrix(nrow= nrow(all_RAD_all_ps), ncol = ncol(all_RAD_all_ps))

for(i in 1:nrow(all_RAD_all_ps_quantiles)) {
  for(j in 1:ncol(all_RAD_all_ps_quantiles)) {
    this_ecdf <-  ecdf(all_RAD_all_ps[,j])
    all_RAD_all_ps_quantiles[i, j] <-this_ecdf(all_RAD_all_ps[i, j])
  }
}


hist(all_RAD_ps)
#== add observed p
abline(v=empirical_p,lty=2, col = 'red')

sads_list <- list()

for(i in 1:nrow(pool_RADs)) {
  sads_list[[i]] <- as.numeric(pool_RADs[i, ])
}

ct = find_ct(sads_list)


empirical_ps <- new_ct_criterion(focal_RAD = toy_empirical, pool_RADs = pool_RADs)[[2]]
