library(feasiblesads)

test <- sample_fs(nsamples = 500, nspec = 5, nind = 30)

system.time(sample_fs(nsamples = 100000, nspec = 16, nind = 4362))

sads <- test[[2]]

sads <- t(sads)

sads <- as.data.frame(sads)
sads <- dplyr::distinct(sads)


sads_list <- list()

for(i in 1:nrow(sads)) {
  sads_list[[i]] <- as.numeric(sads[i, ])
}

ct = find_ct(sads_list)
