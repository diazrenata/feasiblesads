library(feasiblesads)

test <- sample_fs(nsamples = 100000, nspec = 16, nind = 4362)

system.time(sample_fs(nsamples = 100000, nspec = 16, nind = 4362))

sads <- test[[2]]

sads <- t(sads)

sads <- as.data.frame(sads)
sads <- dplyr::distinct(sads)
