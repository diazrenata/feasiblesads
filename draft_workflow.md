draft\_workflow
================
Renata Diaz
4/8/2019

Define S and N
--------------

``` r
S = 5
N = 30
```

Simulate an "empirical" RAD
---------------------------

``` r
empirical_rad <- sample_fs(nsamples = 1, nspec = S, nind = N)[[2]] %>%
  t() %>%
  as.vector()
empirical_rad
```

    ## [1]  3  5  5  6 11

Sample from the feasible set given S and N
------------------------------------------

``` r
sim_rads <- sample_fs(nsamples = 1000, nspec = S, nind = N)[[2]] %>%
  t() %>%
  as.data.frame() %>%
  dplyr::distinct()
head(sim_rads)
```

    ##   V1 V2 V3 V4 V5
    ## 1  3  5  5  6 11
    ## 2  4  4  5  8  9
    ## 3  3  5  6  7  9
    ## 4  4  5  6  7  8
    ## 5  4  5  7  7  7
    ## 6  3  4  5  7 11

``` r
dim(sim_rads)
```

    ## [1] 108   5

Make sure the empirical RAD is in the sampled set
-------------------------------------------------

``` r
sim_rads <- rbind(sim_rads, empirical_rad) %>%
    dplyr::distinct()

dim(sim_rads)
```

    ## [1] 108   5

Identify the central tendency of the (samples of the) feasible set
------------------------------------------------------------------

### Using Locey & White's approach

``` r
lw_ct <- find_ct(sim_rads)
lw_ct
```

    ##    V1 V2 V3 V4 V5
    ## 34  3  4  6  8  9

### Using new criterion

``` r
all_rad_ps <- vector(length = nrow(sim_rads), mode = 'numeric')

for(i in 1:nrow(sim_rads)) {
  all_rad_ps[i] <- new_ct_criterion(focal_RAD = as.integer(sim_rads[i,]), 
                                    pool_RADs = sim_rads)[[1]]
  }

new_ct <- sim_rads[which(all_rad_ps == max(all_rad_ps)),]
new_ct
```

    ##    V1 V2 V3 V4 V5
    ## 34  3  4  6  8  9

Evaluate the overall likelihood of the empirical RAD compared to the feasible set
---------------------------------------------------------------------------------

``` r
empirical_p <- new_ct_criterion(focal_RAD = empirical_rad, 
                                pool_RADs = sim_rads)[[1]]

sim_ecdf <- ecdf(all_rad_ps)
sim_ecdf(empirical_p)
```

    ## [1] 0.6481481

``` r
p_plot <- ggplot(data = as.data.frame(all_rad_ps)) + 
  geom_histogram(binwidth = 0.25, aes(x = as.data.frame(all_rad_ps)$all_rad_ps)) +
    xlim(min(all_rad_ps) - 1, max(all_rad_ps) +1) + 
    geom_vline(xintercept = empirical_p, colour  ='red') +
  theme_bw()

p_plot
```

    ## Warning: Removed 2 rows containing missing values (geom_bar).

![](draft_workflow_files/figure-markdown_github/empirical%20overall%20likelihood-1.png)

Evaluate the empirical likelihood at each rank
----------------------------------------------

``` r
all_rad_all_ps <- matrix(nrow = nrow(sim_rads), ncol = ncol(sim_rads))
colnames(all_rad_all_ps) <- c(1:5)

for(i in 1:nrow(sim_rads)) {
all_rad_all_ps[i, ] <-new_ct_criterion(focal_RAD = as.integer(sim_rads[i,]), pool_RADs = sim_rads)[[2]]
}

all_rad_all_ps <- as.data.frame(all_rad_all_ps) %>%
  stack()

empirical_all_ps <- new_ct_criterion(focal_RAD = empirical_rad,
                                     pool_RADs = sim_rads)[[2]]

rank_p_plot <-ggplot(data = all_rad_all_ps) + 
  geom_jitter(data = all_rad_all_ps, aes(x = all_rad_all_ps$ind, y = all_rad_all_ps$values)) + 
  geom_point(data = as.data.frame(empirical_all_ps), aes(x = c(1:S), y =  as.data.frame(empirical_all_ps)$empirical_all_ps), colour= 'red') +
  theme_bw()

rank_p_plot
```

![](draft_workflow_files/figure-markdown_github/empirical%20rank%20likelihood-1.png)

``` r
empirical_p_quantiles <- vector(length = S)

for(i in 1:S){
  this_rank_ecdf <- ecdf(dplyr::filter(all_rad_all_ps, ind == i)$values)
  empirical_p_quantiles[i] <- this_rank_ecdf(empirical_all_ps[i])
}

empirical_rank_plot <- ggplot(data = as.data.frame(empirical_p_quantiles)) +
  geom_point(data = as.data.frame(empirical_p_quantiles), aes(x = c(1:S), y =  as.data.frame(empirical_p_quantiles)$empirical_p_quantiles)) +
  ylim(c(0, 1)) + 
  theme_bw()
empirical_rank_plot
```

![](draft_workflow_files/figure-markdown_github/empirical%20rank%20likelihood-2.png)
