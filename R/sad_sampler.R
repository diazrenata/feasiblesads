#' @title Draw samples from the feasible set given S and N
#' @description Draw samples from the feasible set of species abundance distributions with `s` species and `n` total individuals
#' @param s the number of species
#' @param n the number of individuals
#' @param nsamples how many samples to draw
#' @param storeyn TRUE/FALSE whether to store the P table. If FALSE, returns the table.
#' @param storepath path to store the P table.
#' @return matrix `nsamples` of ranked abundance distributions. Rows are samples and columns are species, sorted from least to most abundant.
#' @export

sample_fs = function(s, n, nsamples, storeyn = FALSE,
                     storepath = NULL)
{
    ps <- fill_ps(s, n, storeyn = storeyn, storepath = storepath)
    
    # Once you have the ps table you also might as well make a ks table? looking up might be as slow as calculating, but whatever.
    ks <- fill_ks(s, n)
    
    sets <- matrix(NA, nrow = nsamples, ncol = s)
    for(idx in seq_len(nsamples))
    {
        
        this_gnome = vector(length = s, mode = 'integer')
        slots_remaining = s
        n_remaining = n
        
        for(species_slot in 1:s) {
            
            if(s == 1) { # if s = 1 there is only one possible RAD.
                this_gnome[1] = n
            } else {
                if(species_slot == 1) {     # first (lowest abundance) value is a boundary case
                    # for the first one, the values of n range from 1 to ks[s, n+1]
                    ns = 1:ks[s, n+1]
                    
                    # the probability that you select each possible n is
                    # proportional to the # of partitions possible if you
                    # were to choose that n.
                    # divided by the total possible # of partitions across all choices of n.
                    n_parts = vector(length = length(ns), mode = "character")
                    p_ns = vector(length = length(ns), mode = "numeric")
                    for(i in 1:length(n_parts)) {
                        n_parts[i] = ps[s-1, n - (s * ns[i]) + 1]
                    }
                    # total partition is
                    total_parts = 0
                    for(i in 1:length(n_parts)) {
                        total_parts = gmp::add.bigz(gmp::as.bigz(n_parts[i]), total_parts)
                    }
                    # probability of each n choice
                    for(i in 1:length(n_parts)){
                        p_ns[i] = as.numeric(gmp::as.bigq(gmp::as.bigz(n_parts[i]), total_parts))
                    }
                    
                    this_gnome[species_slot] = sample(ns, size = 1, prob = p_ns)
                    n_remaining = n_remaining - (s*this_gnome[species_slot])
                    slots_remaining = slots_remaining - 1
                    
                } else if(species_slot < s){ # recursive probabilistic sampling
                    # for the rest, the values of n range from
                    # 0 to ks[slots_remaining, n_remaining + 1]
                    ns = 0:ks[slots_remaining, n_remaining+1]
                    
                    # the probability that you select each possible n is
                    # proportional to the # of partitions possible if you
                    # were to choose that n.
                    # divided by the total possible # of partitions across all choices of n.
                    n_parts = vector(length = length(ns), mode = "character")
                    p_ns = vector(length = length(ns), mode = "numeric")
                    for(i in 1:length(n_parts)) {
                        n_parts[i] = ps[slots_remaining-1, n_remaining - (slots_remaining * ns[i]) + 1]
                    }
                    # total partition is
                    total_parts = 0
                    for(i in 1:length(n_parts)) {
                        total_parts = gmp::add.bigz(gmp::as.bigz(n_parts[i]), total_parts)
                    }
                    # probability of each n choice
                    for(i in 1:length(n_parts)){
                        p_ns[i] = as.numeric(gmp::as.bigq(gmp::as.bigz(n_parts[i]), total_parts))
                    }
                    
                    this_increment = sample(ns, size = 1, prob = p_ns)
                    this_gnome[species_slot] = this_gnome[species_slot - 1] + this_increment
                    n_remaining = n_remaining - (slots_remaining*this_increment)
                    slots_remaining = slots_remaining - 1
                    
                } else {
                    this_gnome[species_slot] = this_gnome[species_slot - 1] + n_remaining
                }
            }
            
        }
        
        sets[idx, ] = this_gnome
    }
    return(sets)
}

#' @title Filter and tally frequency of distinct RADs
#'
#' @description Filter a matrix of samples from the feasible set to include only unique vectors, and tally the frequency of those vectors.
#'
#' The frequency should be approximately equal for all vectors, but for small sample numbers this will not necessarily happen.
#' @param sets_matrix matrix of draws from a feasible set (output of `sample_fs`)
#' @return data frame of unique species abundance vectors (rows are unique vectors, columns are species) and a column `set_frequency` of the frequency of each vector in the original set of samples.
#'
#' @export
tally_sets <- function(sets_matrix) {
    
    sets_matrix <- as.data.frame(sets_matrix) %>%
        dplyr::group_by_all() %>%
        dplyr::tally() %>%
        dplyr::ungroup() %>%
        dplyr::rename('set_frequency' = 'n')
    
    return(sets_matrix)
}
