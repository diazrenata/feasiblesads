#' @title Draw samples from the feasible set given S and N
#' @description Draw samples from the feasible set of species abundance distributions with `s` species and `n` total individuals
#' @param s the number of species
#' @param n the number of individuals
#' @param nsamples how many samples to draw
#' @param storeyn TRUE/FALSE whether to store the P table.
#' @param storepath path to store the P table.
#' @param p_table optionally pass P table as an argument
#' @return matrix `nsamples` of ranked abundance distributions. Rows are samples and columns are species, sorted from least to most abundant.
#' @export

sample_fs = function(s, n, nsamples, storeyn = FALSE,
                     storepath = NULL, p_table = NULL)
{
  if(is.null(p_table)) {
    ps <- fill_ps(s, n, storeyn = storeyn, storepath = storepath)
  } else {
    ps <- p_table
  }
  # Once you have the ps table you also might as well make a ks table? looking up might be as slow as calculating, but whatever.
  ks <- fill_ks(s, n)

  sets <- matrix(NA, nrow = nsamples, ncol = s)
  for (idx in seq_len(nsamples))
  {
    sets[idx, ] <- send_gnome(s, n, ps, ks)
  }
  return(sets)
}

#' @title Draw a single sample from the feasible set
#'
#' @description Draw a single sample from the feasible set of species abundance
#'   distributions with `s` species and `n` total individuals, uniformly
#'
#' @inheritParams sample_fs
#' @param ps the ps table for number of partitions
#' @param ks the ks table of max partition sizes
#' @return a vector of length `s` whose elements sum up to `n`, in non-
#'   decreasing order
#' @export
send_gnome <- function(s, n, ps, ks)
{
  this_gnome = vector(length = s, mode = 'integer')

  if (s == 1) # if s = 1 there is only one possible RAD.
  {
    this_gnome[1] <- n
    return(this_gnome)
  }

  slots_remaining = s
  n_remaining = n
  current_size <- 0
  for (species_slot in 1:s)
  {
    # final species slot
    if (species_slot == s) {
      this_gnome[species_slot] = current_size + n_remaining
      break
    }

    if(species_slot == 1) {
      # first (lowest abundance) value is a boundary case
      # for the first one, the values of n range from 1 to ks[s, n+1]
      ns = 1:ks[slots_remaining, n_remaining + 1]
    } else {
      # recursive probabilistic sampling
      # for the rest, the values of n range from
      # 0 to ks[slots_remaining, n_remaining + 1]
      ns = 0:ks[slots_remaining, n_remaining+1]
    }

    # the probability that you select each possible n is
    # proportional to the # of partitions possible if you
    # were to choose that n.
    # divided by the total possible # of partitions across all choices of n.
    n_parts = ps[slots_remaining - 1,
                 n_remaining - (slots_remaining * ns) + 1]

    # total partition size
    total_parts = gmp::sum.bigz(gmp::as.bigz(n_parts))

    # probability of each n choice
    p_ns = as.numeric(gmp::as.bigq(gmp::as.bigz(n_parts), total_parts))

    # randomly select the next increment:
    #   update the current_size,
    #   append it to figure out the current species_slot
    #   compute the abundance and slots remaining
    #   save the current increment size
    this_increment = sample(ns, size = 1, prob = p_ns)
    current_size <- current_size + this_increment
    this_gnome[species_slot] = current_size
    n_remaining = n_remaining - (slots_remaining * this_increment)
    slots_remaining = slots_remaining - 1
  }

  return(this_gnome)
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
