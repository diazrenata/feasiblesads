library(gmp)
source('R/build_P_table.R')

sample_fs = function(s, n, nsamples) {
ps <- fill_ps(s, n, storeyn = FALSE, storepath = NULL)

# Once you have the ps table.
# You also might as well make a ks table? looking up might be as slow as calculating, but whatever.
ks <- fill_ks(s, n)

sets = list()

while(length(sets) < nsamples) {

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
        total_parts = add.bigz(as.bigz(n_parts[i]), total_parts)
      }
      # probability of each n choice
      for(i in 1:length(n_parts)){
      p_ns[i] = as.numeric(as.bigq(as.bigz(n_parts[i]), total_parts))
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
        total_parts = add.bigz(as.bigz(n_parts[i]), total_parts)
      }
      # probability of each n choice
      for(i in 1:length(n_parts)){
        p_ns[i] = as.numeric(as.bigq(as.bigz(n_parts[i]), total_parts))
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

sets[[length(sets) + 1]] = this_gnome
}
return(sets)
}

system.time(sample_fs(3, 8, 1))

feasible_sets = sample_fs(3, 8, 1000)
