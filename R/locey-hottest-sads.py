def get_hottest_SAD(unique_SADs):
  """ Find the SAD in a random sample with the greatest average commonness
among its ranked abundance states. This SAD is taken to represent the
central tendency of the set, based on the SAD shape. """

if len(unique_SADs) > 500:
    unique_SADs = random.sample(unique_SADs,500)

    N = sum(unique_SADs[0])
    S = len(unique_SADs[0])
    a1 = 0 # SAD mean
    v1 = 0 # SAD variance
    for rad in unique_SADs:
      in_common = []
    ct1 = 0
    for a in rad: # for each rank
      c = 0
    for sad in unique_SADs:
      if a == sad[ct1]:
      c += 1
    in_common.append(ln(RDF(c)))
    ct1 += 1
    a2 = mean(in_common)
    v2 = variance(in_common)
    if a2 > a1:
      a1 = a2
    v1 = v2
    xRAD = rad
    elif a2 == a1:
      if v2 < v1:
      a1 = a2
    v1 = v2
    xRAD = rad
    #percentile_evar = stats.percentileofscore(sample_evar,obs_evar)
    return xRAD
