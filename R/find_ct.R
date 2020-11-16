#' @title Find the central tendency of a collection of RADs.
#'
#' @description Find the central tendency following `get_hottest_SAD` in Locey & White.
#'
#' @param unique_RADs a matrix unique RADs with the same S and N
#' @return vector of the central tendency of `unique_RADs`.
#'
#' @export

find_ct <- function(unique_RADs){

  unique_RADs <- as.matrix(unique_RADs)

  S = ncol(unique_RADs)
  N = sum(unique_RADs[1, ])


 tally_in_common = function(species_index, focal_rad, sample_rads){
    this_tally = length(which(sample_rads[,species_index] == focal_rad[species_index]))
    this_tally = log(this_tally)
    return(this_tally)
  }

  calculate_in_common = function(focal_rad, sample_rads) {
    this_in_common = vapply(1:S, tally_in_common, focal_rad, sample_rads, FUN.VALUE = 1)
    return(this_in_common)
  }

  in_common = list()
  for(i in 1:nrow(unique_RADs)){
    in_common[[i]] <- calculate_in_common(unique_RADs[i, ], unique_RADs)
  }

  in_common_a <- t(apply(unique_RADs, MARGIN = 1, FUN = calculate_in_common, sample_rads = unique_RADs))

  unique_RADs <- as.data.frame(unique_RADs)
  unique_RADs$a = apply(in_common_a, MARGIN = 1, FUN = mean)
  unique_RADs$v = apply(in_common_a, MARGIN = 1, FUN = stats::var)

  xRAD <- dplyr::filter(unique_RADs, .data$a == max(unique_RADs$a))
  if(nrow(xRAD) > 1) {
    xRAD <- dplyr::filter(xRAD, .data$v == min(xRAD$v))
    if(nrow(xRAD > 1)) {
      xRAD <- xRAD[1, ]
      print("multiple ct RADs")
    }
  }

  xRAD <- as.integer(dplyr::select(xRAD, -.data$a, -.data$v))

  return(xRAD)

}
