#' @title Find the central tendency of a collection of RADs.
#'
#' @description Find the central tendency following `get_hottest_SAD` in Locey & White.
#'
#' @param unique_RADs a list of unique RADs with the same S and N
#' @return vector of the central tendency of `unique_RADs`.
#'
#' @export

find_ct <- function(unique_RADs){

  S = length(unique_RADs[[1]])
  N = sum(unique_RADs[[1]])

  a1 = 0
  v1 = 0

  for (i in 1:length(unique_RADs)) {

    focal_rad = unique_RADs[[i]]

    in_common = vector()

    for (j in 1:S) {
      c = 0

      for (k in 1:length(unique_RADs)) {
        if(focal_rad[j] == unique_RADs[[k]][j]) {
          c = c + 1
        }
      }

      in_common = append(in_common, log(c))
    }

    a2 = mean(in_common)
    v2 = var(in_common)

    if(a2 > a1) {
      a1 = a2
      v1 = v2
      xRAD = focal_rad
    } else if (a2 == a1) {
      if(v2 < v1) {
        a1 = a2
        v1 = v2
        xRAD = focal_rad
      }
    }
    }

  return(xRAD)

}
