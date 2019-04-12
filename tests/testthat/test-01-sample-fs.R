context("Sample feasible sets")

test_that("Check that all feasible sets with s = 3, n = 8 are made", {
    set.seed(42)
    expect_error(output <- sample_fs(3, 8, 20), NA)
    expect_equal(dim(output), c(20, 3))
    expect_true(all(rowSums(output) == 8))
    expect_known_hash(output, "0de117829e")
    
    unique_sads <- dplyr::distinct(as.data.frame(output))
    expect_equal(NROW(unique_sads), 5)
})

test_that("Check that all feasible sets with s = 3, n = 8 are made", {
    set.seed(42)
    sad_freq <- sample_fs(3, 8, 10000) %>%
        as.data.frame() %>%
        dplyr::group_by_all() %>%
        dplyr::tally()
    expect_true(all(sad_freq$n > 1900))
    expect_true(all(sad_freq$n < 2100))
})

test_that("Generate some feasible sets", {
    set.seed(42)
    expect_error(output <- sample_fs(4, 20, 1000), NA)
    expect_equal(dim(output), c(1000, 4))
    expect_true(all(rowSums(output) == 20))
    expect_known_hash(output, "c9253177ed")
    
    unique_sads <- dplyr::distinct(as.data.frame(output))
    expect_equal(NROW(unique_sads), 64)
})
    