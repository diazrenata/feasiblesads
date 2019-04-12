context("Sample feasible sets")

test_that("Generate some feasible sets", {
    set.seed(42)
    expect_error(output <- sample_fs(3, 8, 20), NA)
    expect_equal(dim(output), c(20, 3))
    expect_true(all(rowSums(output) == 8))
    expect_known_hash(output, "0de117829e")
})