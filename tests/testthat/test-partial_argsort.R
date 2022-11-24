N_MAX <- 10000000L
vector_length <- 100000L
n_elements <- 1000L
test_that("partial_argsort works", {
    x <- sample(x = 1:N_MAX, size = vector_length, replace = FALSE)
    res <- partial_argsort(x, n_elements)
    expected_res <- order(x, decreasing = TRUE)[1:n_elements]
    expect_equal(res, expected_res)
})

test_that("partial_argsort provides expected outout when requesting too many elements", {
  x <- sample(x = 1:N_MAX, size = vector_length, replace = FALSE)
  res <- partial_argsort(x, 2L*vector_length)
  expected_res <- order(x, decreasing = TRUE)[1:(2L*vector_length)]
  expect_equal(res, expected_res)
})
