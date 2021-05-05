N_MAX <- 1e+9L
vector_length <- 1e+7L
n_elements <- 1000L
test_that("partial_argsort works",
          {
            x <- sample(x = 1:N_MAX,size = vector_length,replace = FALSE)
            res <- partial_argsort(x,n_elements)
            expected_res <- order(x,decreasing = TRUE)[1:n_elements]
            expect_equal(res,expected_res)
          }  
          )
