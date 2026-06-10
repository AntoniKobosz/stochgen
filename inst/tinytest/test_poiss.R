lambda <- 2
n <- 50
times <- poiss_jump_times(lambda, n)
# test output length
expect_length(times, n)
# test that times are increasing and positive
expect_true(all(times > 0))
expect_true(all(diff(times) > 0))

# no NAs
expect_false(anyNA(times))