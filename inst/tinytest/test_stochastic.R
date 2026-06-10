# ============
# Test vectors
# ===========

# brownian
path <- brownian_vector(dt = 0.01, n = 100)
expect_false(anyNA(path))
expect_length(path, 100)
expect_true(is.numeric(path))

# geometric brownian
path <- geom_brownian_vector(dt = 0.01, n = 100, x0 = 1)
expect_true(all(path > 0))
expect_false(anyNA(path))
expect_length(path, 100)
expect_true(is.numeric(path))

# ornstein-uhlenbeck
path <- ornstein_uhlenbeck_vector(dt = 0.01, n = 100)
expect_false(anyNA(path))
expect_length(path, 100)
expect_true(is.numeric(path))

# =============
# Test matrices
# =============


# brownian
M <- brownian_matrix(dt = 0.01, n_steps = 100, n_traj = 10)
expect_false(anyNA(M))
expect_true(is.matrix(M))
expect_equal(dim(M), c(100, 10))

# geometric brownian
M <- geom_brownian_matrix(dt = 0.01, n_steps = 100, n_traj = 10)
expect_true(all(M>0))
expect_false(anyNA(M))
expect_true(is.matrix(M))
expect_equal(dim(M), c(100, 10))

# ornstein-uhlenbeck
M <- ornstein_uhlenbeck_matrix(dt = 0.01, n_steps = 100, n_traj = 10)
expect_false(anyNA(M))
expect_true(is.matrix(M))
expect_equal(dim(M), c(100, 10))

# =============
# Test end_info
# =============

# brownian
res <- brownian_end_info(dt = 0.01, n_steps = 100, n_traj = 50)
expect_false(anyNA(res$xT))
expect_length(res$xT, 50)
expect_length(res$crossed_upper, 50)
expect_length(res$crossed_lower, 50)
expect_true(is.logical(res$crossed_upper))
expect_true(is.logical(res$crossed_lower))

# Check that thresholds work (Inf cannot be crossed)
res2 <- brownian_end_info(dt = 0.01, n_steps = 100, n_traj = 50,
                          lower_thresh = -Inf, upper_thresh = Inf)
expect_true(all(!res2$crossed_upper))
expect_true(all(!res2$crossed_lower))

# geometric brownian
res <- geom_brownian_end_info(dt = 0.01, n_steps = 100, n_traj = 50)
expect_false(anyNA(res$xT))
expect_length(res$xT, 50)
expect_length(res$crossed_upper, 50)
expect_length(res$crossed_lower, 50)
expect_true(is.logical(res$crossed_upper))
expect_true(is.logical(res$crossed_lower))
expect_true(all(res$T > 0))
# Check that thresholds work (Inf cannot be crossed)
res2 <- brownian_end_info(dt = 0.01, n_steps = 100, n_traj = 50,
                          lower_thresh = -Inf, upper_thresh = Inf)
expect_true(all(!res2$crossed_upper))
expect_true(all(!res2$crossed_lower))

# ornstein-uhlenbeck
res <- ornstein_uhlenbeck_end_info(dt = 0.01, n_steps = 100, n_traj = 50)
expect_false(anyNA(res$xT))
expect_length(res$xT, 50)
expect_length(res$crossed_upper, 50)
expect_length(res$crossed_lower, 50)
expect_true(is.logical(res$crossed_upper))
expect_true(is.logical(res$crossed_lower))

# Check that thresholds work (Inf cannot be crossed)
res2 <- brownian_end_info(dt = 0.01, n_steps = 100, n_traj = 50,
                          lower_thresh = -Inf, upper_thresh = Inf)
expect_true(all(!res2$crossed_upper))
expect_true(all(!res2$crossed_lower))