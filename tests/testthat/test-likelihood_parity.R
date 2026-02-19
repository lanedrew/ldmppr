test_that("full_sc_lhood and full_sc_lhood_fast agree", {
  set.seed(42)

  mk_data <- function(n = 40L, xy = 25) {
    t <- sort(stats::runif(n, 0.02, 0.98))
    x <- stats::runif(n, 0, xy)
    y <- stats::runif(n, 0, xy)
    cbind(time = t, x = x, y = y)
  }

  run_case <- function(data, params, grid) {
    b <- c(1, 25, 25)
    old <- ldmppr:::full_sc_lhood(grid$x, grid$y, grid$t, data[, 1], data, params, b)
    fast <- ldmppr:::full_sc_lhood_fast(grid$x, grid$y, grid$t, data[, 1], data, params, b)
    c(old = old, fast = fast)
  }

  g1 <- list(
    x = seq(0, 25, length.out = 8),
    y = seq(0, 25, length.out = 8),
    t = seq(0.02, 0.98, length.out = 12)
  )
  g2 <- list(
    x = seq(0, 25, length.out = 12),
    y = seq(0, 25, length.out = 12),
    t = sort(c(seq(0.02, 0.98, length.out = 10), 0.5, 0.5, 0.7))
  )

  d1 <- mk_data(35)
  d2 <- mk_data(60)
  d3 <- d1
  d3[2, 1] <- d3[1, 1]
  d3[3, 1] <- d3[1, 1]
  d3 <- d3[order(d3[, 1]), , drop = FALSE]

  p1 <- c(0.1, 6, 0.02, 2.5, 2, 0.6, 1.2, 0.05)
  p2 <- c(-0.5, 3, 0.01, 3.5, 4, 0, 1.5, 0.1)
  p3 <- c(0.2, 10, 0.03, 1.8, 1.2, 0.8, 3.0, 0.0)

  r1 <- run_case(d1, p1, g1)
  r2 <- run_case(d2, p1, g2)
  r3 <- run_case(d2, p2, g1)
  r4 <- run_case(d1, p3, g2)
  r5 <- run_case(d3, p1, g1)

  expect_equal(unname(r1["old"]), unname(r1["fast"]), tolerance = 1e-10)
  expect_equal(unname(r2["old"]), unname(r2["fast"]), tolerance = 1e-10)
  expect_equal(unname(r3["old"]), unname(r3["fast"]), tolerance = 1e-10)
  expect_equal(unname(r4["old"]), unname(r4["fast"]), tolerance = 1e-10)
  expect_equal(unname(r5["old"]), unname(r5["fast"]), tolerance = 1e-10)
})


test_that("full_sc_lhood and full_sc_lhood_fast match duplicate-location edge case", {
  set.seed(7)

  d <- cbind(
    time = sort(stats::runif(25, 0.05, 0.95)),
    x = stats::runif(25, 0, 25),
    y = stats::runif(25, 0, 25)
  )
  d[6, 2:3] <- d[2, 2:3]

  p <- c(0.2, 4, 0.02, 3, 2, 0.7, 2.0, 0.1)
  g <- list(
    x = seq(0, 25, length.out = 10),
    y = seq(0, 25, length.out = 10),
    t = seq(0.05, 0.95, length.out = 15)
  )
  b <- c(1, 25, 25)

  old <- ldmppr:::full_sc_lhood(g$x, g$y, g$t, d[, 1], d, p, b)
  fast <- ldmppr:::full_sc_lhood_fast(g$x, g$y, g$t, d[, 1], d, p, b)

  expect_true(is.infinite(old) && old < 0)
  expect_true(is.infinite(fast) && fast < 0)
})
