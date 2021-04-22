test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})


test_that("different mean differentials functions agree", {
  expect_equal(
    differential_EB(data1, pars1),
    d_mvn_mean(data2, pars2, d_mvn_mean_poisson))
})

test_that("different cov differentials functions agree", {
  expect_equal(
    d_mvn_cov(data2, pars2, d_mvn_cov_poisson),
    differential_SigmaB(data1, pars1)
  )
})


