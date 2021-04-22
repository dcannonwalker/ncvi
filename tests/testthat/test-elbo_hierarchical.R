test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

data("test_hierarchical_poisson")

test_that("elbo_hierarchical agrees with older elbo_ex4", {
  expect_equal(elbo_ex4(data_legacy, init_legacy),
               elbo_hierarchical(data_new, init_new))
})

