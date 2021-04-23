test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

data("data_hierarchical_legacy")
data("data_hierarchical_new")
data("init_hierarchical_legacy")
data("init_hierarchical_new")

test_that("elbo_hierarchical agrees with older elbo_ex4", {
  expect_equal(elbo_legacy(data_hierarchical_legacy, init_hierarchical_legacy),
               elbo_hierarchical(data_hierarchical_new, init_hierarchical_new))
})
elbo_hierarchical(data_hierarchical_new,
                  init_hierarchical_new)
elbo_legacy(data_hierarchical_legacy, init_hierarchical_legacy)
elbo_hierarchical
elbo_extra
