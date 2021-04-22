test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

data("test_poisson")

test_that("new and legacy updates agree", {
  nc_update_mvn()
})
