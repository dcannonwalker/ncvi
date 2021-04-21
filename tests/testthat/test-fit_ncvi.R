test_that("args get passed to update_pars() correctly", {
  data <- list()
  init <- list()
  elbo <- function(data, pars) 1
  update_pars <- function(data, pars, args) {
    print(args)
    pars
  }
  fit_ncvi(data, pars)
})
