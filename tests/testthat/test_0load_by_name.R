
# context("load_by_name")

test_that("Load data on first request",{
  	expect_that(load_by_name("dataset_5_stages"), shows_message("dataset_5_stages already exists in package environment: FALSE"))
	expect_that(ncol(load_by_name("dataset_5_stages"))==6, is_true())
})

test_that("Do not reload data",{
	  	expect_that(load_by_name("dataset_5_stages"), shows_message("dataset_5_stages already exists in package environment: TRUE"))
})
