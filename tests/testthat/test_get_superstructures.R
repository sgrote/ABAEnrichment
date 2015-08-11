
  
# context("get superstructures")

test_that("superstructures from developmental ontology are returned",{
  	expect_that(get_superstructures("Allen:10654"), equals(c("Allen:10153", "Allen:10154", "Allen:10155", "Allen:10653", "Allen:10654")))
	expect_that(length(get_superstructures("Allen:10398")), equals(10))
})

test_that("superstructures from adult ontology are returned",{
  	expect_that(get_superstructures("Allen:4009"), equals(c("Allen:4005", "Allen:4006", "Allen:4007", "Allen:4008", "Allen:4009")))
	expect_that(length(get_superstructures("Allen:4020")), equals(8))
})

test_that("superstructures include structure itself",{
	expect_that("Allen:10194" %in% (get_superstructures("Allen:10194")), is_true())
	expect_that("Allen:4182" %in% (get_superstructures("Allen:4182")), is_true())
})

test_that("multiple structures are not allowed",{
	expect_that(get_superstructures(c("Allen:10194","Allen:10153")), throws_error("Please use a single structure id."))
})

test_that("error message for invalid structure is returned",{
	expect_that(get_superstructures("Allen:123"), throws_error("No data for structure_id Allen:123."))
})
