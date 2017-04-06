

# context("get IDs")

test_that("some IDs given names",{
	
	expect_that(nrow(get_id("")), equals(725))
  	expect_that(get_id("accumbens")$structure_id, equals(c("Allen:4290","Allen:4291","Allen:4292")))
  	expect_that(get_id("telencephalon")$structure_id, equals(c("Allen:10158","Allen:4007")))
  	expect_that(get_id("telencephalon")$ontology, equals(c("developmental","adult")))
  	
})

test_that("error message for invalid structure name",{
	expect_that(get_id("quatsch"), throws_error("No matches found"))
	expect_that(get_id(c("accumbens", "telencephalon")), throws_error("Please use a single brain region name."))
})

