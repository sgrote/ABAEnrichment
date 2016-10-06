
# context("get sampled substructures")

test_that("sampled substructures from developmental ontology are returned",{
  	expect_that(get_sampled_substructures("Allen:13324"), equals("Allen:10252"))
  	expect_that(get_sampled_substructures("Allen:10208"), equals(c("Allen:10209","Allen:10225")))
})

test_that("sampled substructures from adult ontology are returned",{
  	expect_that(get_sampled_substructures("Allen:4011"), equals(c("Allen:4012", "Allen:4013", "Allen:4014", "Allen:4015")))
	expect_that(get_sampled_substructures("Allen:4020"), equals("Allen:4020"))
  	expect_that(get_sampled_substructures("Allen:4016"), equals(c("Allen:4017", "Allen:4018", "Allen:4019", "Allen:4020")))
})

test_that("multiple structures are not allowed",{
	expect_that(get_sampled_substructures(c("Allen:10194","Allen:10153")), throws_error("Please use a single structure id."))
})

test_that("error message for invalid structure is returned",{
	expect_that(get_sampled_substructures("Allen:123"), throws_error("No data for structure_id Allen:123."))
})



## new feature version 1.3.5 - allow also structure IDs without "Allen:"-prefix. copied from above  

test_that("sampled substructures from developmental ontology are returned",{
  	expect_that(get_sampled_substructures(13324), equals("10252"))
  	expect_that(get_sampled_substructures("10208"), equals(c("10209","10225")))
})

test_that("sampled substructures from adult ontology are returned",{
  	expect_that(get_sampled_substructures(4011), equals(c("4012", "4013", "4014", "4015")))
	expect_that(get_sampled_substructures("4020"), equals("4020"))
  	expect_that(get_sampled_substructures(4016), equals(c("4017", "4018", "4019", "4020")))
})

test_that("multiple structures are not allowed",{
	expect_that(get_sampled_substructures(c(10194,10153)), throws_error("Please use a single structure id."))
})

test_that("error message for invalid structure is returned",{
	expect_that(get_sampled_substructures(123), throws_error("No data for structure_id 123."))
})
