

# context("plot expression")

# additional checks concerning dataset and genes are performed by get_expression() inside plot_expression

test_that("normal input results in heatmap-object (list of 9)",{
	expect_that(class(plot_expression(structure_ids=c('Allen:4010'),gene_ids=c(324,8312,673,1029,64764,1499),dataset='adult')), equals("list"))
	expect_that(length(plot_expression(structure_ids=c('Allen:4010'),gene_ids=c(324,8312,673,1029,64764,1499),dataset='adult')), equals(9))
})

test_that("edge cases are handled",{
	# only two genes and two structures
	expect_that(class(plot_expression(structure_ids=c('Allen:10657','Allen:10398'),gene_ids=c(324,8312),dataset="5_stages")),equals("list"))
	# age category as a string
	expect_that(class(plot_expression(structure_ids=c('Allen:13322'),gene_ids=c(324,8312),dataset="5_stages", age_category="3")),equals("list"))
	# only one structure is returned from get_expression (cannot be checked before because is might have multiple sampled substructures)
	expect_that(plot_expression(structure_ids=c('Allen:10655'),gene_ids=c(324,8312),dataset="5_stages"),throws_error("At least two brain structures are needed for plotting."))	
	# only one gene is returned from get_expression
	expect_that(plot_expression(structure_ids=c('Allen:10655','Allen:10398'),gene_ids=c(000,324,"abc"),dataset="5_stages"),throws_error("At least two genes are needed for plotting."))
})

test_that("arguments are checked",{
	# define neither genes nor dataset without having run aba_enrich previously
	expect_that(plot_expression(structure_ids=c('Allen:13322')), throws_error())
	# one gene
	expect_that(plot_expression(structure_ids=c('Allen:13322'),gene_ids=c(8312),dataset="5_stages"), throws_error("Please enter at least two genes."))	
	# wrong age_category
	expect_that(plot_expression(structure_ids=c('Allen:13322'),gene_ids=c(324,8312),dataset="5_stages",age_category=10), throws_error("Please specify an age_category between 1 and 5."))	
	# multiple age categories
	expect_that(plot_expression(structure_ids=c('Allen:13322'),gene_ids=c(324,8312),dataset="5_stages",age_category=1:2), throws_error("Please specify an age_category between 1 and 5."))
})

