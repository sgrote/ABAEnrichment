
# context("get expression")

test_that("normal input results in dataframe (adult, dev_effect) or list( 5_stages)",{
#	expect_that(dim(get_expression(structure_ids=c('Allen:4010'),gene_ids=c(324,8312,673,1029,64764,1499),dataset='adult')), equals(c(8,6)))
	expect_that(dim(get_expression(structure_ids=c('Allen:13322'),gene_ids=c('ENSG00000168036', 'ENSG00000157764', 'ENSG00000163041'),dataset='dev_effect')),equals(c(2,3)))
	expect_that(class(get_expression(structure_ids=c('Allen:13322'),gene_ids=c('NCAPG', 'CACNG2', 'NGFR'),dataset='5_stages')),equals(c("list")))
	expect_that(length(get_expression(structure_ids=c('Allen:13322'),gene_ids=c('NCAPG', 'CACNG2', 'NGFR'),dataset='5_stages')),equals(5))
	# only one gene and one structure
#	expect_that(dim(get_expression(structure_ids=c('Allen:4015'),gene_ids=c(324),dataset='adult')), equals(c(1,1)))
	expect_that(class(get_expression(structure_ids=c('Allen:10657'),gene_ids=c(324),dataset="5_stages")),equals("list"))	
})


test_that("arguments are checked",{
	# define neither genes nor dataset without having run aba_enrich previously
	expect_that(get_expression(structure_ids=c('Allen:13322')), throws_error())
	# define only genes withthout having run aba_enrich previously
	expect_that(get_expression(structure_ids=c('Allen:13322'),gene_ids=c('NCAPG', 'CACNG2', 'NGFR')), throws_error())
	# define neither genes nor dataset without having run aba_enrich previously
	expect_that(get_expression(structure_ids=c('Allen:13322'),dataset="5_stages"), throws_error())
	# dataset is checked
	expect_that(get_expression(structure_ids=c('Allen:10657'),gene_ids=c(324),dataset="stages"), throws_error("Not a valid dataset. Please use 'adult', '5_stages' or 'dev_effect'."))	
	# structures without data lead lead to error (from get_sampled_substructures())
	expect_that(get_expression(structure_ids=c('Allen:123'),gene_ids=c(324),dataset="5_stages"),throws_error("No data for structure_id Allen:123."))
	# genes without data get listed in an error message
	expect_that(get_expression(structure_ids=c('Allen:10657'),gene_ids=c(324,"abc",12345),dataset="5_stages"),gives_warning("No expression data for genes: abc, 12345."))
})

