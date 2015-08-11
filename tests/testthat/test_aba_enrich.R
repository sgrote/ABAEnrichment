
# context("main function aba_enrich")
#aba_enrich=function(genes,dataset="adult",test="hyper",cutoff_quantiles=seq(0.1,0.9,0.1),n_randsets=1000)

genes=rep(0:1,4)
names(genes)=c(324,8312,673,1029,64764,1499,"abc",000)

res1=aba_enrich(genes,dataset='adult',cutoff_quantiles=c(0.2,0.7),n_randsets=5)

test_that("normal input results in list",{
	expect_that(class(res1), equals("list"))
	expect_that(nrow(res1$cutoffs), equals(2))
})

res2=aba_enrich(genes,dataset='5_stages',cutoff_quantiles=c(0.7),n_randsets=5)

test_that("corner case 5_stages with one cutoff works",{
	expect_that(nrow(res2$cutoffs), equals(1))
})

test_that("genes that are not in data are not returned in the genes object",{
	expect_that("abc" %in% names(res2$genes), is_false())	
	expect_that("8312" %in% names(res2$genes), is_true())
    expect_that(length(res2$genes), equals(6))
})
