
# context("main function aba_enrich")
#aba_enrich=function(genes,dataset="adult",test="hyper",cutoff_quantiles=seq(0.1,0.9,0.1),n_randsets=1000, gene_len=FALSE, circ_chrom=FALSE)

set.seed(123)

## genes input
genes=rep(0:1,4)
names(genes)=c(324,8312,673,1029,64764,1499,"abc",000)

res1=aba_enrich(genes,dataset='adult',cutoff_quantiles=c(0.2,0.7),n_randsets=5)

test_that("normal input results in list",{
	expect_that(class(res1), equals("list"))
	expect_that(nrow(res1$cutoffs), equals(2))
	expect_that(length(unlist(strsplit(res1[[1]][1,8],";"))), equals(2))
})

res2=aba_enrich(genes,dataset='adult',cutoff_quantiles=c(0.999,0.2),n_randsets=5)

test_that("cutoffs get sorted, output for successful cutoffs is returned",{
	expect_that(class(res2), equals("list"))
	expect_that(length(unlist(strsplit(res2[[1]][1,8],";"))), equals(1))
})

res3=aba_enrich(genes,dataset='5_stages',cutoff_quantiles=c(0.7),n_randsets=5)

test_that("corner case 5_stages with one cutoff works",{
	expect_that(nrow(res3$cutoffs), equals(1))
})

test_that("genes that are not in data are not returned in the genes object",{
	expect_that("abc" %in% names(res3$genes), is_false())	
	expect_that("8312" %in% names(res3$genes), is_true())
    expect_that(length(res3$genes), equals(6))
})

## genomic region input

# test background < candidate on chrom, in blocks, overlapping input regions
too_small = c(rep(1,4),rep(0,5)) 
names(too_small) = c("X:0-2","2:0-3","2:5-10","13:0-20",  "2:5-10","13:0-10","X:0-1","4:0-10","2:10-12") 
overlap = too_small
names(overlap) = c("X:0-1","2:0-3","2:5-10","13:0-20",  "2:5-20","13:0-30","X:0-1","4:0-10","2:10-12")
overlap2 = too_small
names(overlap2)[2] = "2:0-8"
tight = c(rep(1,4),0)
names(tight) = c("1:104000000-114900000", "3:76500000-90500000", "7:113600000-124700000", "8:54500000-65400000", "5:0-47000000")
no_bg = too_small[1:4]
reverse = too_small
names(reverse)[3] = "2:15-10"

test_that("input_regions are checked - blocks",{
    expect_that(aba_enrich(too_small), throws_error("Candidate regions bigger than any background region:\n  13:0-20")) 
    expect_that(aba_enrich(overlap), throws_error("Background regions overlap: 2:5-20, 2:10-12"))
    expect_that(aba_enrich(overlap2), throws_error("Candidate regions overlap: 2:0-8, 2:5-10"))
    expect_that(aba_enrich(tight), throws_error("Background regions too small."))
    expect_that(aba_enrich(no_bg), throws_error("All values of the 'genes' input are 1. Using chromosomal regions as input requires defining background regions with 0."))
    expect_that(aba_enrich(reverse), throws_error("Invalid regions: 2:15-10.\n  In 'chr:start-stop' start < stop is required."))
})    

test_that("input_regions are checked - circ_chrom",{
    expect_that(aba_enrich(too_small, circ_chrom=TRUE), throws_error("Sum of candidate regions is bigger than sum of background regions on chromosomes: 2, 13, X"))
    expect_that(aba_enrich(overlap, circ_chrom=TRUE), throws_error("Background regions overlap: 2:5-20, 2:10-12"))
    expect_that(aba_enrich(overlap2, circ_chrom=TRUE), throws_error("Candidate regions overlap: 2:0-8, 2:5-10"))      
    expect_that(aba_enrich(tight, circ_chrom=TRUE), throws_error("No background region for chromosomes: 1, 3, 7, 8.\n  With circ_chrom=TRUE only background regions on the same chromosome as a candidate region are used."))
    expect_that(aba_enrich(no_bg, circ_chrom=TRUE), throws_error("All values of the 'genes' input are 1. Using chromosomal regions as input requires defining background regions with 0."))    
    expect_that(aba_enrich(reverse, circ_chrom=TRUE), throws_error("Invalid regions: 2:15-10.\n  In 'chr:start-stop' start < stop is required."))
})    
