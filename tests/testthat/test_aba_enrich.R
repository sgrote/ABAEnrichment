
# context("main function aba_enrich")
#aba_enrich=function(genes,dataset="adult",test="hyper",cutoff_quantiles=seq(0.1,0.9,0.1),n_randsets=1000, gene_len=FALSE, circ_chrom=FALSE)


## genes input
##############

genes=rep(0:1,4)
names(genes)=c(324,8312,673,1029,64764,1499,"abc",000)

set.seed(123)
res1=aba_enrich(genes,dataset='adult',cutoff_quantiles=c(0.2,0.7),n_randsets=100)

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

willi = as.integer(runif(length(genes),1,30))
names(willi) = names(genes)

set.seed(123)
res3=aba_enrich(willi,dataset='5_stages',test="wilcoxon",cutoff_quantiles=c(0.4),n_randsets=100)

test_that("corner case 5_stages with one cutoff works",{
	expect_that(nrow(res3$cutoffs), equals(1))
})

test_that("genes that are not in data are not returned in the genes object",{
	expect_that("abc" %in% names(res3$genes), is_false())	
	expect_that("8312" %in% names(res3$genes), is_true())
    expect_that(length(res3$genes), equals(6))
})


# test warnings and more on handling of no coords and data

test_genes = c(paste(rep('FABP', 5), c(2,4:7), sep=''), c('HSPB2' ,'LINC00239', 'TESTI'))
bg_genes = c('NCAPG', 'NGFR', 'NXPH4', 'C21orf59', 'CACNG2', 'AGTR1', 'ANO1', 'BTBD3', 'MTUS1', 'MIRLET7BHG', 'MUELL')
genes = c(rep(1,length(test_genes)), rep(0,length(bg_genes)))
names(genes) = c(test_genes, bg_genes)

test_that("Warning about test and background genes with no expression.",{
    expect_that(aba_enrich(genes,cutoff_quantiles=c(0.7,0.8),n_randsets=5), gives_warning("No expression data for genes: TESTI, MUELL.\n  These genes were not included in the analysis.")) 
    expect_that(aba_enrich(genes,cutoff_quantiles=c(0.7,0.8),n_randsets=5, gene_len=TRUE), gives_warning("No expression data for genes: TESTI, MUELL.\n  These genes were not included in the analysis."))
})    

test_that("One candidate or background gene without expression data stops.",{
    expect_that(aba_enrich(genes[8:length(genes)],cutoff_quantiles=c(0.7,0.8),n_randsets=5), throws_error("No requested test genes in data."))    
    expect_that(aba_enrich(genes[c(1:8,length(genes))],cutoff_quantiles=c(0.7,0.8),n_randsets=5), throws_error("No requested background genes in data."))    
})   

# 'HSPB2' ,'LINC00239', 'MIRLET7BHG' -> no coordinates
test_that("Warning about test and background genes with no coordinates.",{
    expect_that(aba_enrich(genes,cutoff_quantiles=c(0.7,0.8),n_randsets=5,gene_len=TRUE), gives_warning("No coordinates available for genes: HSPB2, LINC00239, MIRLET7BHG.\n  These genes were not included in the analysis."))    
})

test_that("One candidate or background gene without coordinates stops.",{
    expect_that(aba_enrich(genes[6:length(genes)],cutoff_quantiles=c(0.7,0.8),n_randsets=5,gene_len=T), throws_error("No requested test genes in data."))    
    expect_that(aba_enrich(genes[c(1:8,18,19)],cutoff_quantiles=c(0.7,0.8),n_randsets=5,gene_len=T), throws_error("No requested background genes in data."))    
})   

set.seed(123)
res4 = aba_enrich(genes,n_randsets=100,gene_len=TRUE)
test_that("genes that are not in data or coordinates are not returned in the genes object",{
	expect_that(length(res4$genes), equals(length(genes)-5))
	expect_that("NGFR" %in% names(res4$genes), is_true())
	expect_that("TESTI" %in% names(res4$genes), is_false())	
	expect_that("MUELL" %in% names(res4$genes), is_false())	
	expect_that("HSPB2" %in% names(res4$genes), is_false())	
	expect_that("LINC00239" %in% names(res4$genes), is_false())	
	expect_that("MIRLET7BHG" %in% names(res4$genes), is_false())	

})

# full run using all genes as background  
genes = rep(1,13)
names(genes) = c('NCAPG', 'APOL4', 'NGFR', 'NXPH4', 'C21orf59', 'CACNG2', 'AGTR1', 'ANO1',
'BTBD3', 'MTUS1', 'CALB1', 'GYG1', 'PAX2')
#set.seed(123)
#res5 = aba_enrich(genes,dataset = '5_stages',cutoff_quantiles = c(0.5,0.7,0.9), n_randsets = 100)
set.seed(123)
res6 = aba_enrich(genes,dataset = '5_stages',cutoff_quantiles = c(0.5,0.7,0.9), n_randsets = 100, gene_len = TRUE)

#fwers1 = res1[[1]]
#fwers3 = res3[[1]]
#fwers4 = res4[[1]]
##fwers5 = res5[[1]]
#fwers6 = res6[[1]]
#test_that("check some FWERs - single genes.",{
#	expect_that(round(fwers4[fwers4$structure_id=="Allen:4124","mean_FWER"],7), equals(0.9155556))
#	expect_that(round(mean(fwers4$mean_FWER),7), equals(0.9101996))
#	expect_that(round(mean(fwers4$min_FWER),7), equals(0.7193607))	
#	expect_that(round(mean(fwers3$mean_FWER),7), equals(0.9994815))	
#	expect_that(fwers1[fwers1$structure_id=="Allen:9150","mean_FWER"], equals(0.965))
#	expect_that(round(mean(fwers1$mean_FWER),7), equals(0.9876941))
#	expect_that(round(mean(fwers1$min_FWER),7), equals(0.9753881))
##	expect_that(fwers5[1,8], equals("0.7;0.34;0.24"))
#	expect_that(fwers6[1,8], equals("0.43;0.26;0.15"))
#})


## genomic region input
#######################

# normal runs (examples similar to vignette)
genes1 = c(1, rep(0,6))
names(genes1) = c('8:82000000-83000000', '7:1300000-56800000', '7:74900000-148700000',
 '8:7400000-44300000', '8:47600000-146300000', '9:0-39200000', '9:69700000-140200000')
set.seed(123)
res1 = aba_enrich(genes1, n_randsets=100) 

genes2 = c(1,1,rep(0,4))
names(genes2) = c('8:82000000-83000000','7:113600000-124700000','7:1300000-56800000', '7:74900000-148700000',
 '8:7400000-44300000', '8:47600000-146300000')
set.seed(123)
res2 = aba_enrich(genes2, n_randsets=100, circ_chrom=TRUE) 

#fwers1 = res1[[1]]
#fwers2 = res2[[1]]

#test_that("check some FWERs - genomic regions.",{
#	expect_that(fwers1[fwers1$structure_id=="Allen:4735","mean_FWER"], equals(0.5))
#	expect_that(fwers1[fwers1$structure_id=="Allen:4743","FWERs"], equals("0.9;0.7;0.96;0.57;0.5;0.34;0.4;0.02;0.03"))
#	expect_that(round(mean(fwers1$mean_FWER),7), equals(0.5925165))
#	expect_that(round(mean(fwers1$min_FWER),7), equals(0.05977169))
#	expect_that(fwers2[fwers2$structure_id=="Allen:4598","mean_FWER"], equals(0.93))
#	expect_that(fwers2[fwers2$structure_id=="Allen:4743","FWERs"], equals("1;1;1;1;1;0.99;1;0.81;0.87"))
#	expect_that(round(mean(fwers2$mean_FWER),7), equals(0.9829951))
#	expect_that(round(mean(fwers2$min_FWER),7), equals(0.8800761))
#})


# test checking: background < candidate on chrom, in blocks, overlapping input regions
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
#    expect_that(aba_enrich(tight), throws_error("Background regions too small.")) # testthat fails on bioconductor test (actual value: "basic_string::_M_replace_aux"). error is thrown in blocks.cpp by Rcpp::stop. thestthat failure not reproducible locally (with R_DEFAULT_PACKAGES= LC_COLLATE=C LANGUAGE=EN R) 
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
