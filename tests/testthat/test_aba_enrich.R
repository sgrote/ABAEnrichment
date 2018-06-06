
# context("main function aba_enrich")
#aba_enrich=function(genes,dataset="adult",test="hyper",cutoff_quantiles=seq(0.1,0.9,0.1),n_randsets=1000, gene_len=FALSE, circ_chrom=FALSE, ref_genome="grch37")

## TODO:
	## maybe check that some brain regions are in the top-ten (actual FWER migth be different)
	## check for messages and hyper has a) no candidate and b) no background at a given threshold

## genes input
##############

# NEW: also check that 8312 as candi and bg genes does not throw "Multiple assignment" error
genes=c(rep(0:1,4),0)
names(genes)=c(324,8312,673,1029,64764,1499,"abc",000,8312)

set.seed(123)
res1=aba_enrich(genes,dataset="5_stages",cutoff_quantiles=c(0.2,0.7),n_randsets=5,silent=TRUE)

test_that("normal input results in list",{
	expect_that(class(res1), equals("list"))
	expect_that(nrow(res1$cutoffs), equals(2))
	expect_that(length(unlist(strsplit(res1[[1]][1,8],";"))), equals(2))
})

res2=aba_enrich(genes,dataset="5_stages",cutoff_quantiles=c(0.999,0.2),n_randsets=5,silent=TRUE)

test_that("cutoffs get sorted, output for successful cutoffs is returned",{
	expect_that(class(res2), equals("list"))
	expect_that(length(unlist(strsplit(res2[[1]][1,8],";"))), equals(1))
	expect_that(rownames(res2[[3]]), equals(c("20%","99.9%")))
})

willi = as.integer(runif(8,1,30))
names(willi) = c(324,8312,673,1029,64764,1499,"abc",000)

set.seed(123)
res3=aba_enrich(willi,dataset="5_stages",test="wilcoxon",cutoff_quantiles=c(0.4),n_randsets=5,silent=TRUE)

test_that("corner case 5_stages with one cutoff works",{
	expect_that(nrow(res3$cutoffs), equals(1))
})

test_that("genes that are not in data are not returned in the genes object",{
	expect_that("abc" %in% res3$genes[,1], is_false())	
	expect_that("8312" %in% res3$genes[,1], is_true())
    expect_that(nrow(res3$genes), equals(6))
})

# binomial test
set.seed(123)
high_genes= c("RFFL", "NTS", "LIPE", "GALNT6", "GSN", "BTBD16", "SOX8", "CERS2")
low_genes= c("GDA", "ENC1", "FOXG1", "EGR4", "NPTX1", "VIPR1", "DOC2A", "OASL", "FRY", "NAV3")
count_A = c(sample(20:30, length(high_genes)), sample(5:15, length(low_genes)))
count_B = c(sample(5:15, length(high_genes)), sample(20:30, length(low_genes)))
genes = data.frame(gene=c(high_genes, low_genes), count_A, count_B)
binom = aba_enrich(genes, n_randsets=5, cutoff_quantiles=c(0.2,0.9), test='binomial', silent=TRUE)

test_that("binomial also works",{
	expect_that(class(binom), equals("list"))
	expect_that(length(binom), equals(3))
	expect_that(colnames(binom$genes), equals(c("gene", "count_A", "count_B")))
	expect_that(nrow(binom$cutoffs), equals(2))
})

# contingency test
set.seed(123)
high_substi_genes = c("RFFL", "NTS", "LIPE", "GALNT6", "GSN", "BTBD16", "SOX8", "CERS2")
low_substi_genes = c("GDA", "ENC1", "FOXG1", "EGR4", "NPTX1", "VIPR1", "DOC2A", "OASL", "FRY", "NAV3")
subs_syn = sample(45:55, length(c(high_substi_genes, low_substi_genes)), replace=T)
subs_non_syn = c(sample(15:25, length(high_substi_genes), replace=T), sample(0:10, length(low_substi_genes)))
vari_syn = sample(25:35, length(c(high_substi_genes, low_substi_genes)), replace=T)
vari_non_syn = c(sample(0:10, length(high_substi_genes), replace=T), sample(10:20, length(low_substi_genes)))
genes = data.frame(genes=c(high_substi_genes, low_substi_genes), vari_syn, vari_non_syn, subs_syn, subs_non_syn)
conti = aba_enrich(genes, test='contingency', n_randset=5, cutoff_quantiles=c(0.7,0.8,0.9), silent=TRUE)

test_that("contingency also works",{
	expect_that(class(conti), equals("list"))
	expect_that(length(conti), equals(3))
	expect_that(colnames(conti$genes), equals(c("genes", "vari_syn","vari_non_syn","subs_syn","subs_non_syn")))
	expect_that(nrow(conti$cutoffs), equals(3))
})

# test warnings and more on handling of no coords and data

test_genes = c(paste(rep("FABP", 5), c(2,4:7), sep=""), c("HSPB2" ,"LINC00239", "TESTI"))
bg_genes = c("NCAPG", "NGFR", "NXPH4", "C21orf59", "CACNG2", "AGTR1", "ANO1", "BTBD3", "MTUS1", "MIRLET7BHG", "MUELL")
genes = c(rep(1,length(test_genes)), rep(0,length(bg_genes)))
names(genes) = c(test_genes, bg_genes)

test_that("Warning about test and background genes with no expression.",{
    expect_that(aba_enrich(genes,cutoff_quantiles=c(0.7,0.8),n_randsets=5,silent=TRUE), gives_warning("No expression data for genes: TESTI, MUELL.\n  These genes were not included in the analysis.")) 
    expect_that(aba_enrich(genes,cutoff_quantiles=c(0.7,0.8),n_randsets=5, gene_len=TRUE,silent=TRUE), gives_warning("No expression data for genes: TESTI, MUELL.\n  These genes were not included in the analysis."))
})    

test_that("One candidate or background gene without expression data stops.",{
    expect_that(aba_enrich(genes[8:length(genes)],cutoff_quantiles=c(0.7,0.8),n_randsets=5), throws_error("No requested test genes in data."))    
    expect_that(aba_enrich(genes[c(1:8,length(genes))],cutoff_quantiles=c(0.7,0.8),n_randsets=5), throws_error("No requested background genes in data."))    
})

## TODO: here check warning error for certain cutoff (a: first cutoff fails, no candidate; b: another cutoff fails, background)

# "LINC00239", "MIRLET7BHG" -> no coordinates
test_that("Warning about test and background genes with no coordinates.",{
    expect_that(aba_enrich(genes,cutoff_quantiles=c(0.7,0.8),n_randsets=5,gene_len=TRUE,silent=TRUE), gives_warning("No coordinates available for genes: LINC00239, MIRLET7BHG.\n  These genes were not included in the analysis."))    
})

test_that("One candidate or background gene without coordinates stops.",{
    expect_that(aba_enrich(genes[7:length(genes)],cutoff_quantiles=c(0.7,0.8),n_randsets=5,gene_len=TRUE), throws_error("No requested test genes in data."))    
    expect_that(aba_enrich(genes[c(1:8,18,19)],cutoff_quantiles=c(0.7,0.8),n_randsets=5,gene_len=TRUE), throws_error("No requested background genes in data."))    
})   

set.seed(123)
res4 = aba_enrich(genes,n_randsets=5,gene_len=TRUE,silent=TRUE)
test_that("genes that are not in data or coordinates are not returned in the genes object",{
	expect_that(nrow(res4$genes), equals(length(genes)-4))
	expect_that("NGFR" %in% res4$genes[,1], is_true())
	expect_that("HSPB2" %in% res4$genes[,1], is_true()) # not in hgnc-symbol
	expect_that("TESTI" %in% res4$genes[,1], is_false())	
	expect_that("MUELL" %in% res4$genes[,1], is_false())	
	expect_that("LINC00239" %in% res4$genes[,1], is_false())	
	expect_that("MIRLET7BHG" %in% res4$genes[,1], is_false())	

})

# grch38
set.seed(123)
res4b = aba_enrich(genes,n_randsets=5,gene_len=TRUE,ref_genome="grch38",silent=TRUE)
test_that("genes that are not in data or coordinates are not returned in the genes object",{
	expect_that(nrow(res4b$genes), equals(length(genes)-2))
	expect_that("NGFR" %in% res4b$genes[,1], is_true())
	expect_that("HSPB2" %in% res4b$genes[,1], is_true()) # not in hgnc-symbol
	expect_that("TESTI" %in% res4b$genes[,1], is_false())	
	expect_that("MUELL" %in% res4b$genes[,1], is_false())	
	expect_that("LINC00239" %in% res4b$genes[,1], is_true()) # not in hg19	
	expect_that("MIRLET7BHG" %in% res4b$genes[,1], is_true())# not in hg19
})

## test unused arguments
test_that("Warning about unused arguments.",{
    expect_that(aba_enrich(genes,test="wilcoxon",cutoff_quantiles=c(0.7,0.8),n_randsets=5,gene_len=TRUE), throws_error("Argument 'gene_len = TRUE' can only be used with 'test = 'hyper''."))    
    expect_that(aba_enrich(genes,cutoff_quantiles=c(0.7,0.8),n_randsets=5,ref_genome="grch38",silent=TRUE), gives_warning("Unused argument: 'ref_genome = grch38'."))    
})

genes_willi = 1:10
names(genes_willi) = c("NCAPG", "APOL4", "NGFR", "NGFR", "NXPH4", "C21orf59", "AGTR1", "CACNG2", "AGTR1", "ANO1")
test_that("multiple assignment of scores to genes throws error",{
    expect_that(aba_enrich(genes_willi, test="wilcoxon"), throws_error("Genes with multiple assignment in input: AGTR1, NGFR"))
})


# test checking: background < candidate on chrom, in blocks, overlapping input regions
too_small = c(rep(1,4),rep(0,5)) 
names(too_small) = c("X:0-2","2:0-3","2:5-10","13:0-20",  "2:5-15","13:0-10","X:0-3","4:0-10","2:10-12") 
overlap = too_small
names(overlap) = c("X:0-1","2:0-3","2:5-10","13:0-20",  "2:5-20","13:0-30","X:0-3","4:0-10","2:10-12")
overlap2 = too_small
names(overlap2)[2] = "2:0-8"
tight = c(rep(1,4),0)
names(tight) = c("1:104000000-114900000", "3:76500000-90500000", "7:113600000-124700000", "8:54500000-65400000", "5:0-4700000")
no_bg = too_small[1:4]
reverse = too_small
names(reverse)[3] = "2:15-10"
no_can_genes = c(1, rep(0,6))
names(no_can_genes) = c("1:10-20", "7:1300000-56800000", "7:74900000-148700000", "8:7400000-44300000", "8:47600000-146300000", "9:0-39200000", "9:69700000-140200000")
no_bg_genes = c(1,0)
names(no_bg_genes) = c("8:82000000-83000000", "21:1-3000000")

test_that("input_regions are checked - blocks",{
    expect_that(aba_enrich(overlap), throws_error("Background regions overlap: 2:5-20, 2:10-12"))
    expect_that(aba_enrich(overlap2), throws_error("Candidate regions overlap: 2:0-8, 2:5-10"))
    expect_that(aba_enrich(tight), throws_error("Background regions too small."))
    expect_that(aba_enrich(no_bg), throws_error("All values of the 'genes' input are 1. Using chromosomal regions as input requires defining background regions with 0."))
	expect_that(aba_enrich(reverse), throws_error("Invalid regions: 2:15-10.\n  In 'chr:start-stop' start < stop is required."))
	expect_that(aba_enrich(no_can_genes), throws_error("Candidate regions do not contain protein-coding genes."))
	expect_that(aba_enrich(no_bg_genes), throws_error("Background regions do not contain protein-coding genes."))
})    

test_that("input_regions are checked - circ_chrom",{
    expect_that(aba_enrich(overlap, circ_chrom=TRUE), throws_error("Background regions overlap: 2:5-20, 2:10-12"))
    expect_that(aba_enrich(overlap2, circ_chrom=TRUE), throws_error("Candidate regions overlap: 2:0-8, 2:5-10"))      
    expect_that(aba_enrich(no_bg, circ_chrom=TRUE), throws_error("All values of the 'genes' input are 1. Using chromosomal regions as input requires defining background regions with 0."))    
    expect_that(aba_enrich(reverse, circ_chrom=TRUE), throws_error("Invalid regions: 2:15-10.\n  In 'chr:start-stop' start < stop is required."))
})
