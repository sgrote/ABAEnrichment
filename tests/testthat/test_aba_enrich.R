
# context("main function aba_enrich")
#aba_enrich=function(genes,dataset="adult",test="hyper",cutoff_quantiles=seq(0.1,0.9,0.1),n_randsets=1000, gene_len=FALSE, circ_chrom=FALSE, ref_genome="grch37")


## genes input
##############

genes=rep(0:1,4)
names(genes)=c(324,8312,673,1029,64764,1499,"abc",000)

set.seed(123)
res1=aba_enrich(genes,dataset="adult",cutoff_quantiles=c(0.2,0.7),n_randsets=30)

test_that("normal input results in list",{
	expect_that(class(res1), equals("list"))
	expect_that(nrow(res1$cutoffs), equals(2))
	expect_that(length(unlist(strsplit(res1[[1]][1,8],";"))), equals(2))
})

res2=aba_enrich(genes,dataset="adult",cutoff_quantiles=c(0.999,0.2),n_randsets=5)

test_that("cutoffs get sorted, output for successful cutoffs is returned",{
	expect_that(class(res2), equals("list"))
	expect_that(length(unlist(strsplit(res2[[1]][1,8],";"))), equals(1))
	expect_that(rownames(res2[[3]]), equals(c("20%","99.9%")))
})

willi = as.integer(runif(length(genes),1,30))
names(willi) = names(genes)

set.seed(123)
res3=aba_enrich(willi,dataset="5_stages",test="wilcoxon",cutoff_quantiles=c(0.4),n_randsets=30)

test_that("corner case 5_stages with one cutoff works",{
	expect_that(nrow(res3$cutoffs), equals(1))
})

test_that("genes that are not in data are not returned in the genes object",{
	expect_that("abc" %in% names(res3$genes), is_false())	
	expect_that("8312" %in% names(res3$genes), is_true())
    expect_that(length(res3$genes), equals(6))
})


# test warnings and more on handling of no coords and data

test_genes = c(paste(rep("FABP", 5), c(2,4:7), sep=""), c("HSPB2" ,"LINC00239", "TESTI"))
bg_genes = c("NCAPG", "NGFR", "NXPH4", "C21orf59", "CACNG2", "AGTR1", "ANO1", "BTBD3", "MTUS1", "MIRLET7BHG", "MUELL")
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

# "LINC00239", "MIRLET7BHG" -> no coordinates
test_that("Warning about test and background genes with no coordinates.",{
    expect_that(aba_enrich(genes,cutoff_quantiles=c(0.7,0.8),n_randsets=5,gene_len=TRUE), gives_warning("No coordinates available for genes: LINC00239, MIRLET7BHG.\n  These genes were not included in the analysis."))    
})

test_that("One candidate or background gene without coordinates stops.",{
    expect_that(aba_enrich(genes[7:length(genes)],cutoff_quantiles=c(0.7,0.8),n_randsets=5,gene_len=T), throws_error("No requested test genes in data."))    
    expect_that(aba_enrich(genes[c(1:8,18,19)],cutoff_quantiles=c(0.7,0.8),n_randsets=5,gene_len=T), throws_error("No requested background genes in data."))    
})   

set.seed(123)
res4 = aba_enrich(genes,n_randsets=5,gene_len=TRUE)
test_that("genes that are not in data or coordinates are not returned in the genes object",{
	expect_that(length(res4$genes), equals(length(genes)-4))
	expect_that("NGFR" %in% names(res4$genes), is_true())
	expect_that("HSPB2" %in% names(res4$genes), is_true()) # not in hgnc-symbol
	expect_that("TESTI" %in% names(res4$genes), is_false())	
	expect_that("MUELL" %in% names(res4$genes), is_false())	
	expect_that("LINC00239" %in% names(res4$genes), is_false())	
	expect_that("MIRLET7BHG" %in% names(res4$genes), is_false())	

})

# NEW: grch38
set.seed(123)
res4b = aba_enrich(genes,n_randsets=5,gene_len=TRUE,ref_genome="grch38")
test_that("genes that are not in data or coordinates are not returned in the genes object",{
	expect_that(length(res4b$genes), equals(length(genes)-2))
	expect_that("NGFR" %in% names(res4b$genes), is_true())
	expect_that("HSPB2" %in% names(res4$genes), is_true()) # not in hgnc-symbol
	expect_that("TESTI" %in% names(res4b$genes), is_false())	
	expect_that("MUELL" %in% names(res4b$genes), is_false())	
	expect_that("LINC00239" %in% names(res4b$genes), is_true()) # not in hg19	
	expect_that("MIRLET7BHG" %in% names(res4b$genes), is_true())# not in hg19
})

## NEW test unused arguments
test_that("Warning about unuesed arguments.",{
    expect_that(aba_enrich(genes,test="wilcoxon",cutoff_quantiles=c(0.7,0.8),n_randsets=5,gene_len=TRUE), gives_warning("Unused argument: 'gene_len = TRUE'."))    
    expect_that(aba_enrich(genes,cutoff_quantiles=c(0.7,0.8),n_randsets=5,ref_genome="grch38"), gives_warning("Unused argument: 'ref_genome = grch38'."))    
})

genes = c(rep(1,6),rep(0,4))
genes_willi = 1:10
names(genes) = c("NCAPG", "APOL4", "NGFR", "NXPH4", "C21orf59", "CACNG2", "AGTR1", "ANO1","APOL4", "NXPH4")
names(genes_willi) = c("NCAPG", "APOL4", "NGFR", "NGFR", "NXPH4", "C21orf59", "AGTR1", "CACNG2", "AGTR1", "ANO1")
test_that("multiple assignment of scores to genes throws error",{
    expect_that(aba_enrich(genes), throws_error("Genes with multiple assignment in input: APOL4, NXPH4"))
    expect_that(aba_enrich(genes_willi, test="wilcoxon"), throws_error("Genes with multiple assignment in input: AGTR1, NGFR"))
})



## full run using all genes as background  
#genes = rep(1,13)
#names(genes) = c("NCAPG", "APOL4", "NGFR", "NXPH4", "C21orf59", "CACNG2", "AGTR1", "ANO1",
#"BTBD3", "MTUS1", "CALB1", "GYG1", "PAX2")
#set.seed(123)
##res5 = aba_enrich(genes,dataset = "5_stages",cutoff_quantiles = c(0.5,0.7,0.9), n_randsets = 100)
#set.seed(123)
#res6 = aba_enrich(genes,dataset = "5_stages",cutoff_quantiles = c(0.5,0.7,0.9), n_randsets = 100, gene_len = TRUE)

#fwers1 = res1[[1]]
#fwers3 = res3[[1]]
#fwers4 = res4[[1]]
#fwers4b = res4b[[1]] #NEW
##fwers5 = res5[[1]]
#fwers6 = res6[[1]]
#test_that("check some FWERs - single genes.",{
#	expect_that(round(fwers4[fwers4$structure_id=="Allen:4124","mean_FWER"],7), equals(0.8155556))
#	expect_that(round(mean(fwers4$mean_FWER),7), equals(0.7959293))
#	expect_that(round(mean(fwers4$min_FWER),7), equals(0.4330594))	
	
#	# NEW
#	expect_that(round(fwers4b[fwers4b$structure_id=="Allen:4124","mean_FWER"],7), equals(0.8744444))
#	expect_that(round(mean(fwers4b$mean_FWER),7), equals(0.8469711))
#	expect_that(round(mean(fwers4b$min_FWER),7), equals(0.5332725))	
	
#	expect_that(round(mean(fwers3$mean_FWER),7), equals(0.9994815))	
	
#	expect_that(fwers1[fwers1$structure_id=="Allen:9150","mean_FWER"], equals(0.965))
#	expect_that(round(mean(fwers1$mean_FWER),7), equals(0.9876941))
#	expect_that(round(mean(fwers1$min_FWER),7), equals(0.9753881))

##	expect_that(fwers5[1,8], equals("0.7;0.34;0.24")) # ("0.7;0.35;0.24" in standard testing environment)
#	expect_that(fwers6[1,8], equals("0.37;0.25;0.11"))
#})


## genomic region input
#######################

### normal runs (examples similar to vignette)
#genes1 = c(1, rep(0,6))
#names(genes1) = c("8:81000000-83000000", "7:1300000-56800000", "7:74900000-148700000",
# "8:7400000-44300000", "8:47600000-146300000", "9:0-39200000", "9:69700000-140200000")
#set.seed(123)
#res1 = aba_enrich(genes1, n_randsets=100) 

#genes2 = c(1,1,rep(0,4))
#names(genes2) = c("8:81000000-83000000","7:113600000-124700000","7:1300000-56800000", "7:74900000-148700000",
# "8:7400000-44300000", "8:47600000-146300000")
#set.seed(123)
#res2 = aba_enrich(genes2, n_randsets=100, circ_chrom=TRUE) 

### NEW: hg20
#set.seed(123)
#res1b = aba_enrich(genes1, n_randsets=100, ref_genome="grch38") 
#set.seed(123)
#res2b = aba_enrich(genes2, n_randsets=100, circ_chrom=TRUE, ref_genome="grch38") 

#fwers1 = res1[[1]]
#fwers2 = res2[[1]]
#fwers1b = res1b[[1]]
#fwers2b = res2b[[1]]

#test_that("check some FWERs - genomic regions.",{
#	expect_that(round(fwers1[fwers1$structure_id=="Allen:4735","mean_FWER"],7), equals(0.5277778))
#	expect_that(round(fwers1b[fwers1b$structure_id=="Allen:4735","mean_FWER"],7), equals(0.6166667))
#	expect_that(fwers1[fwers1$structure_id=="Allen:4743","FWERs"], equals("0.98;0.74;0.88;0.46;0.43;0.34;0.37;0;0.1"))
#	expect_that(fwers1b[fwers1b$structure_id=="Allen:4743","FWERs"], equals("0.97;0.84;0.99;0.68;0.61;0.48;0.64;0.07;0.14"))
#	expect_that(round(mean(fwers1$mean_FWER),7), equals(0.6313699))
#	expect_that(round(mean(fwers1b$mean_FWER),7), equals(0.7154997))
#	expect_that(round(mean(fwers1$min_FWER),7), equals(0.0472603))
#	expect_that(round(mean(fwers1b$min_FWER),7), equals(0.1219635))
#	expect_that(round(fwers2[fwers2$structure_id=="Allen:4598","mean_FWER"],7), equals(0.9022222))
#	expect_that(round(fwers2b[fwers2b$structure_id=="Allen:4598","mean_FWER"],7), equals(0.9166667))
#	expect_that(fwers2[fwers2$structure_id=="Allen:4743","FWERs"], equals("1;1;1;1;0.98;0.97;0.99;0.65;0.84"))
#	expect_that(fwers2b[fwers2b$structure_id=="Allen:4743","FWERs"], equals("1;1;1;1;1;0.98;1;0.85;0.9"))
#	expect_that(round(mean(fwers2$mean_FWER),7), equals(0.9648537))
#	expect_that(round(mean(fwers2b$mean_FWER),7), equals(0.9824505))
#	expect_that(round(mean(fwers2$min_FWER),7), equals(0.7839726))
#	expect_that(round(mean(fwers2b$min_FWER),7), equals(0.8764079))
#})


# test checking: background < candidate on chrom, in blocks, overlapping input regions
too_small = c(rep(1,4),rep(0,5)) 
names(too_small) = c("X:0-2","2:0-3","2:5-10","13:0-20",  "2:5-10","13:0-10","X:0-1","4:0-10","2:10-12") 
overlap = too_small
names(overlap) = c("X:0-1","2:0-3","2:5-10","13:0-20",  "2:5-20","13:0-30","X:0-1","4:0-10","2:10-12")
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
