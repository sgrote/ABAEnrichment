
# context("get genes annotated to brain regions")

#get_annotated_genes = function(res, fwer_threshold=0.05, background=FALSE, structure_ids=NULL, cutoff_quantiles=NULL, dataset=NULL, genes=NULL)


## aba_enrich result as input
#############################

# adult, hyper, hgnc
test_genes = paste(rep('FABP', 5), c(2,4:7), sep='')
bg_genes = c('NCAPG', 'NGFR', 'NXPH4', 'C21orf59', 'CACNG2', 'AGTR1', 'ANO1', 'BTBD3', 'MTUS1', 'MIRLET7BHG')
genes = c(rep(1,length(test_genes)), rep(0,length(bg_genes)))
names(genes) = c(test_genes, bg_genes)
set.seed(123)
res_adult = aba_enrich(genes,n_randsets=100, cutoff_quantiles=c(0.9,0.3,0.5,0.7,0.9999))

test_that("warning and return(NULL) when no node below fwer-threshold",{
    expect_that(get_annotated_genes(res_adult), gives_warning("No significantly enriched brain regions at FWER-threshold 0.05"))
    expect_that(get_annotated_genes(res_adult), equals(NULL))
}) 

x = get_annotated_genes(res_adult, fwer_threshold=0.5)
test_that("result for res_adult as expected",{
    expect_that(unique(x[,1]), equals(factor(5)))
    expect_that(unique(x[,2]), equals(c("Allen:4263","Allen:4264")))
    expect_that(unique(x[,3]), equals(0.5))
    expect_that(unique(x[,4]), equals(c("FABP5","FABP6","FABP7")))
    expect_that(unique(x[,5]), equals(0.49))
    expect_that(unique(x[,6]), equals(1))
}) 

y = get_annotated_genes(res_adult, background=TRUE, fwer_threshold=0.5)
test_that("result for res_adult as expected - background=TRUE",{
    expect_that(nrow(y), equals(12))
    expect_that(unique(y[,6]), equals(c(0,1)))
}) 

# 5-stages, wilcox, entrez
set.seed(123)
genes = sample(1:50,15)
names(genes) = c(324,8312,673,1029,64764,1499,3021,3417,3418,8085,3845,9968,5290,5727,5728)
res_5 = aba_enrich(genes, dataset="5_stages", test="wilcoxon", cutoff_quantiles=c(0.2,0.5,0.8,0.9,0.95), n_randsets=100)

test_that("warning and return(NULL) when no node below fwer-threshold - 5_stages",{
    expect_that(get_annotated_genes(res_5), gives_warning("No significantly enriched brain regions at FWER-threshold 0.05"))
    expect_that(get_annotated_genes(res_5), equals(NULL))
}) 

z = get_annotated_genes(res_5, fwer_threshold=0.55)
test_that("result for res_adult as expected",{
    expect_that(as.numeric(unique(z[,1])), equals(c(4,1)))
    expect_that(unique(z[,3]), equals(0.8))
    expect_that(unique(z[,5]), equals(c(0.42, 0.53)))
    expect_that(unique(z[,6]), equals(c(3,15,24,46,43,49)))
}) 


## structures, dataset, cutoffs defined by user
###############################################
x0 = get_annotated_genes(structure_ids=unique(x$structure_id)[1:2], cutoff_quantiles=0.5, dataset="adult") 
x1 = get_annotated_genes(structure_ids=unique(x$structure_id)[1:2], cutoff_quantiles=0.5, dataset="adult",genes=c(test_genes)) 
y1 = get_annotated_genes(structure_ids=unique(x$structure_id)[1:2], cutoff_quantiles=0.5, dataset="adult", genes=c(test_genes, bg_genes)) 
z1 = get_annotated_genes(structure_ids=unique(z$structure_id), cutoff_quantiles=c(0.2,0.8), dataset="5_stages", genes=names(genes))

test_that("check nrow of result with user-defined structures, cutoffs, dataset, and genes",{
    expect_that(nrow(x0), equals(14175))
    expect_that(all(x1==x[,1:4]), is_true())
    expect_that(all(table(y1$anno_gene)==table((y$anno_gene))), is_true())
    expect_that(unique(z1[,3]), equals(c(0.2,0.8)))
}) 












