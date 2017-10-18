
# context("get genes annotated to brain regions")

#get_annotated_genes = function(res, fwer_threshold=0.05, background=FALSE, structure_ids=NULL, cutoff_quantiles=NULL, dataset=NULL, genes=NULL)


## aba_enrich result as input
#############################

# dev, hyper, hgnc
test_genes = paste(rep('FABP', 5), c(2,4:7), sep='')
bg_genes = c('NCAPG', 'NGFR', 'NXPH4', 'C21orf59', 'CACNG2', 'AGTR1', 'ANO1', 'BTBD3', 'MTUS1', 'MIRLET7BHG')
genes = c(rep(1,length(test_genes)), rep(0,length(bg_genes)))
names(genes) = c(test_genes, bg_genes)
set.seed(123)
res_dev = aba_enrich(genes,n_randsets=5, dataset='5_stages', cutoff_quantiles=c(0.3,0.9999), silent=TRUE)

test_that("warning and return(NULL) when no node below fwer-threshold",{
    expect_that(get_annotated_genes(res_dev), gives_warning("No significantly enriched brain regions at FWER-threshold 0.05"))
    expect_that(get_annotated_genes(res_dev), equals(NULL))
})

## TODO: example with working FWER-threshold (at least that it's a dataframe with x columns)

## TODO: independent input





