### ABAEnrichment: Gene Expression Enrichment in Human Brain Regions 

[_ABAEnrichment_](https://www.bioconductor.org/packages/release/bioc/html/ABAEnrichment.html)
is an R-package designed to test user-defined genes for expression enrichment in different human brain regions.
The package integrates the expression of the input gene set and the structural information of the brain using an ontology, both provided by the [_Allen Brain Atlas project_](http://www.brain-map.org/).
The statistical analysis is performed by the core function `aba_enrich` which interfaces with the ontology enrichment software [_FUNC_](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1800870/).
Additional functions provided in this package are `get_expression`, `plot_expression`, `get_name`, `get_id`, `get_sampled_substructures`, `get_superstructures` and `get_annotated_genes` supporting the exploration and visualization of the expression data.

#### Installation
A stable release version can be obtained from [_Bioconductor_](https://www.bioconductor.org/packages/release/bioc/html/ABAEnrichment.html).

+ Installation from Bioconductor
```r
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("ABAEnrichment")
```

The developmental (this) version can be obtained from the ['devel' version of Bioconductor](https://bioconductor.org/developers/how-to/useDevel/) or directly from
GitHub:

+ Installation from GitHub

```r
if (!requireNamespace("devtools", quietly=TRUE))
    install.packages("devtools")
devtools::install_github("sgrote/ABAEnrichment")
```



#### Usage  
See the package's [vignette](https://www.bioconductor.org/packages/release/bioc/vignettes/ABAEnrichment/inst/doc/ABAEnrichment.html) for a tutorial.  
Also refer to the man-pages for single functions, e.g.

```r
library(ABAEnrichment)
?aba_enrich
?get_expression
```


