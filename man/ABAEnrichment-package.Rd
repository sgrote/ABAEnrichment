\name{ABAEnrichment-package}
\alias{ABAEnrichment-package}
\alias{ABAEnrichment}
\docType{package}
\title{
Gene expression enrichment in human brain regions
}
\description{
The package ABAEnrichment is designed to test for enrichment of user defined candidate genes in the set of expressed genes in different human brain regions. The package integrates the expression of the candidate gene set (averaged across donors) and the structural information of the brain using an ontology, both provided by the Allan Brain Atlas project [1-4]. The statistical analysis is performed by the core function \code{\link{aba_enrich}} which interfaces the ontology enrichment software FUNC [5]. Additional functions provided in this package like \code{\link{get_expression}} and \code{\link{plot_expression}} facilitate exploring the expression data.
}
\details{
\tabular{ll}{
Package: \tab ABAEnrichment\cr
Type: \tab Package\cr
Version: \tab 0.99.4\cr
Date: \tab 2015-08-06\cr
License: \tab GPL (>= 2)\cr
}

For details see \code{vignette("ABAEnrichment",package="ABAEnrichment")}\cr

% eventuell noch die passage aus der Documentation.pdf einfuegen einfuegen die so beginnt:

%The statistical tests are performed by using either a hypergeometric test or a wilcoxon rank test [5]. In addition to gene expression the enrichment may refer to a developmental effect score, which describes how much a gene's expression changes over time. The user can choose between three different datasets: first, the developmental effect score, second, microarray data from adult donors and ....
}

\author{
Steffi Grote\cr
Maintainer: Steffi Grote <steffi_grote@eva.mpg.de>
}
\references{
[1] Hawrylycz, M.J. et al. (2012) An anatomically comprehensive atlas of the adult human brain transcriptome, Nature 489: 391-399. \url{doi:10.1038/nature11405}\cr
[2] Miller, J.A. et al. (2014) Transcriptional landscape of the prenatal human brain, Nature 508: 199-206. \url{doi:10.1038/nature13185}\cr
[3] Allen Institute for Brain Science. Allen Human Brain Atlas [Internet]. Available from:
\url{http://human.brain-map.org/}\cr
[4] Allen Institute for Brain Science. BrainSpan Atlas of the Developing Human Brain
[Internet]. Available from: \url{http://brainspan.org/}\cr
[5] Pruefer, K. et al. (2007) FUNC: A package for detecting significant associations between gene
sets and ontological, BMC Bioinformatics 8: 41. \url{doi:10.1186/1471-2105-8-41}
}
\keyword{ package }
\keyword{ htest }
\seealso{
\code{vignette("ABAEnrichment",package="ABAEnrichment")}\cr
\code{vignette("ABAData",package="ABAData")}\cr
\code{\link{aba_enrich}}\cr
\code{\link{get_expression}}\cr 
\code{\link{plot_expression}}\cr
\code{\link{get_name}}\cr
\code{\link{get_sampled_substructures}}\cr
\code{\link{get_superstructures}}\cr
}

\examples{

#### Perform gene expression enrichment analysis on 13 test genes in five developmental stages of the human brain using the hypergeometric test implemented in FUNC[5]   
## create input vector with test genes 
genes=rep(1,13)
names(genes)=c('NCAPG', 'APOL4', 'NGFR', 'NXPH4', 'C21orf59', 'CACNG2', 'AGTR1', 'ANO1', 'BTBD3', 'MTUS1', 'CALB1', 'GYG1', 'PAX2')
## run enrichment analysis
res=aba_enrich(genes,dataset='5_stages',cutoff_quantiles=c(0.5,0.7,0.9))

## get expression data (list) for the first 10 brain structures of the enrichment analysis output 
get_expression(res$results$structure_id[1:10],background=FALSE)
## plot this expression with dendrogram
plot_expression(res$results$structure_id[1:10],dendro=TRUE,age_category=1)
## plot expression from the 10 highest scoring brain structures in age category 3 without dendrogram
plot_expression(res$results[res$results$age_category==3,'structure_id'][1:10],dendro=FALSE,age_category=3)
## get the full name of a structure given the structure id
get_name('Allen:10657')
## get ids of the substructers of 'Allen:10208' that have expression data annotated
get_sampled_substructures('Allen:10208')
## get ids of the superstructers of 'Allen:10208' 
get_superstructures('Allen:10208')

#### Perform gene expression enrichment analysis on 15 test genes in the adult human brain using the wilcoxon rank test implemented in FUNC[5] 
## create input vector with random scores associated with the test genes (using entrez-Ids)
genes=sample(1:50,15)
names(genes)=c(324,8312,673,1029,64764,1499,3021,3417,3418,8085,3845,9968,5290,5727,5728)
## run enrichment analysis
res=aba_enrich(genes,dataset='adult',test='wilcoxon',cutoff_quantiles=c(0.2,0.5,0.8))

## plot expression from the 10 highest scoring brain structures without dendrogram
plot_expression(res$results$structure_id[1:10],dendro=FALSE)


#### Explore expression data idependently from enrichment analysis
## get expression data for six genes in the brain structure 'Allen:4010' 
get_expression(structure_ids=c('Allen:4010'),gene_ids=c(324,8312,673,1029,64764,1499),dataset='adult')
## plot this expression
plot_expression(structure_ids=c('Allen:4010'),gene_ids=c(324,8312,673,1029,64764,1499),dataset='adult')
## 'Allen:4010' apparently was not sampled directly - get the names of the substructures that provide expression data for 'Allen:4010'
get_name(get_sampled_substructures('Allen:4010'))
## plot expression data of six genes in two brain regions from children (3-11 yrs, age_category 3)
plot_expression(structure_ids=c('Allen:10657','Allen:10208'), gene_ids=c('ENSG00000168036', 'ENSG00000157764', 'ENSG00000182158', 'ENSG00000147889', 'ENSG00000103126', 'ENSG00000134982'),dataset='5_stages',dendro=FALSE, age_category=3)
## get developmental effect scores of six genes in two brain regions
get_expression(structure_ids=c('Allen:10657','Allen:10208'), gene_ids=c('ENSG00000168036', 'ENSG00000157764', 'ENSG00000182158', 'ENSG00000147889'),dataset='dev_effect')

}
