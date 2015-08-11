\name{get_expression}
\alias{get_expression}
\title{
Get expression data for given genes and brain structure ids
}
\description{
Expression data obtained from the Allen Brain Atlas project [1-4].   
}
\usage{
get_expression(structure_ids, gene_ids = NA, dataset = NA, background = FALSE)
}
\arguments{
  \item{structure_ids}{
vector of brain structure ids, e.g. 'Allen:10208'.
}
  \item{gene_ids}{
vector of gene identifiers, either Entrez-ID, Ensembl-ID or HGNC-symbol. If not defined, genes from previous enrichment analysis with \code{\link{aba_enrich}} are used.
}
  \item{dataset}{
'adult' for the microarray dataset of adult human brains; '5_stages' for RNA-seq expression data of the developing human brain, grouped into 5 developmental stages;
'dev_effect' for a developmental effect score. If not defined, dataset from last enrichment analysis with \code{\link{aba_enrich}} are used.
}
  \item{background}{
logical indicating whether expression from background genes should be included. Only used when \code{gene_ids} and \code{dataset} are \code{NA} and \code{test} from preceding \code{\link{aba_enrich}} call is 'hyper'.  
}
}
\details{
Get gene expression in defined brain regions from adult or developing humans, or a developmental effect score for the developing human brain. Expression data is obtained from the Allen Brain Atlas project [1-4], averaged across donors, and for the developing human brain divided into five major age categories. The developmental effect score is based on expression data of the developing human brain. If \code{gene_ids} and \code{dataset} are not specified, the genes and dataset from the last enrichment analysis with \code{\link{aba_enrich}} are used, since it may be a common case to first run the enrichment analysis and then look at the expression data. If a requested brain region has no expression data annotated, data from sampled substructures of this region is returned.\cr 
Please refer to the \code{\link[ABAData]{ABAData}} package vignette for details on the datasets.\cr
}
\value{

A matrix with expression values or developmental effect scores per brain region (rows) and gene (columns).\cr
For expression data from the developing human brain ('5_stages') it is a list with an expression matrix for each of the 5 developmental stages.\cr 


}
\references{
[1] Hawrylycz, M.J. et al. (2012) An anatomically comprehensive atlas of the adult human brain transcriptome, Nature 489: 391-399. \url{doi:10.1038/nature11405}\cr
[2] Miller, J.A. et al. (2014) Transcriptional landscape of the prenatal human brain, Nature 508: 199-206. \url{doi:10.1038/nature13185}\cr
[3] Allen Institute for Brain Science. Allen Human Brain Atlas [Internet]. Available from:
\url{http://human.brain-map.org/}\cr
[4] Allen Institute for Brain Science. BrainSpan Atlas of the Developing Human Brain
[Internet]. Available from: \url{http://brainspan.org/}

}
\author{
Steffi Grote
}


\seealso{
\code{vignette("ABAEnrichment",package="ABAEnrichment")}\cr
\code{vignette("ABAData",package="ABAData")}\cr
\code{\link{plot_expression}}\cr
\code{\link{aba_enrich}}\cr
\code{\link{get_name}}\cr
\code{\link{get_sampled_substructures}}\cr
}

\examples{ 
## get expression data for six genes in the brain structure 'Allen:4010'
get_expression(structure_ids=c('Allen:4010'),gene_ids=c(324,8312,673,1029,64764,1499),dataset='adult')
## get expression data of six genes in two brain regions from developing human brain, each of the five list elements corresponds to an age category 
get_expression(structure_ids=c('Allen:10657','Allen:10208'), gene_ids=c('ENSG00000168036', 'ENSG00000157764', 'ENSG00000182158', 'ENSG00000147889'),dataset='5_stages')
}

