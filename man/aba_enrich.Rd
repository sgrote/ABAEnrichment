\name{aba_enrich}
\alias{aba_enrich}
\title{
Test genes for expression enrichment in human brain regions
}
\description{
Tests for enrichment of user defined candidate genes in the set of expressed protein-coding genes in different human brain regions. It integrates the expression of the candidate gene set (averaged across donors) and the structural information of the brain using an ontology, both provided by the Allen Brain Atlas project [1-4]. The statistical analysis is performed using the ontology enrichment software FUNC [5].
}
\usage{
aba_enrich(genes, dataset = 'adult', test = 'hyper', 
  cutoff_quantiles = seq(0.1, 0.9, 0.1), n_randsets = 1000, 
  gene_len = FALSE, circ_chrom = FALSE, ref_genome = 'grch37')
}
\arguments{
  \item{genes}{
If \code{test='wilcoxon'} a numeric vector of scores.
If \code{test='hyper'} (default) a binary vector with 1 for candidate genes and 0 for background
genes. If no background genes are defined, all remaining protein-coding genes are used as background.
The names of the vector are the gene identifiers: either Entrez-ID, Ensembl-ID or gene-symbol. 
For \code{test='hyper'} the names of the vector can also describe chromosomal regions ('chr:start-stop').
}
  \item{dataset}{
'adult' for the microarray dataset of adult human brains; '5_stages' for RNA-seq expression data for different stages of the developing human brain, grouped into 5 developmental stages;
'dev_effect' for a developmental effect score. For details see \code{vignette("ABAData", package="ABAData")}.
}
  \item{test}{
'hyper' (default) for the hypergeometric test or 'wilcoxon' for the Wilcoxon rank test.
}
  \item{cutoff_quantiles}{
the FUNC enrichment analyses will be performed for the sets of expressed genes at given expression quantiles defined in this vector [0,1].
}
  \item{n_randsets}{
integer defining the number of random sets created to compute the FWER.
}  
  \item{gene_len}{
logical. If \code{test='hyper'} the probability of a background gene to be chosen as a candidate gene in a random set is dependent on the gene length. 
}
  \item{circ_chrom}{
logical. When \code{genes} defines chromosomal regions, \code{circ_chrom=TRUE} uses background regions from the same chromosome and allows randomly chosen blocks to overlap multiple background regions. Only if \code{test='hyper'}.
}
  \item{ref_genome}{
'grch37' (default) or 'grch38'. Defines the reference genome used when genomic regions are provided as input or when \code{gene_len=TRUE}.
}

}
\details{
For details please refer to \code{vignette("ABAEnrichment", package="ABAEnrichment")}.
}
\value{
A list with components
\item{results}{
a dataframe with the FWERs from the enrichment analyses per brain region and age category, ordered by 'age_category', 'times_FWER_under_0.05','min_FWER' and 'mean_FWER'; with 'min_FWER' for example denoting the minimum FWER for expression enrichment of the candidate genes in this brain region across all expression cutoffs. 'FWERs' is a semicolon separated string with the single FWERs for all cutoffs. 'equivalent_structures' is a semicolon separated string that lists structures with identical expression data due to lack of independent expression measurements in all regions.}
\item{genes}{
a vector of the requested genes, excluding those genes for which no expression data is available and which therefore were not included in the enrichment analysis.}    
\item{cutoffs}{
a dataframe with the expression values that correspond to the requested cutoff quantiles.}   
}

\references{
[1] Hawrylycz, M.J. et al. (2012) An anatomically comprehensive atlas of the adult human brain transcriptome, Nature 489: 391-399. \doi{10.1038/nature11405}\cr
[2] Miller, J.A. et al. (2014) Transcriptional landscape of the prenatal human brain, Nature 508: 199-206. \doi{10.1038/nature13185}\cr
[3] Allen Institute for Brain Science. Allen Human Brain Atlas [Internet]. Available from:
\url{http://human.brain-map.org/}\cr
[4] Allen Institute for Brain Science. BrainSpan Atlas of the Developing Human Brain
[Internet]. Available from: \url{http://brainspan.org/}\cr
[5] Pruefer, K. et al. (2007) FUNC: A package for detecting significant associations between gene
sets and ontological, BMC Bioinformatics 8: 41. \doi{10.1186/1471-2105-8-41}
}
\author{
Steffi Grote
}

\seealso{
\code{vignette("ABAEnrichment", package="ABAEnrichment")}\cr
\code{vignette("ABAData", package="ABAData")}\cr
\code{\link{get_expression}}\cr
\code{\link{plot_expression}}\cr
\code{\link{get_name}}\cr
\code{\link{get_id}}\cr
\code{\link{get_sampled_substructures}}\cr
\code{\link{get_superstructures}}\cr
\code{\link{get_annotated_genes}}\cr

}
\examples{

#### Note that arguments 'cutoff_quantiles' and 'n_randsets' are reduced 
#### to lower computational time in the examples. 
#### Using the default values is recommended.

#### Perform an enrichment analysis for the developing brain
#### with defined background genes
#### and with random sets dependent on gene length
genes = rep(c(1,0),each=7)
names(genes) = c('PENK', 'COCH', 'PDYN', 'CA12', 'SYNDIG1L', 'MME', 
  'ANO3', 'KCNJ6', 'ELAVL4', 'BEAN1', 'PVALB', 'EPN3', 'PAX2', 'FAB12')
res = aba_enrich(genes, dataset='5_stages', cutoff_quantiles=c(0.5,0.9), 
  n_randsets=100, gene_len=TRUE)
## see results for the brain regions with highest enrichment for children (age_category 3)
fwers = res[[1]]
head(fwers[fwers$age_category==3,])
## see the input genes vector (only genes with expression data available) 
res[2]
## see the expression values that correspond to the requested cutoff quantiles
res[3]


# For more examples please refer to the package vignette.

}

\keyword{ htest }
