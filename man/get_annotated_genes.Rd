\name{get_annotated_genes}
\alias{get_annotated_genes}
\title{
Get genes that are expressed in enriched or user-defined brain regions
}
\description{
Uses an object returned from \code{\link{aba_enrich}} as input and returns the brain regions that are significantly (given a FWER-threshold) enriched in expressed candidate genes (or genes with high scores if \code{test='wilcoxon'} in \code{\link{aba_enrich}}) together with the genes that are expressed in those brain regions (i.e. are 'annotated' to the brain regions).
Alternatively, also user-defined brain regions, dataset and expression cutoffs can be used as input.
}

\usage{
get_annotated_genes(res, fwer_threshold = 0.05, background = FALSE, 
  structure_ids = NULL, dataset = NULL, cutoff_quantiles = NULL, genes = NULL)
}
\arguments{
  \item{res}{
an object returned from \code{\link{aba_enrich}} (list of 3 elements). If not defined, \code{structure_ids}, \code{dataset} and \code{cutoff_quantiles} have to be specified. 
}
  \item{fwer_threshold}{
numeric defining the FWER-threshold. Given \code{res} as input, \code{get_annotated_genes} extracts all brain-region/expression-cutoff combinations from \code{res} that have a FWER < \code{fwer_threshold} and adds the (candidate) genes that are annotated to those brain regions at the given expression cutoffs.  
}
  \item{background}{
logical indicating whether background genes should be included. Only used when \code{res} is defined and contains the results from a hypergeometric test (which is the default in \code{\link{aba_enrich}}). 
}
  \item{structure_ids}{
vector of brain structure IDs, e.g. 'Allen:10208'. If \code{res} is not defined, \code{structure_ids} specifies the brain regions for which annotated genes at the given \code{cutoff_quantiles} will be returned.
}
  \item{dataset}{
'adult' for the microarray dataset of adult human brains; '5_stages' for RNA-seq expression data of the developing human brain, grouped into 5 developmental stages; 'dev_effect' for a developmental effect score.
Only used when \code{res} is not defined.
}
  \item{cutoff_quantiles}{
vector of numeric values between 0 and 1. They define the expression quantiles (across all genes) which are used as cutoffs to decide whether a gene counts as expressed (and gets annotated to a brain region) or not.
Only used when \code{res} is not defined.
}
  \item{genes}{
optional vector of gene identifiers, either Entrez-ID, Ensembl-ID or gene-symbol. If defined, only annotations of those genes are returned. If not defined, all expressed genes from Allen Brain Atlas are returned. Only used when \code{res} is not defined.
}
}

\details{
Genes get annotated to a brain region when their expression value in that brain region, which is provided by the Allen Brain Atlas, exceeds a certain cutoff. Multiple cutoffs can be used. They are defined as quantiles of gene expression across all genes and brain regions. An expression cutoff of e.g. 0.8 means that only genes with expression levels higher than 80\% of all measured values count as 'expressed' and get annotated to the corresponding brain region.
Note that those annotations are inherited by all superstructures (parent nodes).\cr
When the 'dev_effect' dataset is used as \code{dataset}-argument, the expression value is replaced by a measure of how much a gene's expression changes during the development (prenatal to adult). Genes with a developmental score higher than the cutoff then get annotated to the brain regions. 
}

\value{
A data frame with the age category, the brain region ID, the expression cutoff quantile and the annotated genes for significantly enriched or user-defined brain-region/expression-cutoff combinations. If \code{res} was given as input, two additional columns are added: the FWER and the score which was used for the genes in the \code{\link{aba_enrich}} input (1/0 for candidate and background genes for the hypergeometric test or scores for the wilcoxon rank sum test). Given \code{res} as input, the output is ordered by the FWER, else it is ordered by the expression cutoff.\cr 
}

\references{
[1] Hawrylycz, M.J. et al. (2012) An anatomically comprehensive atlas of the adult human brain transcriptome, Nature 489: 391-399. \doi{10.1038/nature11405}\cr
[2] Miller, J.A. et al. (2014) Transcriptional landscape of the prenatal human brain, Nature 508: 199-206. \doi{10.1038/nature13185}\cr
[3] Allen Institute for Brain Science. Allen Human Brain Atlas [Internet]. Available from:
\url{http://human.brain-map.org/}\cr
[4] Allen Institute for Brain Science. BrainSpan Atlas of the Developing Human Brain
[Internet]. Available from: \url{http://brainspan.org/}

}
\author{
Steffi Grote
}


\seealso{
\code{browseVignettes("ABAEnrichment")}\cr
\code{browseVignettes("ABAData")}\cr
\code{\link{aba_enrich}}\cr
\code{\link{get_name}}\cr
\code{\link{get_id}}\cr
}

\examples{ 
#### Note that arguments 'cutoff_quantiles' and 'n_randsets' are reduced 
#### to lower computational time in the examples.

## perform an enrichment analysis
## for expression of 7 candidate and 7 background genes
## and get candidate genes annotated to brain regions that have a FWER < 0.05
set.seed(123)
gene_ids = c('PENK', 'COCH', 'PDYN', 'CA12', 'SYNDIG1L', 'MME', 
    'ANO3', 'KCNJ6', 'ELAVL4', 'BEAN1', 'PVALB', 'EPN3', 'PAX2', 'FAB12')
is_candidate = rep(c(1,0), each=7)
genes = data.frame(gene_ids, is_candidate)
res = aba_enrich(genes, dataset='5_stages', cutoff_quantiles=c(0.3,0.5,0.7,0.9),
    n_randset=100)
anno = get_annotated_genes(res, fwer_threshold=0.05)
head(anno)

## find out which of the above genes have expression above 
## the 50% and 70% expression-cutoff, respectively,
## in the Cerebellar Cortex of the developing human brain (Allen:10657)
get_annotated_genes(structure_ids="Allen:10657", dataset="5_stages",
    cutoff_quantiles=c(0.5,0.7), genes=gene_ids)
}

