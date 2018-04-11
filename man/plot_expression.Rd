\name{plot_expression}
\alias{plot_expression}
\title{
    Plot expression data for given genes and brain structure IDs
}
\description{
    The function produces a heatmap (\code{\link{heatmap.2}} from package \code{gplots}) of gene expression in defined brain regions from adult or developing humans, or a developmental effect score for the developing human brain.
    Expression data is obtained from the Allen Brain Atlas project [1-4], averaged across donors, and for the developing human brain divided into five major age categories.
    If \code{gene_ids} and \code{dataset} are not specified, the genes and dataset from the last enrichment analysis with \code{\link{aba_enrich}} are used, since it may be a common case to first run the enrichment analysis and then look at the expression data.
    If a requested brain region has no expression data annotated, data from sampled substructures of this region is returned.
}

\usage{
    plot_expression(structure_ids, gene_ids = NA, dataset = NA, background = FALSE, 
    dendro = TRUE, age_category = 1)
}

\arguments{
    \item{structure_ids}{
        vector of brain structure ids, e.g. "Allen:10208".
    }
    \item{gene_ids}{
        vector of gene identifiers, either Entrez-ID, Ensembl-ID or gene-symbol. If not defined, genes from previous enrichment analysis with \code{\link{aba_enrich}} are used.
    }
    \item{dataset}{
        'adult' for the microarray dataset of adult human brains; '5_stages' for RNA-seq expression data of the developing human brain, grouped into 5 developmental stages;
        'dev_effect' for a developmental effect score. If not defined, dataset from last enrichment analysis with \code{\link{aba_enrich}} are used.
    }
    \item{background}{
        logical indicating whether expression from background genes should be included. Only used when \code{gene_ids} and \code{dataset} are \code{NA} so that genes from the last enrichment analysis with \code{\link{aba_enrich}} are used and when this analysis was performed using the hypergeometric test.  
    }
    \item{dendro}{
        logical indicating whether rows and columns should be rearranged with a dendrodram based on row/column means (using \code{\link{hclust}}).
        If \code{FALSE} and if \code{gene_ids} and \code{dataset} are \code{NA} so that genes from the last enrichment analysis with \code{\link{aba_enrich}} are used, a side bar indicates the gene-associated variables from the enrichment analysis:
        for a hypergeometric test the genes are grouped into candidate (red) and background genes (black);
        for the Wilcoxon rank-sum test, the side bar indicates the user-defined scores; for the binomial test the side bar shows *A/(A+B)* and for the 2x2 contingency table test *((A+1)/(B+1)) / ((C+1)/(D+1))* (+1 added to prevent division by 0, this is just a visual indication of the proportion of the ratios and not the real odds ratio from the 2x2 contingency table test)
    }
    \item{age_category}{
        an integer between 1 and 5 indicating the age category if \code{dataset='5_stages'}.
    }
}

\value{
    Invisibly, a list with components
    \item{rowInd}{row index permutation vector as returned by \code{\link{order.dendrogram}}}
    \item{colInd}{column index permutation vector.}    
    \item{call}{the matched call}
    \item{carpet}{reordered 'x' values used to generate the main 'carpet'}
    \item{rowDendrogram}{row dendrogram, if present}
    \item{colDendrogram}{column dendrogram, if present}
    \item{breaks}{values used for color break points}
    \item{col}{colors used}
    \item{colorTable}{A three-column data frame providing the lower and upper bound and color for each bin}
}

\references{
    [1] Hawrylycz, M.J. et al. (2012) An anatomically comprehensive atlas of the adult human brain transcriptome, Nature 489: 391-399. \doi{10.1038/nature11405}\cr
    [2] Miller, J.A. et al. (2014) Transcriptional landscape of the prenatal human brain, Nature 508: 199-206. \doi{10.1038/nature13185}\cr
    [3] Allen Institute for Brain Science. Allen Human Brain Atlas. Available from:
    \url{http://human.brain-map.org/}\cr
    [4] Allen Institute for Brain Science. BrainSpan Atlas of the Developing Human Brain. Available from: \url{http://brainspan.org/}
}
\author{
    Steffi Grote
}


\seealso{
    \code{browseVignettes("ABAEnrichment")}\cr
    \code{browseVignettes("ABAData")}\cr
    \code{\link{get_expression}}\cr
    \code{\link{aba_enrich}}\cr
    \code{\link{get_name}}\cr
    \code{\link{get_id}}\cr
    \code{\link{get_sampled_substructures}}\cr
    \code{\link{get_annotated_genes}}\cr
    \code{\link{heatmap.2}}\cr
    \code{\link{hclust}}\cr
}

\examples{
## plot expression data of six genes in two brain regions
## from children (age_category 3) without dendrogram
plot_expression(structure_ids=c("Allen:10657","Allen:10208"), 
    gene_ids=c("ENSG00000168036","ENSG00000157764","ENSG00000182158",
    "ENSG00000147889"), dataset="5_stages", dendro=FALSE, age_category=3)

## see the package vignette for more examples
}

