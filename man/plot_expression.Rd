\name{plot_expression}
\alias{plot_expression}
\title{
    Plot a heatmap given expression data
}
\description{
    The function produces a (\code{\link{heatmap.2}} from package \code{gplots}). To be used with gene-expression matrices obtained from \code{\link{get_expression}}.
}

\usage{
    plot_expression(expr_mat, dendro = TRUE, gene_vars = NULL, main = "")
}

\arguments{
    \item{expr_mat}{
        matrix of expression values with genes in columns and brain regions in rows, like returned from \code{\link{get_expression}}
        (rownames are expected to be brain region identifiers, and their acronym is added automatically)
    }
    \item{dendro}{
        logical indicating whether rows and columns should be rearranged with a dendrogram based on row/column means (using \code{\link{hclust}}).
    }
    \item{gene_vars}{
        optional data.frame of gene-associated variables to be plotted as a colored sidebar.
        This can be the same input as for an enrichment analysis, i.e. a data.frame with genes in the first columns and with additional columns depending on the test used in \code{\link{aba_enrich}}.
        For a hypergeometric test the genes are grouped into candidate (red) and background genes (black);
        for the Wilcoxon rank-sum test, the side bar indicates the user-defined scores;
        for the binomial test the side bar shows *A/(A+B)* and for the 2x2 contingency table test *((A+1)/(B+1)) / ((C+1)/(D+1))* (+1 added to prevent division by 0).
    }
    \item{main}{
        optional character to use as a title for the heatmap.
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
## from children (age_category 3)
expr = get_expression(structure_ids=c("Allen:10657","Allen:10208"),
    gene_ids=c('RFFL', 'NTS', 'LIPE', 'GALNT6', 'GSN', 'BTBD16', 'CERS2'),
    dataset="5_stages")
expr_age_3 = expr[[3]]
plot_expression(expr_age_3)

## see the package vignette for more examples
}

