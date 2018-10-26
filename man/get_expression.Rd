\name{get_expression}
\alias{get_expression}
\title{
    Get expression data for given genes and brain structure IDs
}
\description{
    Expression data obtained from the Allen Brain Atlas project [1-4].   
}
\usage{
    get_expression(structure_ids, gene_ids, dataset = "adult")
}
\arguments{
    \item{structure_ids}{
        vector of brain structure IDs, e.g. 'Allen:10208'.
    }
    \item{gene_ids}{
        vector of gene identifiers, either Entrez-ID, Ensembl-ID or gene-symbol.
    }
    \item{dataset}{
        'adult' (default) for the microarray dataset of adult human brains; '5_stages' for RNA-seq expression data of the developing human brain, grouped into 5 developmental stages;
        'dev_effect' for a developmental effect score.
    }
}
\details{
    Get gene expression in defined brain regions from adult or developing humans, or a developmental effect score for the developing human brain.
    Expression data is obtained from the Allen Brain Atlas project [1-4], averaged across donors, and for the developing human brain divided into five major age categories.
    The developmental effect score is based on expression data of the developing human brain.
    If a requested brain region has no expression data annotated, data from sampled substructures of this region is returned.\cr 
    Please refer to the \code{ABAData} package vignette for details on the datasets.\cr
}

\value{
    A matrix with expression values or developmental effect scores per brain region (rows) and gene (columns).\cr
    For expression data from the developing human brain ('5_stages') it is a list with an expression matrix for each of the 5 developmental stages.\cr 
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
    \code{\link{plot_expression}}\cr
    \code{\link{aba_enrich}}\cr
    \code{\link{get_name}}\cr
    \code{\link{get_id}}\cr
    \code{\link{get_sampled_substructures}}\cr
    \code{\link{get_annotated_genes}}\cr
}

\examples{ 
## get expression data of six genes in two brain regions 
## from the developing human brain, 
## each of the five list elements corresponds to an age category 
get_expression(structure_ids=c('Allen:10657','Allen:10208'),
    gene_ids=c('ENSG00000168036','ENSG00000157764','ENSG00000182158',
    'ENSG00000147889'), dataset='5_stages')

## see the package vignette for more examples
}

