\name{ABAEnrichment-package}
\alias{ABAEnrichment-package}
\alias{ABAEnrichment}
\docType{package}
\title{
Gene expression enrichment in human brain regions
}
\description{
The package ABAEnrichment is designed to test for enrichment of user defined candidate genes in the set of expressed genes in different human brain regions. The package integrates the expression of the candidate gene set (averaged across donors) and the structural information of the brain using an ontology, both provided by the Allen Brain Atlas project [1-4]. The statistical analysis is performed by the core function \code{\link{aba_enrich}} which interfaces the ontology enrichment software FUNC [5]. Additional functions provided in this package like \code{\link{get_expression}} and \code{\link{plot_expression}} facilitate exploring the expression data.
}
\details{
\tabular{ll}{
Package: \tab ABAEnrichment\cr
Type: \tab Package\cr
Version: \tab 1.7.2\cr
Date: \tab 2017-06-25\cr
License: \tab GPL (>= 2)\cr
}

For details see \code{vignette("ABAEnrichment", package="ABAEnrichment")}\cr

}

\author{
Steffi Grote\cr
Maintainer: Steffi Grote <steffi_grote@eva.mpg.de>
}
\references{
[1] Hawrylycz, M.J. et al. (2012) An anatomically comprehensive atlas of the adult human brain transcriptome, Nature 489: 391-399. \doi{10.1038/nature11405}\cr
[2] Miller, J.A. et al. (2014) Transcriptional landscape of the prenatal human brain, Nature 508: 199-206. \doi{10.1038/nature13185}\cr
[3] Allen Institute for Brain Science. Allen Human Brain Atlas [Internet]. Available from:
\url{http://human.brain-map.org/}\cr
[4] Allen Institute for Brain Science. BrainSpan Atlas of the Developing Human Brain
[Internet]. Available from: \url{http://brainspan.org/}\cr
[5] Pruefer, K. et al. (2007) FUNC: A package for detecting significant associations between gene
sets and ontological annotations, BMC Bioinformatics 8: 41. \doi{10.1186/1471-2105-8-41}
}
\keyword{ package }
\keyword{ htest }
\seealso{
\code{vignette("ABAEnrichment", package="ABAEnrichment")}\cr
\code{vignette("ABAData", package="ABAData")}\cr
\code{\link{aba_enrich}}\cr
\code{\link{get_expression}}\cr 
\code{\link{plot_expression}}\cr
\code{\link{get_name}}\cr
\code{\link{get_id}}\cr
\code{\link{get_sampled_substructures}}\cr
\code{\link{get_superstructures}}\cr
\code{\link{get_annotated_genes}}\cr
}


