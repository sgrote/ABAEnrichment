\name{get_sampled_substructures}
\alias{get_sampled_substructures}
\title{
Return sampled substructures of a given brain region
}
\description{
The function returns for a given brain structure ID all its substructures with available expression data, potentially including the structure itself.
}
\usage{
get_sampled_substructures(structure_id)
}
\arguments{
  \item{structure_id}{
a brain structure ID, e.g. 'Allen:10657' or '10657'
}
}
\details{
The ontology enrichment analysis in \code{\link{aba_enrich}} tests all brain regions for which data is available, although the region might not have been sampled directly. In this case the region inherits the expression data from its substructures with available expression data. The function \code{get_sampled_substructures} helps to explore where the expression data for a brain region came from.

}
\value{
vector of brain structure IDs that contains all substructures of the requested brain region that were sampled.
}
\references{
[1] Allen Institute for Brain Science. Allen Human Brain Atlas [Internet]. Available from:
\url{http://human.brain-map.org/} \cr
[2] Allen Institute for Brain Science. BrainSpan Atlas of the Developing Human Brain
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
\code{\link{get_superstructures}}\cr
}
\examples{
## get the brain structures from which the brain structures 'Allen:4010' and 'Allen:10208' 
## inherit their expression data
get_sampled_substructures('Allen:4010')
get_sampled_substructures('Allen:10208')

}

