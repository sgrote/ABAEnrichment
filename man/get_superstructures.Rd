\name{get_superstructures}
\alias{get_superstructures}
\title{
    Returns all superstructures of a brain region using the Allen Brain Atlas ontology
}
\description{
    Returns all superstructures of a brain region and the brain region itself given a structure ID, e.g. 'Allen:10657' as used throughout the ABAEnrichment package.
    The output vector contains the superstructures according to the hierarchy provided by the Allen Brain Atlas ontology [1,2] beginning with the root ('brain' or 'neural plate') and ending with the requested brain region. 
}

\usage{
    get_superstructures(structure_id)
}

\arguments{
    \item{structure_id}{
        a brain structure ID, e.g. 'Allen:10657' or '10657'
    }
}

\value{
    vector of brain structure IDs that contains all superstructures of the requested brain region and the brain region itself.
    The order of the brain regions follows the hierarchical organization of the brain.
}

\references{
    [1] Allen Institute for Brain Science. Allen Human Brain Atlas.
    Available from: \url{http://human.brain-map.org/} \cr
    [2] Allen Institute for Brain Science. BrainSpan Atlas of the Developing Human Brain.
    Available from: \url{http://brainspan.org/} 
}
\author{
    Steffi Grote
}
\note{
    The ontologies for the adult and the developing human brain are different.
}


\seealso{
    \code{\link{get_name}}\cr
    \code{\link{get_id}}\cr
    \code{\link{get_sampled_substructures}}\cr
}

\examples{
## Get the IDs of the superstructures of the precentral gyrus
## (adult brain ontology)
get_superstructures('Allen:4010')

## Get the IDs and the names of the superstructures
## of the dorsolateral prefrontal cortex 
## (developing brain ontology) 
data.frame(superstructure=get_name(get_superstructures("Allen:10173")))
}

