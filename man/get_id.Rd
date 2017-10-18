\name{get_id}
\alias{get_id}
\title{
Get the structure ID of a brain region given its name
}
\description{
Returns brain regions given (part of) their name, together with their structure IDs from the ontologies for the adult and for the developing brain (e.g. 'Allen:10657' as used throughout the ABAEnrichment package).

}
\usage{
get_id(structure_name)
}
\arguments{
  \item{structure_name}{
(partial) name of a brain structure, e.g. 'telencephalon' 
}
}

\value{
a data frame with the full names of the brain structures that contain \code{structure_name}; together with the ontology ('developmental' or 'adult') and the structure IDs.
}
\references{
[1] Allen Institute for Brain Science. Allen Human Brain Atlas [Internet]. Available from:
\url{http://human.brain-map.org/}\cr
[2] Allen Institute for Brain Science. BrainSpan Atlas of the Developing Human Brain
[Internet]. Available from: \url{http://brainspan.org/}
}
\author{
Steffi Grote
}


\seealso{
\code{\link{get_name}}\cr
\code{\link{get_sampled_substructures}}\cr
\code{\link{get_superstructures}}\cr
\code{browseVignettes("ABAData")}\cr

}
\examples{
## get structure IDs of brain regions that contain 'accumbens' in their names
get_id('accumbens')
## get structure IDs of brain regions that contain 'telencephalon' in their names
get_id('telencephalon')
## get all brain regions that have direct or indirect expression data together with thier IDs
all_regions = get_id('')
head(all_regions)
}



