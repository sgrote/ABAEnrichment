\name{get_name}
\alias{get_name}
\title{
Get the full name of a brain region given structure ids
}
\description{
Returns the full name of brain regions given the structure IDs, e.g. 'Allen:10657' as used throughout the ABAEnrichment package. The full name is composed of an acronym and the name as used by the Allen Brain Atlas project [1-2]. 
}
\usage{
get_name(structure_ids)
}
\arguments{
  \item{structure_ids}{
a vector of brain structure IDs, e.g. c('Allen:10657','Allen:10173')
}
}

\value{
vector of the full names of the brain structures; composed of acronym, underscore and name.
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
\note{
The acronym is added because the names alone are not unique.
}


\seealso{
\code{\link{get_sampled_substructures}}\cr
\code{\link{get_superstructures}}\cr
}
\examples{
## get the full names of the brain structures 'Allen:10657' and 'Allen:10225'
get_name(c('Allen:10657','Allen:10225'))
}



