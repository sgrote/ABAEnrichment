\name{aba_enrich}
\alias{aba_enrich}
\title{
Test genes for expression enrichment in human brain regions
}
\description{
Tests for enrichment of user defined candidate genes in the set of expressed protein coding genes in different human brain regions. It integrates the expression of the candidate gene set (averaged across donors) and the structural information of the brain using an ontology, both provided by the Allen Brain Atlas project [1-4]. The statistical analysis is performed by interfacing the ontology enrichment software FUNC [5]. 
}
\usage{
aba_enrich(genes, dataset = 'adult', test = 'hyper', 
  cutoff_quantiles = seq(0.1, 0.9, 0.1), n_randsets = 1000, gene_len = FALSE, circ_chrom = FALSE)
}
\arguments{
  \item{genes}{
If \code{test = 'wilcoxon'} a numeric vector of scores.
If \code{test = 'hyper'} (default) a binary vector with 1 for candidate genes and 0 for background
genes. If no background genes are defined, all remaining protein coding genes are used as background.
The names of the vector are the gene identifiers: either Entrez-ID, Ensembl-ID or HGNC-symbol. 
For \code{test = 'hyper'} the names of the vector can also describe chromosomal regions ('chr:start-stop').
}
  \item{dataset}{
'adult' for the microarray dataset of adult human brains; '5_stages' for RNA-seq expression data for different stages of the developing human brain, grouped into 5 developmental stages;
'dev_effect' for a developmental effect score. For details see \code{vignette("ABAData",package="ABAData")}.
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
logical. If \code{test = 'hyper'} the probability of a background gene to be chosen as a candidate gene in a random set is dependent on the gene length. 
}
  \item{circ_chrom}{
logical. When \code{genes} defines chromosomal regions, \code{circ_chrom = TRUE} uses background regions from the same chromosome and allows randomly chosen blocks to overlap multiple background regions. Only if \code{test = 'hyper'}.
}
}
\details{
The function \code{aba_enrich} performs enrichment analyses of candidate genes within expressed protein coding genes in human brain regions. The brain regions are categorized using an ontology. Enrichment of candidate genes is tested using the hypergeometric or the Wilcoxon rank test of the ontology enrichment software FUNC [5].\cr
The hypergeometric test evaluates the enrichment of expressed candidate genes compared to a set of expressed background genes for each brain region. The background genes can be defined explicitly like the candidate genes or, as default, consist of all protein coding genes from the dataset, which are not candidate genes.\cr
To account for multiple testing the FWER is computed using random permutations of candidate and background genes (see package vignette for details). By default each gene is chosen with the same probablity as a random candidate gene. If \code{gene_len = TRUE} the probability is dependent on the gene length, i.e. a gene that is twice as long as another gene is also twice as likely to be chosen as a random candidate gene.

Instead of defining candindate and background genes explicitly in the \code{genes} input vector, it is also possible to define entire chromosomal regions as candidate and background regions. The expression enrichment is then tested for all protein coding genes located in or overlapping the candidate region on the plus or the minus strand. The gene coordinates used to identify those genes were obtained from \url{http://grch37.ensembl.org/biomart/martview/}. For the random permutations used to compute the FWER, blocks as long as candidate regions are chosen from the background regions and genes contained in these blocks are considered candidate genes. The output of \code{aba_enrich} is identical to the one that is produced for single genes.\cr
To define chromosomal regions in the input vector, the names of the 1/0 vector have to be of the form \code{chr:start-stop}, where 'start' always has to be smaller than 'stop'. Note that this option requires the input of background regions. If multiple candidate regions are provided, in the randomsets they are placed randomly, but non-overlapping into the background regions. If the background regions are relatively small, it can happen that the remaining background regions available (after a candiate region has been placed there) are too small for the next candidate region to fit entirely and non-overlapping. In this case the random selection of candidate regions inside the background regions is restarted. If this fails 10 times \code{aba_enrich} quits.\cr
An alternative method to choose random blocks from the background regions can be used with the option \code{circ_chrom=TRUE}. Every candidate region is then compared to background regions on the same chromosome. And in contrast to the default \code{circ_chrom=FALSE}, randomly chosen blocks do not have to be located inside a single background region, but are allowed to overlap multiple background regions. This means that a randomly chose block can consist of the end of the last background region and the beginning of the first background region on a given chromosome.

The Wilcoxon rank test does not compare candidate and background genes, but the user defined scores associated with the candidate genes, i.e. it compares the ranks of the scores of expressed genes in a given brain region to the ranks of all candidate genes that are expressed somewhere in the brain.\cr 
In addition to gene expression the enrichment may refer to a developmental effect score, which describes how much a gene's expression changes over time. Three different \code{datasets} can be used with \code{aba_enrich}: first, the developmental effect score, second, microarray data from adult donors and third, RNA-seq data from donors of five different developmental stages (prenatal, infant, child, adolescent, adult). In the latter case the analyses are performed independently for each developmental stage.\cr 
The expression definition for genes is variable. Different quantiles of expression over all genes are used (e.g. the lowest 40\% of gene expression are 'not expressed' and the upper 60\% are 'expressed' for a quantile of 0.4). These cutoffs are set with the parameter \code{cutoff_quantiles} and an analysis is run for every cutoff separately.
}

\value{

A list with components
\item{results}{a dataframe with the FWERs from the enrichment analyses per brain region and age category, ordered by 'age_category', 'times_FWER_under_0.05', 'mean_FWER' and 'min_FWER'; with 'min_FWER' for example denoting the minimum FWER for expression enrichment of the candidate genes in this brain region across all expression cutoffs. 'FWERs' is a semicolon separated string with the single FWERs for all cutoffs. 'equivalent_structures' is a semicolon separated string that lists structures with identical expression data due to lack of independent expression measurements in all regions.}
\item{genes}{a vector of the requested genes, excluding those genes for which no expression data is available and which therefore were not included in the enrichment analysis.}    
\item{cutoffs}{a dataframe with the expression values that correspond to the requested cutoff quantiles.}   
}
\references{
[1] Hawrylycz, M.J. et al. (2012) An anatomically comprehensive atlas of the adult human brain transcriptome, Nature 489: 391-399. \url{doi:10.1038/nature11405}\cr
[2] Miller, J.A. et al. (2014) Transcriptional landscape of the prenatal human brain, Nature 508: 199-206. \url{doi:10.1038/nature13185}\cr
[3] Allen Institute for Brain Science. Allen Human Brain Atlas [Internet]. Available from:
\url{http://human.brain-map.org/}\cr
[4] Allen Institute for Brain Science. BrainSpan Atlas of the Developing Human Brain
[Internet]. Available from: \url{http://brainspan.org/}\cr
[5] Pruefer, K. et al. (2007) FUNC: A package for detecting significant associations between gene
sets and ontological, BMC Bioinformatics 8: 41. \url{doi:10.1186/1471-2105-8-41}
}
\author{
Steffi Grote
}

\seealso{
\code{vignette("ABAEnrichment",package="ABAEnrichment")}\cr
\code{vignette("ABAData",package="ABAData")}\cr
\code{\link{get_expression}}\cr
\code{\link{plot_expression}}\cr
\code{\link{get_name}}\cr
\code{\link{get_sampled_substructures}}\cr
\code{\link{get_superstructures}}\cr

}
\examples{

#### Note that arguments 'cutoff_quantiles' and 'n_randsets' are reduced to lower computational time in the examples. Using the default values is recommended.

#### Perform gene expression enrichment analysis on 13 candidate genes in five developmental 
#### stages of the human brain using the hypergeometric test implemented in FUNC[5]   
## create input vector with candidate genes (HGNC-symbols)
genes=rep(1,13)
names(genes)=c('NCAPG', 'APOL4', 'NGFR', 'NXPH4', 'C21orf59', 'CACNG2', 'AGTR1', 'ANO1',
  'BTBD3', 'MTUS1', 'CALB1', 'GYG1', 'PAX2')
## run enrichment analysis
res=aba_enrich(genes,dataset='5_stages',cutoff_quantiles=c(0.5,0.7,0.9), n_randsets=100)
## get FWERs for enrichment of candidate genes among expressed genes
fwers=res[[1]]
## see results for the brain regions with highest enrichment for children (age_category 3)
head(fwers[fwers[,1]==3,])
## see the input genes vector (only genes with expression data available) 
res[2]
## see the expression values that correspond to the requested cutoff quantiles
res[3]

#### Perform the same analysis, but with random sets dependent on gene length
res=aba_enrich(genes,dataset='5_stages',cutoff_quantiles=c(0.5,0.7,0.9), n_randsets=100, gene_len=TRUE)

#### Perform gene expression enrichment analysis for a chromosomal region 
#### for the adult human brain using the hypergeometric test implemented in FUNC[5] 
## create input vector with a candidate regions on chromosome 3 and background regions on chromosome 3, 4 and 5
genes = c(1,rep(0,6))
names(genes) = c('3:76500000-90500000', '3:0-91600000', '3:92500000-198000000','4:3600000-50300000', '4:51400000-191100000', '5:0-47000000', '5:48600000-180700000')
## run enrichment analysis for the 'adult' dataset
res = aba_enrich(genes,dataset='adult', cutoff_quantiles=c(0.5,0.7,0.9), n_randsets=100)
## look at he results from the enrichment analysis
head(res[[1]])
## see which genes are located in the candidate regions
input_genes = res[[2]]
candidate_genes = input_genes[input_genes==1]
candidate_genes

#### Perform gene expression enrichment analysis on 15 candidate genes in the  
#### adult human brain using the wilcoxon rank test implemented in FUNC[5] 
## create input vector with random scores associated with the candidate genes (Entrez-Ids)
genes=sample(1:50,15)
names(genes)=c(324,8312,673,1029,64764,1499,3021,3417,3418,8085,3845,9968,5290,5727,5728)
## run enrichment analysis
res=aba_enrich(genes,dataset='adult',test='wilcoxon',cutoff_quantiles=c(0.2,0.5,0.8),
  n_randsets=100)
## see results for the brain regions with highest enrichment 
head(res[[1]])
## see the input genes vector and the expression values that correspond 
## to the requested cutoff quantiles
res[2:3]

}

\keyword{ htest }
