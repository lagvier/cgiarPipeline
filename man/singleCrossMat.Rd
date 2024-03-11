\name{singleCrossMat}
\alias{singleCrossMat}

\title{
Single-Cross Marker Matrix Built.
}
\description{
Calculates the marker matrix for single crosses given that markers for parental lines are provided.
}
\usage{
singleCrossMat( 
  object= NULL,
  hybridBatch=1000,
  allHybrids=FALSE,
  verbose=TRUE
)
}

\arguments{
  \item{object}{
The data object structure produced from bioflow which contains data, metadata, modifications and tables of results (metrics, predictions, modeling, status).
}
\item{hybridBatch}{
The number of hybrids to build in each batch. Building all hybrids at once can be computationally intensive and unnecesary. The default value is 1000 hybrids per batch.
}
\item{allHybrids}{
an indication to know if only the marker-profiles for the hybrids present in the phenotypic dataset should be computed or marker-profiles for all cross combinations should be created. This should be used carefully when the number of males and females in the pedigree file is big..
}
\item{verbose}{
Should the logs or warninig messages of the analysis be printed?
}

}
\details{

Hybrid breeding based on single-cross of lines allows the genotyping of only the parental lines and the consequent formation of single-cross marker matrices. The idea is that the user provides the marker information from the parental lines and has uploaded a pedigree of the phenotypic dataset. The marker profiles of the possible cross combinations will be computed based on the availability of the marker information..

}
\value{
If arguments are properly specified the function returns:
\describe{
\item{result}{the same input object with additional marker matrix for the hybrids.}
}
}

\references{

Nishio M and Satoh M. 2014. Including Dominance Effects in the Genomic BLUP Method for Genomic Evaluation. Plos One 9(1), doi:10.1371/journal.pone.0085792

Su G, Christensen OF, Ostersen T, Henryon M, Lund MS. 2012. Estimating Additive and Non-Additive Genetic Variances and Predicting Genetic Merits Using Genome-Wide Dense Single Nucleotide Polymorphism Markers. PLoS ONE 7(9): e45293. doi:10.1371/journal.pone.0045293

}

\examples{
####=========================================####

data(DT_example)

result$status
result$metadata$pheno
# mock up the genotype data from parents
M <- result$data$geno
parents <- unique(c(result$data$pedigree$Hybrid_Parent1, result$data$pedigree$Hybrid_Parent2))
Mparents <- M[1:length(parents),]
rownames(Mparents) <- parents
result$data$geno <- Mparents
dim(result$data$geno)

# compute the hybrid genotypes
final <- singleCrossMat( # single cross matrix function
  object= result,
  hybridBatch=100,
  allHybrids=FALSE,
  verbose=TRUE
)

dim(final$data$geno)
intersect(rownames(final$data$geno), parents) # parents available before
setdiff(rownames(final$data$geno), parents) # hybrids available now

}

