\name{traitTransformation}
\alias{traitTransformation}

\title{
Trait Transformation
}
\description{
Perform trait transformation for some available functions.
}
\usage{
traitTransformation(
    object= NULL,
    trait=NULL,
    transformation = NULL,
    verbose=FALSE
)
}

\arguments{
  \item{object}{
The data object structure produced from bioflow which contains data, metadata, modifications and tables of results (metrics, predictions, modeling, status).
}
\item{trait}{
A set of columns or traits that are considered to be the transformed
}
\item{transformation}{
Transformation to be transformed.
}
\item{verbose}{
Should the logs or warninig messages of the analysis be printed?
}

}
\details{

This is a very simple function that transforms traits.

}
\value{
If arguments are properly specified the function returns:
\describe{
\item{result}{the same input object with an additional traits.}
}
}

\references{

None

}

\examples{
####=========================================####

data(DT_example)

result$status
result$metadata$pheno

final <- traitTransformation(
    object= result,
    trait="Plant_Height_cm", # per trait
    transformation = c("log"),
    verbose=FALSE
)
  
tail(final$data$pheno)

}

