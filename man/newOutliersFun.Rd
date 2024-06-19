\name{newOutliersFun}
\alias{newOutliersFun}

\title{
Outlier identification
}
\description{
Perform the outlier identification using threshold values for the interquantile ranges.
}
\usage{
newOutliersFun(myObject, trait, outlierCoefOutqPheno, 
               traitLBOutqPheno = NULL, traitUBOutqPheno = NULL)
}

\arguments{
  \item{myObject}{
The data object structure produced from bioflow which contains data, metadata, modifications and tables of results (metrics, predictions, modeling, status), in order to pick the necessary information to calculate an index.
}
\item{trait}{
A character vector with the traits that should be scanned for outliers given the next parameters.
}
\item{outlierCoefOutqPheno}{
this determines how far the plot whiskers extend out from the box. If coef is positive, the whiskers extend to the most extreme data point which is no more than coef times the length of the box away from the box. A value of zero causes the whiskers to extend to the data extremes (and no outliers be returned).
}
\item{traitLBOutqPheno}{
Lower bound threshold determining as an outlier any value smaller than this.
}
\item{traitUBOutqPheno}{
Upper bound threshold determining as an outlier any value greater than this.
}

}
\details{

The first step in genetic evaluation is to ensure that input phenotypic records are of good quality. This option aims to allow users to select outliers based on plot whiskers and absolute values.

}
\value{
If arguments are properly specified the function returns:
\describe{
\item{outleirs}{a table of modifications for the phenotypic dataset givent the requested parameters.}
}
}


\examples{
####=========================================####
data(DT_example)

result$status
result$metadata$pheno

final <- newOutliersFun(myObject=result, trait=c("Plant_Height_cm"), 
                        outlierCoefOutqPheno=2.5, traitLBOutqPheno=0, 
                        traitUBOutqPheno=Inf)

tail(final)

}

