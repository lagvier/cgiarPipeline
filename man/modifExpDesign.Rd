\name{modifExpDesign}
\alias{modifExpDesign}

\title{
Modify Experimental Design
}
\description{
Perform modifications to the experimental design based on a data frame which specifies which factors should be silenced by adding them to the modifications table of the data object.
}
\usage{
modifExpDesign(object, df)
}

\arguments{
  \item{object}{
The data object structure produced from bioflow which contains data, metadata, modifications and tables of results (metrics, predictions, modeling, status), in order to use the data and identify the records to be silenced for a given column.
}
\item{df}{
A table of rows specifying environments and columns specifying experimental design factors in order to specify with a zero value the factors that should be silenced for a given environment. The records matching will be added to the modifications table for the phenotype dataset to be silenced during the analyticial modules if the module is enabled to do so.
}

}
\details{

This is a very simple function that identifies which records match the specified table and keep track of those to be later silenced if requested.

}
\value{
If arguments are properly specified the function returns:
\describe{
\item{result}{the same input object with an additional modifications.}
}
}

\references{

None

}

\examples{
####=========================================####
data(DT_example)

# create a table of environments by factors
dtProv <- result$data$pheno
dtProv$environment <- dtProv$Env
head(dtProv)
fieldNames <- as.character(unique(dtProv[,"environment"]))
spD <- split(dtProv,dtProv[,"environment"])
presentFactors <- c("Replicate", "Block","Range", "Pass")
presentFactorsPerField <- lapply(spD, function(x){
  apply(x[,presentFactors, drop=FALSE], 2, function(y){length(unique(y))})
})
presentFactorsPerField <- do.call(rbind, presentFactorsPerField)
dtProvTable = as.data.frame(presentFactorsPerField);  rownames(dtProvTable) <- fieldNames
dtProvTable
# set to zero one of the factors in one environment
dtProvTable[2,3]=0
dtProvTable
# apply the function to create the modifications table
xx <- modifExpDesign(object=result, df=dtProvTable)
tail(xx$modifications$pheno)
xx$status

}

