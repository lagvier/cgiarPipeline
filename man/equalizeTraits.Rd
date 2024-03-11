\name{equalizeTraits}
\alias{equalizeTraits}

\title{
Trait Equalizing
}
\description{
Perform trait equalizing when traits are named in a different way but the are in the same units.
}
\usage{
equalizeTraits(object, traits, newName=NULL)
}

\arguments{
  \item{object}{
The data object structure produced from bioflow which contains data, metadata, modifications and tables of results (metrics, predictions, modeling, status), in order to pick the necessary information to calculate a new trait.
}
\item{traits}{
A set of columns or traits that are considered to be the same but just named differently
}
\item{newName}{
A new name to add to the trait
}

}
\details{

This is a very simple function that averages all traits that supposed to be the same trait and creates a new trait. Please do not use this feature if your traits are in different units.

}
\value{
If arguments are properly specified the function returns:
\describe{
\item{result}{the same input object with an additional trait merging all suggested traits.}
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

object2 <- equalizeTraits(object=result, traits=c("Plant_Height_cm","Ear_Height_cm"), newName="newTrait")

tail(object2$data$pheno)
object2$metadata$pheno

}

